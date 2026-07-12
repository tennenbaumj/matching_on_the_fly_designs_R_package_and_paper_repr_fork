#' Negative Binomial Regression Inference for Count Responses
#'
#' Fits a negative binomial regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountNegBin$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountNegBin = R6::R6Class("InferenceCountNegBin",
	lock_objects = FALSE,
	inherit = InferenceCountLikelihood,
	public = list(

		#' @description Initialize a negative binomial regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose               Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart optimizer start values by default.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			X_data = private$get_X()
			X_full = if (is.null(X_data) || ncol(X_data) == 0) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					df = as.data.frame(X_fit[, -1, drop = FALSE])
					df$y = as.numeric(private$y)
					ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					fit = tryCatch(
						suppressWarnings(
							MASS::glm.nb(
								y ~ .,
								data = df,
								weights = row_weights,
								start = if (!is.null(ws_args$start_params) && length(ws_args$start_params) >= ncol(X_fit)) ws_args$start_params[1:ncol(X_fit)] else NULL,
								init.theta = if (!is.null(ws_args$start_params) && length(ws_args$start_params) >= (ncol(X_fit) + 1L)) exp(ws_args$start_params[ncol(X_fit) + 1L]) else NULL
							)
						),
						error = function(e) NULL
					)
					if (is.null(fit)) return(NULL)
					coef_vec = stats::coef(fit)[colnames(X_fit)]
					if (length(coef_vec) != ncol(X_fit) || any(!is.finite(coef_vec))) return(NULL)
					list(
						b = as.numeric(coef_vec),
						theta_hat = fit$theta,
						fisher_information = tryCatch(solve(stats::vcov(fit)), error = function(e) NULL),
						ssq_b_j = NA_real_
					)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && length(mod$b) >= 2L && is.finite(mod$b[2L])
				}
			)
			private$cached_mod = attempt$fit
			if (is.null(attempt$fit) || is.null(attempt$fit$b) || length(attempt$fit$b) < 2L || !is.finite(attempt$fit$b[2L])) {
				fallback = tryCatch(
					{
						df_fb = as.data.frame(X_full[, -1, drop = FALSE])
						df_fb$y = as.numeric(private$y)
						fit_fb = stats::glm(
							y ~ .,
							data = df_fb,
							family = stats::poisson(link = "log"),
							weights = row_weights
						)
						list(
							b = as.numeric(stats::coef(fit_fb)[colnames(X_full)]),
							fisher_information = tryCatch(solve(stats::vcov(fit_fb)), error = function(e) NULL)
						)
					},
					error = function(e) NULL
				)
				if (is.null(fallback) || is.null(fallback$b) || length(fallback$b) < 2L || !is.finite(fallback$b[2L])) {
					private$cached_values$beta_hat_T = NA_real_
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$df = NA_real_
					return(NA_real_)
				}
				private$cached_mod = fallback
				private$clear_nonestimable_state()
				private$cached_values$beta_hat_T = as.numeric(fallback$b[2L])
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				private$set_fit_warm_start(as.numeric(fallback$b), "beta", fisher = fallback$fisher_information, force_pd = TRUE)
				return(private$cached_values$beta_hat_T)
			}
			private$clear_nonestimable_state()
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[2L])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$set_fit_warm_start(
				c(as.numeric(attempt$fit$b), log(as.numeric(attempt$fit$theta_hat))),
				"params",
				fisher = attempt$fit$fisher_information,
				force_pd = TRUE
			)
			private$cached_values$beta_hat_T
		},
		#' @description Negative-binomial delete-one refits are unstable for
		#'   jackknife inference; report explicit non-estimability.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_estimate = function(unit = "auto"){
			private$cache_nonestimable_se("negbin_jackknife_not_supported")
			NA_real_
		},
		#' @description Non-estimable jackknife bias-corrected estimate.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_corrected_estimate = function(unit = "auto"){
			self$compute_jackknife_estimate(unit = unit)
		},
		#' @description Non-estimable jackknife bias estimate.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_bias_estimate = function(unit = "auto"){
			private$cache_nonestimable_se("negbin_jackknife_not_supported")
			NA_real_
		},
		#' @description Non-estimable jackknife standard error.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_std_error = function(unit = "auto"){
			private$cache_nonestimable_se("negbin_jackknife_not_supported")
			NA_real_
		},
		#' @description Alias for \code{compute_jackknife_std_error()}.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_standard_error = function(unit = "auto"){
			self$compute_jackknife_std_error(unit = unit)
		},
		#' @description Non-estimable jackknife Wald two-sided p-value.
		#' @param delta Null treatment-effect value. Default 0.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_wald_two_sided_pval = function(delta = 0, unit = "auto"){
			private$cache_nonestimable_se("negbin_jackknife_not_supported")
			NA_real_
		},
		#' @description Non-estimable jackknife Wald confidence interval.
		#' @param alpha Significance level. Default 0.05.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			private$cache_nonestimable_se("negbin_jackknife_not_supported")
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		best_X_colnames = NULL,
		negbin_X_full_cache = NULL,
		negbin_w_cache = NULL,
		get_complexity_tier = function() "heavy",
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_X_colnames
			X_data = private$get_X()
			if (length(X_cols) == 0L){
				X = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(1, treatment = private$w, X_cov)
			}
			ws_args = private$get_backend_warm_start_args(ncol(X) + 1L)
			res = tryCatch(
				fast_neg_bin_cpp(
					X = X, y = as.integer(private$y),
					warm_start_params = ws_args$start_params,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(c(as.numeric(res$b), log(as.numeric(res$theta_hat))), "params", fisher = res$fisher_information, force_pd = TRUE)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_lik_ratio_param_bootstrap = function(){
			FALSE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		get_supported_testing_types_impl = function(){
			c("wald", "score", "gradient")
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			stop(class(self)[1], " does not support likelihood-ratio p-values.", call. = FALSE)
		},
		compute_lik_ratio_confidence_interval_impl = function(alpha){
			stop(class(self)[1], " does not support likelihood-ratio confidence intervals.", call. = FALSE)
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null   = as.numeric(null_fit$b)
			theta    = as.numeric(null_fit$theta_hat)
			if (!is.finite(theta) || theta <= 0) return(NULL)
			mu       = pmax(exp(as.numeric(spec$X %*% b_null)), 0)
			y_sim    = as.integer(rnbinom(length(mu), size = theta, mu = mu))
			X_fit    = spec$X
			j        = spec$j
			full_res = tryCatch(
				fast_neg_bin_cpp(
					X = X_fit, y = y_sim,
					smart_cold_start = private$smart_cold_start_default,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || length(full_res$b) < j || !is.finite(full_res$b[j])) return(NULL)
			full_fit_boot = list(
				b          = as.numeric(full_res$b),
				theta_hat  = full_res$theta_hat,
				neg_loglik = -as.numeric(full_res$logLik)
			)
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = tryCatch(
						fast_neg_bin_cpp(
							X = X_fit, y = y_sim,
							warm_start_params = start %||% c(as.numeric(full_res$b), log(as.numeric(full_res$theta_hat))),
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					res
				},
				neg_loglik = function(fit){
					val = fit$neg_loglik %||% fit$neg_log_lik %||% fit$neg_ll
					if (!is.null(val)) return(as.numeric(val)[1L])
					loglik = fit$logLik %||% fit$log_lik
					if (!is.null(loglik)) return(-as.numeric(loglik)[1L])
					NA_real_
				}
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					fast_neg_bin_cpp(
						X = X_fit,
						y = as.integer(y),
						warm_start_params = start %||% ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_cold_start = private$smart_cold_start_default,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){
					c(as.numeric(fit$b), log(as.numeric(fit$theta_hat)))
				},
				score = function(fit){
					get_negbin_regression_score_cpp(X_fit, y, c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))
				},
				observed_information = function(fit){
					-get_negbin_regression_hessian_cpp(X_fit, y, c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))
				},
				fisher_information = function(fit){
					-get_negbin_regression_hessian_cpp(X_fit, y, c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))
				},
				information = function(fit){
					-get_negbin_regression_hessian_cpp(X_fit, y, c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))
				},
				neg_loglik = function(fit){
					val = fit$neg_loglik %||% fit$neg_log_lik %||% fit$neg_ll
					if (!is.null(val)) return(as.numeric(val)[1L])
					loglik = fit$logLik %||% fit$log_lik
					if (!is.null(loglik)) return(-as.numeric(loglik)[1L])
					NA_real_
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			if (is.null(private$negbin_X_full_cache) || !identical(private$w, private$negbin_w_cache)) {
				X_data = private$get_X()
				private$negbin_X_full_cache = if (is.null(X_data) || ncol(X_data) == 0) {
					cbind(`(Intercept)` = 1, treatment = private$w)
				} else {
					cbind(`(Intercept)` = 1, treatment = private$w, X_data)
				}
				private$negbin_w_cache = private$w
			}
			X_full = private$negbin_X_full_cache
			
			if (!private$harden) {
				ws_args = private$get_backend_warm_start_args(ncol(X_full) + 1L)
				if (estimate_only) {
					res = tryCatch(
						fast_neg_bin_cpp(
							X = X_full, y = as.integer(private$y),
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					res$beta_hat_T = as.numeric(res$b[2L])
					res$ssq_b_2 = NA_real_
					res$params = c(as.numeric(res$b), log(as.numeric(res$theta_hat)))
					res$neg_log_lik = -as.numeric(res$logLik)
					res$fisher_information = res$fisher_information
					res$j_treat = 2L
				} else {
					res = tryCatch(
						fast_neg_bin_with_var_cpp(
							X = X_full, y = as.integer(private$y),
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					res$j_treat = 2L
					res$beta_hat_T = as.numeric(res$b[2L])
					hess = res$hess_fisher_info_matrix
					vcov = tryCatch(solve(hess), error = function(e) matrix(NA_real_, nrow(hess), ncol(hess)))
					res$ssq_b_2 = if (res$j_treat <= ncol(X_full)) as.numeric(vcov[res$j_treat, res$j_treat]) else NA_real_
					res$params = c(as.numeric(res$b), log(as.numeric(res$theta_hat)))
					res$neg_log_lik = -as.numeric(res$logLik)
					res$fisher_information = res$hess_fisher_info_matrix
				}
				private$best_X_colnames = setdiff(colnames(X_full), c("(Intercept)", "treatment"))
				private$clear_nonestimable_state()
				private$cached_values$likelihood_test_context = list(X = X_full, j_treat = 2L)
				return(res)
			}

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					if (estimate_only) {
						res = tryCatch(
							fast_neg_bin_cpp(
								X = X_fit, y = as.integer(private$y),
								warm_start_params = ws_args$start_params,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								optimization_alg = private$optimization_alg
							),
							error = function(e) NULL
						)
						if (is.null(res) || !isTRUE(res$converged)) return(NULL)
						list(b = as.numeric(res$b), ssq_b_2 = NA_real_, j_treat = j_treat,
						     theta_hat = res$theta_hat, neg_loglik = -as.numeric(res$logLik),
						     fisher_information = res$fisher_information)
					} else {
						res = tryCatch(
							fast_neg_bin_with_var_cpp(
								X = X_fit, y = as.integer(private$y),
								warm_start_params = ws_args$start_params,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								optimization_alg = private$optimization_alg
							),
							error = function(e) NULL
						)
						if (is.null(res) || !isTRUE(res$converged)) return(NULL)
						hess = res$hess_fisher_info_matrix
						vcov = tryCatch(solve(hess), error = function(e) matrix(NA_real_, nrow(hess), ncol(hess)))
						ssq_b_j = if (j_treat <= ncol(X_fit)) as.numeric(vcov[j_treat, j_treat]) else NA_real_
						list(b = as.numeric(res$b), ssq_b_j = ssq_b_j, j_treat = j_treat,
						     theta_hat = res$theta_hat, neg_loglik = -as.numeric(res$logLik),
						     fisher_information = res$hess_fisher_info_matrix)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j %||% mod$ssq_b_2) && (mod$ssq_b_j %||% mod$ssq_b_2) > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$clear_nonestimable_state()
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = attempt$fit$j_treat
				)
				# Important: pack parameters for the base class shared logic
				attempt$fit$params = c(as.numeric(attempt$fit$b), log(as.numeric(attempt$fit$theta_hat)))
				attempt$fit$ssq_b_2 = attempt$fit$ssq_b_2 %||% attempt$fit$ssq_b_j
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		}
	)
)
