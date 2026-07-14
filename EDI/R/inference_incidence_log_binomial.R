#' Log-Binomial Regression Inference for Incidence Responses
#'
#' Fits a log-binomial regression for binary (incidence) responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidLogBinomial$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidLogBinomial = R6::R6Class("InferenceIncidLogBinomial",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(

		#' @description Initialize a log-binomial regression inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose               Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		#' @param harden                Whether to apply robustness measures.
		#' @param max_abs_reasonable_coef Cap for reasonable log-binomial coefficients.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL, harden = TRUE, max_abs_reasonable_coef = 25){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, harden = harden, smart_cold_start_default = smart_cold_start_default)
			private$max_abs_reasonable_coef = max_abs_reasonable_coef
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = match(2L, keep)
					res = tryCatch(
						fast_log_binomial_regression_weighted_cpp(
							X = X_fit,
							y = as.numeric(private$y),
							weights = as.numeric(row_weights),
							warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					res$j_treat = j_treat
					ssq_b_j = NA_real_
					if (!estimate_only && !is.null(res$fisher_information) &&
					    is.matrix(res$fisher_information) && nrow(res$fisher_information) >= j_treat) {
						inv_fi = tryCatch(solve(res$fisher_information), error = function(e) NULL)
						if (!is.null(inv_fi) && is.finite(inv_fi[j_treat, j_treat]) && inv_fi[j_treat, j_treat] > 0) {
							ssq_b_j = inv_fi[j_treat, j_treat]
						}
					}
					res$ssq_b_j = ssq_b_j
					res$ssq_b_2 = ssq_b_j
					res
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = match(2L, keep)
					if (!isTRUE(private$is_log_binomial_fit_reasonable(mod, X_fit, j_treat))) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j %||% mod$ssq_b_2)
				}
			)
			private$cached_mod = attempt$fit
			j_treat = match(2L, attempt$keep)
			if (!isTRUE(private$is_log_binomial_fit_reasonable(attempt$fit, attempt$X, j_treat))){
				private$cache_nonestimable_estimate("log_binomial_weighted_fit_unavailable")
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$set_fit_warm_start(attempt$fit$b, "beta", fisher = attempt$fit$fisher_information, force_pd = TRUE)
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[j_treat])
			ssq = attempt$fit$ssq_b_j %||% attempt$fit$ssq_b_2
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$beta_hat_T
		},
		compute_score_confidence_interval = function(alpha = 0.05){
			ci = tryCatch(
				super$compute_score_confidence_interval(alpha = alpha),
				error = function(e){
					msg = if (length(e$message) == 0L) "" else e$message
					if (grepl("'names' attribute", msg, fixed = TRUE) ||
					    grepl("must be the same length as the vector", msg, fixed = TRUE)) {
						private$cache_nonestimable_se("score_confidence_interval_unavailable")
						return(c(NA_real_, NA_real_))
					}
					stop(e)
				}
			)
			if (length(ci) < 2L || !all(is.finite(ci[1:2]))) {
				private$cache_nonestimable_se("score_confidence_interval_unavailable")
				return(c(NA_real_, NA_real_))
			}
			ci
		},
		compute_gradient_confidence_interval = function(alpha = 0.05){
			ci = tryCatch(
				super$compute_gradient_confidence_interval(alpha = alpha),
				error = function(e){
					msg = if (length(e$message) == 0L) "" else e$message
					if (grepl("'names' attribute", msg, fixed = TRUE) ||
					    grepl("must be the same length as the vector", msg, fixed = TRUE)) {
						private$cache_nonestimable_se("gradient_confidence_interval_unavailable")
						return(c(NA_real_, NA_real_))
					}
					stop(e)
				}
			)
			if (length(ci) < 2L || !all(is.finite(ci[1:2]))) {
				private$cache_nonestimable_se("gradient_confidence_interval_unavailable")
				return(c(NA_real_, NA_real_))
			}
			ci
		}
	),
	private = list(
		best_X_colnames = NULL,
		logbin_X_full_cache = NULL,
		logbin_w_cache = NULL,
		max_abs_reasonable_coef = 25,
		is_log_binomial_fit_reasonable = function(mod, X_fit = NULL, j_treat = 2L){
			if (is.null(mod) || is.null(mod$b)) return(FALSE)
			j_treat = as.integer(j_treat %||% mod$j_treat %||% 2L)
			if (length(j_treat) != 1L || !is.finite(j_treat) || j_treat < 1L) return(FALSE)
			b = as.numeric(mod$b)
			if (length(b) < j_treat || any(!is.finite(b))) return(FALSE)
			if (any(abs(b) > private$max_abs_reasonable_coef)) return(FALSE)
			if (!is.null(mod$converged) && !isTRUE(mod$converged)) return(FALSE)
			if (!is.null(X_fit)) {
				eta = tryCatch(as.numeric(as.matrix(X_fit) %*% b), error = function(e) NA_real_)
				if (any(!is.finite(eta))) return(FALSE)
				if (any(eta > 1e-6)) return(FALSE)
			} else if (!is.null(mod$mu_hat)) {
				mu = as.numeric(mod$mu_hat)
				if (any(!is.finite(mu)) || any(mu < -1e-10) || any(mu > 1 + 1e-10)) return(FALSE)
			}
			TRUE
		},
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
			ws_args = private$get_backend_warm_start_args(ncol(X))
			res = tryCatch(
				fast_log_binomial_regression_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)

			if (!isTRUE(private$is_log_binomial_fit_reasonable(res, X, 2L))){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$fisher_information, force_pd = TRUE)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_lik_ratio_param_bootstrap = function(){
			TRUE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		compute_gradient_confidence_interval_impl = function(alpha){
			ci = private$invert_gradient_ci_uniroot(alpha)
			if (length(ci) >= 2L && all(is.finite(ci[1:2]))) return(ci)

			# Log-binomial fixed-effect profiles can be infeasible on one side of
			# the treatment coefficient.  Use the finite likelihood-ratio profile
			# interval rather than returning an NA gradient CI.
			ci_fallback = tryCatch(private$invert_lik_ratio_ci_newton(alpha), error = function(e) c(NA_real_, NA_real_))
			if (length(ci_fallback) >= 2L && all(is.finite(ci_fallback[1:2]))) {
				names(ci_fallback) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci_fallback)
			}
			ci
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null     = as.numeric(null_fit$b)
			mu         = pmin(pmax(exp(as.numeric(spec$X %*% b_null)), 0), 1)
			y_sim      = as.numeric(rbinom(length(mu), 1L, mu))
			X_fit      = spec$X
			j          = spec$j

			# Parametric bootstrap: use observed fit as anchor
			ws_args = private$get_backend_warm_start_args(ncol(X_fit))
			full_fit_b = tryCatch(
				fast_log_binomial_regression_cpp(
					X = X_fit, y = y_sim,
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (!isTRUE(private$is_log_binomial_fit_reasonable(full_fit_b, X_fit, j))) return(NULL)
			list(
				full_fit = full_fit_b,
				fit_null = function(d, start = NULL){
					ws_args_null = private$get_backend_warm_start_args(ncol(X_fit))
					res = tryCatch(
						fast_log_binomial_regression_cpp(
							X = X_fit, y = y_sim,
							warm_start_beta = start %||% full_fit_b$b,
							warm_start_fisher_info = ws_args_null$warm_start_fisher_info,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE
						),
						error = function(e) NULL
					)
					if (!isTRUE(private$is_log_binomial_fit_reasonable(res, X_fit, j))) return(NULL)
					res
				},
				neg_loglik = function(fit){
					eta_f  = as.numeric(X_fit %*% as.numeric(fit$b))
					mu_fit = exp(eta_f)
					-sum(y_sim * log(pmax(mu_fit, 1e-15)) + (1 - y_sim) * log(pmax(1 - mu_fit, 1e-15)))
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
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					res = tryCatch(
						fast_log_binomial_regression_cpp(
							X_fit, y,
							warm_start_beta = start %||% ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							fixed_idx = j_treat, fixed_values = delta,
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)

					if (!isTRUE(private$is_log_binomial_fit_reasonable(res, X_fit, j_treat))) return(NULL)
					res
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_log_binomial_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_log_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_log_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				information = function(fit){
					-get_log_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					mu = exp(eta)
					-sum(y * log(pmax(mu, 1e-15)) + (1 - y) * log(pmax(1 - mu, 1e-15)))
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			if (is.null(private$logbin_X_full_cache) || !identical(private$w, private$logbin_w_cache)) {
				X_data = private$get_X()
				private$logbin_X_full_cache = if (is.null(X_data) || ncol(X_data) == 0) {
					cbind(`(Intercept)` = 1, treatment = private$w)
				} else {
					cbind(`(Intercept)` = 1, treatment = private$w, X_data)
				}
				private$logbin_w_cache = private$w
			}
			X_full = private$logbin_X_full_cache
			
			if (!private$harden) {
				ws_args = private$get_backend_warm_start_args(ncol(X_full))
				if (estimate_only) {
					res = tryCatch(
						fast_log_binomial_regression_cpp(
							X_full, private$y,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = TRUE
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					res$beta_hat_T = as.numeric(res$b[2L])
					res$ssq_b_j = NA_real_
					res$ssq_b_2 = NA_real_
				} else {
					res = tryCatch(
						fast_log_binomial_regression_with_var_cpp(
							X_full, private$y, j = 2L,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					res$j_treat = 2L
					res$beta_hat_T = as.numeric(res$b[2L])
					res$ssq_b_2 = res$ssq_b_j
				}
				if (!isTRUE(private$is_log_binomial_fit_reasonable(res, X_full, 2L))) {
					private$cache_nonestimable_estimate("log_binomial_fit_unavailable")
					private$cached_values$likelihood_test_context = NULL
					return(NULL)
				}
				private$best_X_colnames = setdiff(colnames(X_full), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = X_full,
					j_treat = 2L,
					full_neg_loglik = res$neg_ll
				)
				return(res)
			}

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L, # intercept and treatment
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					if (estimate_only) {
						res = tryCatch(
							fast_log_binomial_regression_cpp(
								X = X_fit, y = private$y,
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							),
							error = function(e) NULL
						)
						if (is.null(res)) return(NULL)
						list(b = res$b, beta_hat_T = as.numeric(res$b[j_treat]), ssq_b_j = NA_real_, j_treat = j_treat, fisher_information = res$fisher_information, neg_ll = res$neg_ll, converged = res$converged, mu_hat = res$mu_hat)
					} else {
						res = tryCatch(
							fast_log_binomial_regression_with_var_cpp(
								X = X_fit, y = private$y, j = j_treat,
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							),
							error = function(e) NULL
						)
						if (is.null(res)) return(NULL)
						res$j_treat = j_treat
						res$beta_hat_T = as.numeric(res$b[j_treat])
						res$ssq_b_2 = res$ssq_b_j
						res
					}
				},

				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (!isTRUE(private$is_log_binomial_fit_reasonable(mod, X_fit, j_treat))) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j %||% mod$ssq_b_2)
				}
			)
			if (!isTRUE(private$is_log_binomial_fit_reasonable(attempt$fit, attempt$X, attempt$fit$j_treat %||% match(2L, attempt$keep)))) {
				private$cache_nonestimable_estimate("log_binomial_fit_unavailable")
				private$cached_values$likelihood_test_context = NULL
				return(NULL)
			}
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = attempt$fit$j_treat,
					full_neg_loglik = attempt$fit$neg_ll
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},
		build_design_matrix = function(){
			X_data = private$get_X()
			if (is.null(X_data) || ncol(X_data) == 0) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
		}
	)
)
