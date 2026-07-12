#' Beta Regression Inference for Proportion Responses
#'
#' Fits a beta regression for proportion responses (constrained to (0, 1)) using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropBetaRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropBetaRegr = R6::R6Class("InferencePropBetaRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a beta-regression inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		#' @param optimization_alg Character scalar specifying the optimization algorithm. 
		#'   Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
				assertFormula(model_formula, null.ok = TRUE)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					y_fit = sanitize_beta_response(as.numeric(private$y))
					res = tryCatch({
						n_params = ncol(X_fit) + 1L
						ws_args = private$get_backend_warm_start_args(n_params)
						start_params = ws_args$start_params
						start_beta = if (!is.null(start_params) && length(start_params) >= ncol(X_fit)) {
							start_params[seq_len(ncol(X_fit))]
						} else {
							ws_args$start_beta
						}
						start_phi = if (!is.null(start_params) && length(start_params) >= n_params) {
							exp(start_params[n_params])
						} else {
							10
						}
						fast_beta_regression_weighted_cpp(
							X_sexp = X_fit,
							y_sexp = y_fit,
							weights_sexp = row_weights,
							warm_start_beta = start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							start_phi = start_phi,
							estimate_only = estimate_only,
							optimization_alg = private$optimization_alg
						)
					}, error = function(e) NULL)
					if (!is.null(res)) {
						coef_vec = as.numeric(res$coefficients)
						if (all(is.finite(coef_vec))) {
							ssq_b_2 = NA_real_
							if (!estimate_only && !is.null(res$fisher_information) &&
							    is.matrix(res$fisher_information) && nrow(res$fisher_information) >= 2L) {
								inv_fi = tryCatch(solve(res$fisher_information), error = function(e) NULL)
								if (!is.null(inv_fi) && is.finite(inv_fi[2L, 2L]) && inv_fi[2L, 2L] > 0) {
									ssq_b_2 = inv_fi[2L, 2L]
								}
							}
							return(list(
								b = coef_vec,
								phi = as.numeric(res$phi),
								fisher_information = res$fisher_information,
								ssq_b_2 = ssq_b_2
							))
						}
					}
					lm_fit = tryCatch(
						stats::lm.wfit(x = X_fit, y = logit(y_fit), w = row_weights),
						error = function(e) NULL
					)
					if (is.null(lm_fit) || length(lm_fit$coefficients) < 2L) return(NULL)
					list(
						b = as.numeric(lm_fit$coefficients),
						phi = NA_real_,
						fisher_information = NULL,
						ssq_b_2 = NA_real_
					)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && length(mod$b) >= 2L && is.finite(mod$b[2L])
				}
			)
			private$cached_mod = attempt$fit
			if (is.null(attempt$fit) || is.null(attempt$fit$b) || length(attempt$fit$b) < 2L || !is.finite(attempt$fit$b[2L])) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[2L])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$set_fit_warm_start(
				c(as.numeric(attempt$fit$b), if (is.finite(attempt$fit$phi)) log(as.numeric(attempt$fit$phi)) else 0),
				"params",
				fisher = attempt$fit$fisher_information
			)
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		best_X_colnames = NULL,
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
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			n_params = ncol(X) + 1L
			ws_args = private$get_backend_warm_start_args(n_params)
			res = fast_beta_regression_cpp(
				X = X, y = sanitize_beta_response(as.numeric(private$y)),
				warm_start_beta = ws_args$start_beta,
				warm_start_fisher_info = ws_args$warm_start_fisher_info,
				compute_std_errs = FALSE,
				smart_cold_start = private$smart_cold_start_default,
				optimization_alg = private$optimization_alg
			)

			if (is.null(res) || !is.finite(res$coefficients[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(c(as.numeric(res$coefficients), log(as.numeric(res$phi))), "params", fisher = res$fisher_information)
			as.numeric(res$coefficients[2])
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
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$b)
			phi    = as.numeric(null_fit$phi)
			if (!is.finite(phi) || phi <= 0) return(NULL)
			mu     = pmin(pmax(plogis(as.numeric(spec$X %*% b_null)), 1e-8), 1 - 1e-8)
			y_sim  = rbeta(length(mu), shape1 = mu * phi, shape2 = (1 - mu) * phi)
			y_sim  = pmin(pmax(y_sim, 1e-8), 1 - 1e-8)
			X_fit  = spec$X
			j      = spec$j
			
			# Parametric bootstrap: use observed fit as anchor for the full fit
			ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
			full_res = tryCatch(
				fast_beta_regression_cpp(
					X = X_fit, y = y_sim,
					warm_start_beta = ws_args$start_beta,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || length(full_res$coefficients) < j || !is.finite(full_res$coefficients[j])) return(NULL)
			full_fit_boot = list(
				b          = as.numeric(full_res$coefficients),
				phi        = as.numeric(full_res$phi),
				neg_loglik = as.numeric(full_res$neg_loglik)
			)
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					ws_args_null = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					res = tryCatch(
						fast_beta_regression_cpp(
							X = X_fit, y = y_sim,
							warm_start_beta = (if (!is.null(start)) start[seq_len(ncol(X_fit))] else full_fit_boot$b),
							warm_start_fisher_info = ws_args_null$warm_start_fisher_info,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					list(b = as.numeric(res$coefficients), phi = as.numeric(res$phi), neg_loglik = as.numeric(res$neg_loglik))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = sanitize_beta_response(as.numeric(private$y))
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					res = fast_beta_regression_cpp(
						X_fit, y,
						fixed_idx = j_treat, fixed_values = delta,
						warm_start_beta = if (!is.null(start)) start[1:ncol(X_fit)] else ws_args$start_beta,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						optimization_alg = private$optimization_alg
					)
					if (is.null(res)) return(NULL)
					list(b = as.numeric(res$coefficients), phi = res$phi, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
				},
				extract_start = function(fit){
					c(as.numeric(fit$b), log(as.numeric(fit$phi)))
				},
				score = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
					get_beta_regression_score_cpp(X_fit, y, params)
				},
				observed_information = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
					-get_beta_regression_hessian_cpp(X_fit, y, params)
				},
				fisher_information = function(fit){
					fit$fisher_information %||% {
						params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
						-get_beta_regression_hessian_cpp(X_fit, y, params)
					}
				},
				information = function(fit){
					fit$information %||% fit$fisher_information %||% {
						params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
						-get_beta_regression_hessian_cpp(X_fit, y, params)
					}
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			y_san = sanitize_beta_response(as.numeric(private$y))
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit){
					n_params = ncol(X_fit) + 1L
					ws_args = private$get_backend_warm_start_args(n_params)
					if (estimate_only) {
						res = fast_beta_regression_cpp(
							X_fit, y_san,
							warm_start_beta = ws_args$start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							compute_std_errs = FALSE,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res)) return(NULL)
						list(b = res$coefficients, ssq_b_2 = NA_real_, phi = res$phi, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
					} else {
						res = fast_beta_regression_with_var_cpp(
							X_fit, y_san,
							warm_start_beta = ws_args$start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res)) return(NULL)
						list(b = res$coefficients,
						     ssq_b_2 = if (!is.null(res$vcov) && nrow(res$vcov) >= 2L) res$vcov[2L, 2L] else NA_real_,
						     phi = res$phi, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
					}
				},

				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && length(mod$b) >= 2L && is.finite(mod$b[2L])
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$set_fit_warm_start(c(as.numeric(attempt$fit$b), log(as.numeric(attempt$fit$phi))), "params", fisher = attempt$fit$fisher_information)
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = which(attempt$keep == 2L)
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		}
	)
)
