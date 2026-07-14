#' Weibull AFT Inference for Survival Responses
#'
#' Fits a Weibull Accelerated Failure Time (AFT) model for survival responses
#' using the treatment indicator and, optionally, all recorded covariates as
#' predictors. The treatment effect is reported on the log-time-ratio scale.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalWeibullRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' \donttest{
#' inf$set_seed(1)
#' inf$compute_lik_ratio_bootstrap_two_sided_pval(delta = 0, B = 9, show_progress = FALSE)
#' }
#' @export
InferenceSurvivalWeibullRegr = R6::R6Class("InferenceSurvivalWeibullRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a Weibull-regression inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart optimizer start values by default.
		#' @param optimization_alg Character scalar specifying the optimization algorithm. 
		#'   Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFormula(model_formula, null.ok = TRUE)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the Weibull AFT treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) return(beta_hat_T)
			}
			X_fit = private$build_design_matrix()[, -1, drop = FALSE]
			colnames(X_fit)[1L] = "treatment"
			fit = weighted_weibull_bootstrap_surrogate_fit(private$y, private$dead, X_fit, row_weights)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		best_X_colnames = NULL,
		get_complexity_tier = function() "light",
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
			
			res = fast_weibull_regression_cpp(
				y = private$y, dead = private$dead, X = X,
				warm_start_params = ws_args$start_params,
				warm_start_fisher_info = ws_args$warm_start_fisher_info,
				smart_cold_start = private$smart_cold_start_default,
				estimate_only = TRUE, optimization_alg = private$optimization_alg
			)
			if (is.null(res) || !isTRUE(res$converged) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(c(as.numeric(res$b), as.numeric(res$log_sigma)), "params", fisher = res$fisher_information)
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
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null    = as.numeric(null_fit$b)
			log_sigma = as.numeric(null_fit$log_sigma)
			X_fit    = spec$X
			j        = spec$j
			sim_data = private$simulate_param_boot_weibull_observed(
				X = X_fit,
				b_null = b_null,
				log_sigma = log_sigma,
				y_obs = private$y,
				dead_obs = private$dead
			)
			if (is.null(sim_data)) return(NULL)
			y_sim = sim_data$y
			dead_sim = sim_data$dead
			full_res = tryCatch(
				fast_weibull_regression_cpp(
					y = y_sim, dead = dead_sim, X = X_fit,
					smart_cold_start = private$smart_cold_start_default,
					estimate_only = FALSE, optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			p_fit = ncol(X_fit)
			if (is.null(full_res) || !isTRUE(full_res$converged) || length(full_res$params) < j || !is.finite(full_res$params[j])) return(NULL)
			full_fit_boot = list(
				b          = as.numeric(full_res$params[seq_len(p_fit)]),
				log_sigma  = as.numeric(full_res$params[p_fit + 1L]),
				neg_loglik = as.numeric(full_res$neg_ll)
			)
			list(
				worker_data = list(y = y_sim, dead = dead_sim),
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = tryCatch(
						fast_weibull_regression_cpp(
							y = y_sim, dead = dead_sim, X = X_fit,
							warm_start_params = start %||% as.numeric(full_res$params),
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE,
							estimate_only = FALSE, optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$params[seq_len(p_fit)]), log_sigma = as.numeric(res$params[p_fit + 1L]), neg_loglik = as.numeric(res$neg_ll))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			dead = as.numeric(private$dead)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					res = tryCatch(
						fast_weibull_regression_cpp(
							y = y, dead = dead, X = X_fit,
							warm_start_params = start %||% ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							fixed_idx = j_treat, fixed_values = delta,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$params[seq_len(ncol(X_fit))]), log_sigma = as.numeric(res$params[ncol(X_fit) + 1L]), neg_loglik = as.numeric(res$neg_ll), fisher_information = res$information)
				},
				extract_start = function(fit){
					c(as.numeric(fit$b), as.numeric(fit$log_sigma))
				},
				score = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					get_weibull_regression_score_cpp(X_fit, y, dead, params)
				},
				observed_information = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					-get_weibull_regression_hessian_cpp(X_fit, y, dead, params)
				},
				fisher_information = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					-get_weibull_regression_hessian_cpp(X_fit, y, dead, params)
				},
				information = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					-get_weibull_regression_hessian_cpp(X_fit, y, dead, params)
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit){
					n_params = ncol(X_fit) + 1L
					ws_args = private$get_backend_warm_start_args(n_params)
					res = fast_weibull_regression_cpp(
						y = private$y, dead = private$dead, X = X_fit,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					p = ncol(X_fit)
					b_vals = as.numeric(if (estimate_only) res$b[seq_len(p)] else res$params[seq_len(p)])
					ls_val = as.numeric(if (estimate_only) res$log_sigma else res$params[p + 1L])
					list(
						b                  = b_vals,
						params             = c(b_vals, ls_val),
						log_sigma          = ls_val,
						fisher_information = if (estimate_only) NULL else res$information,
						neg_loglik         = as.numeric(res$neg_ll),
						ssq_b_2            = if (estimate_only || is.null(res$vcov)) NA_real_ else {
							if (nrow(res$vcov) >= 2L) res$vcov[2L, 2L] else NA_real_
						}
					)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2]) || abs(mod$b[2]) > 5) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$params, "params", fisher = attempt$fit$fisher_information)
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
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
