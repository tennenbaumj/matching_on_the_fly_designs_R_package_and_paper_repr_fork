#' Abstract class for all-subject modified-Poisson inference in KK designs
#'
#' @keywords internal
InferenceAbstractKKModifiedPoisson = R6::R6Class("InferenceAbstractKKModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKMarginalIncid,
	public = list(
		#' @description Compute treatment estimate
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the KK marginal incidence estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Numeric vector. Row weights for bootstrap.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$beta_hat_T = beta_hat_T
					private$cached_values$s_beta_hat_T = NA_real_
					return(private$cached_values$beta_hat_T)
				}
			}
			X = private$build_design_matrix()
			if (is.null(X)) return(NA_real_)
			X = as.matrix(X)
			ok = is.finite(row_weights) & row_weights > 0 & is.finite(as.numeric(private$y))
			if (!any(ok)) return(NA_real_)
			X_fit = X[ok, , drop = FALSE]
			y_fit = as.numeric(private$y[ok])
			w_fit = as.numeric(row_weights[ok])
			mod = tryCatch(
				fast_poisson_regression_weighted_cpp(
					X = X_fit,
					y = y_fit,
					weights = w_fit,
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
					smart_cold_start = private$smart_cold_start_default,
					optimization_alg = private$optimization_alg %||% "irls"
				),
				error = function(e) NULL
			)
			beta_hat_T = if (is.null(mod) || length(mod$b) < 2L) NA_real_ else as.numeric(mod$b[2L])
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Compute asymp confidence interval
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute asymp two sided pval for treatment effect
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		max_abs_reasonable_coef = 25,
		best_X_colnames = NULL,
		cached_mod = NULL,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			# Use the same design matrix structure as the original fit
			X_cols = private$best_X_colnames
			X_data = private$get_X()
			
			if (length(X_cols) == 0L){
				# Univariate case
				X = cbind("(Intercept)" = 1, treatment = private$w)
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind("(Intercept)" = 1, treatment = private$w, X_cov)
			}
			fit = tryCatch(private$fit_modified_poisson(X, j_treat = 2L, estimate_only = estimate_only), error = function(e) NULL)
			if (is.null(fit) || !is.finite(fit$beta_hat)){
				return(NA_real_)
			}
			as.numeric(fit$beta_hat)
		},
		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},
		set_failed_fit_cache = function(){
			private$cache_nonestimable_estimate("kk_modified_poisson_fit_unavailable")
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},
		coefficients_are_usable = function(coef_hat){
			length(coef_hat) > 0L &&
				all(is.finite(coef_hat)) &&
				max(abs(coef_hat), na.rm = TRUE) <= private$max_abs_reasonable_coef
		},
		fit_modified_poisson = function(X_fit, j_treat, estimate_only = FALSE){
			mod = tryCatch(
				fast_poisson_regression_cpp(
					X = X_fit, 
					y = as.numeric(private$y),
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
				),
				error = function(e) NULL
			)
			if (is.null(mod)){
				return(NULL)
			}
			coef_hat = as.numeric(mod$b)
			if (length(coef_hat) != ncol(X_fit) || !private$coefficients_are_usable(coef_hat)){
				return(NULL)
			}
			if (estimate_only){
				return(list(
					beta_hat = coef_hat[j_treat],
					se = NA_real_,
					coefficients = coef_hat,
					vcov = NULL,
					summary_table = NULL,
					mod = mod
				))
			}
			mu_hat = as.numeric(exp(X_fit %*% coef_hat))
			if (length(mu_hat) != nrow(X_fit) || any(!is.finite(mu_hat)) || any(mu_hat <= 0)){
				return(NULL)
			}
			post_fit = tryCatch(
				glm_cluster_sandwich_post_fit_cpp(
					X_fit = X_fit,
					y = as.numeric(private$y),
					coef_hat = coef_hat,
					mu_hat = mu_hat,
					working_weights = mu_hat,
					cluster_id = private$get_cluster_ids(),
					j_treat = j_treat
				),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				return(NULL)
			}
			coef_names = colnames(X_fit)
			names(coef_hat) = coef_names
			vcov_robust = post_fit$vcov
			colnames(vcov_robust) = rownames(vcov_robust) = coef_names
			std_err = post_fit$std_err
			names(std_err) = coef_names
			z_vals = post_fit$z_vals
			names(z_vals) = coef_names
			list(
				beta_hat = post_fit$beta_hat,
				se = post_fit$se,
				coefficients = coef_hat,
				vcov = vcov_robust,
				mod = mod,
				summary_table = cbind(
					Value = coef_hat,
					`Std. Error` = std_err,
					`z value` = z_vals,
					`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
				)
			)
		},
		get_standard_error = function(){
			private$shared()
			se = private$cached_values$s_beta_hat_T
			if (is.null(se) || length(se) == 0L) {
				return(NA_real_)
			}
			as.numeric(se)[1L]
		},
		get_degrees_of_freedom = function(){
			private$shared()
			private$cached_values$df %||% Inf
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		supports_lik_ratio_param_bootstrap = function(){
			TRUE
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null     = as.numeric(null_fit$b)
			mu         = pmax(exp(as.numeric(spec$X %*% b_null)), 0)
			y_sim      = as.numeric(rpois(length(mu), mu))
			X_fit      = spec$X
			j          = spec$j
			full_fit_b = tryCatch(
				fast_poisson_regression_cpp(
					X = X_fit, y = y_sim,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(full_fit_b) || length(full_fit_b$b) < j || !is.finite(full_fit_b$b[j])) return(NULL)
			list(
				full_fit = full_fit_b,
				fit_null = function(d, start = NULL){
					tryCatch(
						fast_poisson_regression_with_var_cpp(
							X = X_fit, y = y_sim, j = j,
							warm_start_beta = start %||% full_fit_b$b,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit){
					eta_f = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y_sim * eta_f - exp(eta_f) - lgamma(y_sim + 1))
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
					fast_poisson_regression_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						warm_start_beta = start %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_cold_start = private$smart_cold_start_default
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					as.numeric(fit$score %||% get_poisson_regression_score_cpp(X_fit, y, as.numeric(fit$b)))
				},
				observed_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y * eta - exp(eta) - lgamma(y + 1))
				}
			)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			if (is.null(colnames(X_full))) {
				colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)
			}
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}
			fit = private$fit_modified_poisson(X_fit, j_treat, estimate_only = estimate_only)
			if (!is.null(fit)) {
				private$best_X_colnames = setdiff(colnames(X_fit), c("(Intercept)", "treatment"))
			}
			if (private$harden && is.null(fit) && ncol(X_full) > 2L){
				reduced = private$reduce_design_matrix_preserving_treatment(X_full[, 1:2, drop = FALSE])
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (!is.null(X_fit) && is.finite(j_treat) && nrow(X_fit) > ncol(X_fit)){
					fit = private$fit_modified_poisson(X_fit, j_treat, estimate_only = estimate_only)
				}
			}
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}
			private$cached_mod = fit$mod %||% fit
			# Save warm start
			private$set_fit_warm_start(private$cached_mod$b, "beta", fisher = private$cached_mod$fisher_information)
			
			private$cached_values$likelihood_test_context = list(
				X = X_fit,
				j_treat = j_treat
			)
			private$cached_values$beta_hat_T = fit$beta_hat
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$full_vcov = fit$vcov
			private$cached_values$summary_table = fit$summary_table
		}
	)
)
#' G-Computation Risk-Difference Inference for KK Designs with Binary Responses
#'
#' Fits an all-subject logistic working model for a KK incidence outcome using
#' treatment and, optionally, all recorded covariates, then estimates the marginal
#' risk difference by standardizing predicted risks under all-treated and
#' all-control assignments over the empirical covariate distribution. Matched
#' pairs are treated as clusters and reservoir subjects are treated as singletons
#' when computing the sandwich covariance.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKGCompRiskDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKGCompRiskDiff = R6::R6Class("InferenceIncidKKGCompRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidKKGCompAbstract,
	public = list(
	),
	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		},
		get_estimand_type = function() "RD"
	)
)
#' G-Computation Risk-Ratio Inference for KK Designs with Binary Responses
#'
#' Fits a all-subject logistic working model for a KK incidence outcome using
#' treatment and, optionally, all recorded covariates, then estimates the marginal
#' risk ratio by standardizing predicted risks under all-treated and all-control
#' assignments over the empirical covariate distribution. Matched pairs are
#' treated as clusters and reservoir subjects are treated as singletons when
#' computing the sandwich covariance.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKGCompRiskRatio$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKGCompRiskRatio = R6::R6Class("InferenceIncidKKGCompRiskRatio",
	lock_objects = FALSE,
	inherit = InferenceIncidKKGCompAbstract,
	public = list(
	),
	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		},
		get_estimand_type = function() "RR"
	)
)
#' Modified-Poisson Inference for KK Designs with Binary Responses
#'
#' Fits an all-subject modified-Poisson working model for incidence outcomes under
#' a KK matching-on-the-fly design using treatment and, optionally, all recorded
#' covariates as predictors. Matched pairs are treated as clusters and reservoir
#' subjects are treated as singleton clusters when computing the sandwich
#' covariance.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKModifiedPoisson$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKModifiedPoisson = R6::R6Class("InferenceIncidKKModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKModifiedPoisson,
	public = list(
	),
	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
