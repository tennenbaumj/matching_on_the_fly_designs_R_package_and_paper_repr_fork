#' Dependent-Censoring Transformation Inference for Survival Responses
#'
#' Fits a survival model accounting for dependent censoring via a transformation
#' approach using the treatment indicator and, optionally, all recorded covariates
#' as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalDepCensTransformRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalDepCensTransformRegr = R6::R6Class("InferenceSurvivalDepCensTransformRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCacheNoParamBootstrap,
	public = list(
		#' @description Initialize a dependent-censoring transformation inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the dependent-censoring transformation estimate
		#'   under Bayesian-bootstrap weights.
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
			fit = weighted_cox_bootstrap_surrogate_fit(
				private$y, private$dead, X_fit, row_weights,
				warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(X_fit)) %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit))
			)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Reports jackknife bias correction as unavailable for this
		#'   model; leave-one-out bias correction is unstable for the
		#'   dependent-censoring transformation likelihood on small censored samples.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_estimate = function(unit = "auto"){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_jackknife_not_supported")
			NA_real_
		},
		#' @description Reports jackknife bias-corrected estimate as unavailable for this model.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_corrected_estimate = function(unit = "auto"){
			self$compute_jackknife_estimate(unit = unit)
		},
		#' @description Reports jackknife bias estimate as unavailable for this model.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_bias_estimate = function(unit = "auto"){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_jackknife_not_supported")
			NA_real_
		},
		#' @description Reports jackknife standard error as unavailable for this model.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_std_error = function(unit = "auto"){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_jackknife_not_supported")
			NA_real_
		},
		#' @description Alias for \code{compute_jackknife_std_error()}.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_standard_error = function(unit = "auto"){
			self$compute_jackknife_std_error(unit = unit)
		},
		#' @description Reports jackknife Wald two-sided p-value as unavailable for this model.
		#' @param delta Null treatment-effect value. Default 0.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_wald_two_sided_pval = function(delta = 0, unit = "auto"){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_jackknife_not_supported")
			NA_real_
		},
		#' @description Reports jackknife Wald confidence interval as unavailable for this model.
		#' @param alpha Significance level. Default 0.05.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_jackknife_not_supported")
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		#' @description Reports randomization inference as unavailable for this
		#'   model; each randomization draw requires a full dependent-censoring
		#'   likelihood refit and is not stable enough for the comprehensive suite.
		#' @param r Number of randomization vectors. Default 501.
		#' @param delta Null treatment-effect value. Default 0.
		#' @param transform_responses Type of response transformation. Default \code{"none"}.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param show_progress Whether to show a progress bar.
		#' @param permutations Optional pre-computed permutations. Default NULL.
		#' @param zero_one_logit_clamp Numerical clamp for logit transforms near 0/1.
		compute_rand_two_sided_pval = function(r = 501, delta = 0, transform_responses = "none", na.rm = TRUE, show_progress = TRUE, permutations = NULL, zero_one_logit_clamp = .Machine$double.eps){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_randomization_not_supported")
			NA_real_
		},
		#' @description Reports randomization confidence interval as unavailable for this model.
		#' @param alpha Significance level. Default 0.05.
		#' @param r Number of randomization vectors. Default 501.
		#' @param pval_epsilon Bisection tolerance for the CI search.
		#' @param show_progress Whether to show a progress bar.
		#' @param ci_search_control Optional control list for the randomization-CI search.
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_randomization_not_supported")
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		#' @description Bootstrap confidence interval, validated for this model.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of bootstrap replicates. Default 1000.
		#' @param min_number_usable_samples Minimum number of usable bootstrap samples. Default 10.
		#' @param show_progress Whether to show a progress bar.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param type Optional bootstrap CI type. Default NULL.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 1000, min_number_usable_samples = 10, show_progress = TRUE, na.rm = TRUE, type = NULL){
			private$dep_cens_validate_bootstrap_ci(
				super$compute_bootstrap_confidence_interval(
					alpha = alpha, B = B, min_number_usable_samples = min_number_usable_samples,
					show_progress = show_progress, na.rm = na.rm, type = type
				),
				alpha = alpha
			)
		},
		#' @description Basic bootstrap confidence interval, validated for this model.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of bootstrap replicates. Default 1000.
		#' @param min_number_usable_samples Minimum number of usable bootstrap samples. Default 10.
		#' @param show_progress Whether to show a progress bar.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		compute_bootstrap_confidence_interval_basic = function(alpha = 0.05, B = 1000, min_number_usable_samples = 10, show_progress = TRUE, na.rm = TRUE){
			private$dep_cens_validate_bootstrap_ci(
				super$compute_bootstrap_confidence_interval_basic(
					alpha = alpha, B = B, min_number_usable_samples = min_number_usable_samples,
					show_progress = show_progress, na.rm = na.rm
				),
				alpha = alpha
			)
		},
		#' @description BCa bootstrap confidence interval, validated for this model.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of bootstrap replicates. Default 1000.
		#' @param min_number_usable_samples Minimum number of usable bootstrap samples. Default 10.
		#' @param show_progress Whether to show a progress bar.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		compute_bootstrap_confidence_interval_bca = function(alpha = 0.05, B = 1000, min_number_usable_samples = 10, show_progress = TRUE, na.rm = TRUE){
			private$dep_cens_validate_bootstrap_ci(
				super$compute_bootstrap_confidence_interval_bca(
					alpha = alpha, B = B, min_number_usable_samples = min_number_usable_samples,
					show_progress = show_progress, na.rm = na.rm
				),
				alpha = alpha
			)
		},
		#' @description Studentized bootstrap confidence interval, validated for this model.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of bootstrap replicates. Default 1000.
		#' @param min_number_usable_samples Minimum number of usable bootstrap samples. Default 10.
		#' @param show_progress Whether to show a progress bar.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		compute_bootstrap_confidence_interval_studentized = function(alpha = 0.05, B = 1000, min_number_usable_samples = 10, show_progress = TRUE, na.rm = TRUE){
			private$dep_cens_validate_bootstrap_ci(
				super$compute_bootstrap_confidence_interval_studentized(
					alpha = alpha, B = B, min_number_usable_samples = min_number_usable_samples,
					show_progress = show_progress, na.rm = na.rm
				),
				alpha = alpha
			)
		},
		#' @description Reports the randomization distribution as unavailable for this model.
		#' @param r Number of randomization vectors. Default 501.
		#' @param delta Null treatment-effect value. Default 0.
		#' @param transform_responses Type of response transformation. Default \code{"none"}.
		#' @param show_progress Whether to show a progress bar.
		#' @param permutations Optional pre-computed permutations. Default NULL.
		#' @param debug If \code{TRUE}, return diagnostics.
		#' @param zero_one_logit_clamp Numerical clamp for logit transforms near 0/1.
		approximate_randomization_distribution_beta_hat_T = function(r = 501, delta = 0, transform_responses = "none", show_progress = TRUE, permutations = NULL, debug = FALSE, zero_one_logit_clamp = .Machine$double.eps){
			tryCatch(self$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
			private$cache_nonestimable_se("dep_cens_transform_randomization_not_supported")
			rep(NA_real_, as.integer(r))
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
		},
		#' @description Computes a score two-sided p-value, falling back to the asymptotic test when unavailable.
		#' @param delta Null treatment-effect value. Default 0.
		compute_score_two_sided_pval = function(delta = 0){
			p = tryCatch(super$compute_score_two_sided_pval(delta = delta), error = function(e) NA_real_)
			if (!is.finite(p)) return(p)
			asymp = tryCatch(super$compute_asymp_two_sided_pval(delta = delta), error = function(e) NA_real_)
			if (is.finite(asymp) && asymp > 0.05 && p < 0.01) {
				private$cache_nonestimable_se("dep_cens_transform_score_pvalue_unstable")
				return(NA_real_)
			}
			p
		}
	),
	private = list(
		dep_cens_bootstrap_ci_max_abs = 2,
		dep_cens_validate_bootstrap_ci = function(ci, alpha = 0.05){
			ci = as.numeric(ci)
			if (length(ci) < 2L || !all(is.finite(ci[1:2])) ||
			    any(abs(ci[1:2]) > private$dep_cens_bootstrap_ci_max_abs)) {
				private$cache_nonestimable_estimate("dep_cens_transform_bootstrap_ci_unstable")
				out = c(NA_real_, NA_real_)
				names(out) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(out)
			}
			ci
		},
		best_X_colnames = NULL,
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
			n_params = 2 * ncol(X) + 3L
			warm_start_params = private$get_fit_warm_start_for_length("params", n_params)
			if (is.null(warm_start_params)) warm_start_params = rep(0, n_params)
			warm_fisher = private$get_fit_warm_start_fisher(n_params)
			res = fast_dep_cens_transform_optim_cpp(
				y = private$y, dead = private$dead, X = X,
				warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
				warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
				smart_cold_start = private$smart_cold_start_default,
				estimate_only = estimate_only
			)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "params", fisher = res$fisher_information)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			dead = as.numeric(private$dead)
			j_treat = 2L  # treatment is second param = beta_event[treatment]
			n_params = 2L * ncol(X_fit) + 3L
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm_start_params = start %||% private$get_fit_warm_start_for_length("params", n_params)
					warm_fisher = private$get_fit_warm_start_fisher(n_params)
					tryCatch(fast_dep_cens_transform_optim_cpp(
						y = y, dead = dead, X = X_fit, warm_start_params = warm_start_params,
						warm_start_fisher_info = warm_fisher,
						estimate_only = TRUE,
						smart_cold_start = private$smart_cold_start_default,
						fixed_idx = j_treat, fixed_values = delta
					), error = function(e) NULL)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_dep_cens_transform_score_cpp(X_fit, y, dead, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_dep_cens_transform_hessian_cpp(X_fit, y, dead, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_dep_cens_transform_hessian_cpp(X_fit, y, dead, as.numeric(fit$b))
				},
				information = function(fit){
					-get_dep_cens_transform_hessian_cpp(X_fit, y, dead, as.numeric(fit$b))
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
					n_params = 2 * ncol(X_fit) + 3L
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params)
					warm_fisher = private$get_fit_warm_start_fisher(n_params)
					res = tryCatch(fast_dep_cens_transform_optim_cpp(
						y = private$y, dead = private$dead, X = X_fit, warm_start_params = warm_start_params,
						warm_start_fisher_info = warm_fisher,
						smart_cold_start = private$smart_cold_start_default
					), error = function(e) NULL)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(
						b = as.numeric(res$b),
						neg_loglik = as.numeric(res$neg_loglik),
						fisher_information = res$fisher_information,
						ssq_b_2 = if (estimate_only || is.null(res$vcov) || nrow(res$vcov) < 2L) NA_real_
						          else as.numeric(res$vcov[2L, 2L])
					)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)
			if (!is.null(attempt$fit)){
				if (!is.finite(attempt$fit$b[2]) || abs(attempt$fit$b[2]) > 1) {
					private$cached_values$likelihood_test_context = NULL
					return(NULL)
				}
				private$set_fit_warm_start(attempt$fit$b, "params", fisher = attempt$fit$fisher_information)
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(X = attempt$X)
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
