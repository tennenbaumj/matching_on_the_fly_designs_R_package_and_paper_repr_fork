#' Count-Specific Likelihood Inference
#'
#' @name InferenceCountLikelihood
#' @description Intermediate base class for count-based likelihood families
#' (Poisson, Negative Binomial, Zero-Inflated, Hurdle). This class centralizes
#' count-specific parameter packing, warm starts, and likelihood dispatch.
#'
#' @keywords internal
inference_count_likelihood_public = list(
		#' @description Computes the treatment estimate using the underlying model.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an asymptotic confidence interval using the configured test.
		#' @param alpha Significance level 1 - \code{alpha}. Default 0.05.
		#' @return A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			if (private$testing_type == "wald") {
				private$shared(estimate_only = FALSE)
				if (is.finite(private$cached_values$s_beta_hat_T %||% NA_real_)) {
					return(private$compute_z_or_t_ci_from_s_and_df(alpha))
				}
			}
			super$compute_asymp_confidence_interval(alpha)
		},
		#' @description Computes an asymptotic two-sided p-value using the configured test.
		#' @param delta Null treatment effect to test against. Default 0.
		#' @return The asymptotic p-value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			if (private$testing_type == "wald") {
				private$shared(estimate_only = FALSE)
				if (is.finite(private$cached_values$s_beta_hat_T %||% NA_real_)) {
					return(private$compute_z_or_t_two_sided_pval_from_s_and_df(delta))
				}
			}
			super$compute_asymp_two_sided_pval(delta)
		},
		#' @description Computes the Wald two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_wald_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_wald_two_sided_pval(delta)
		},
		#' @description Computes the Wald confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_wald_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_wald_confidence_interval(alpha)
		},
		#' @description Computes the score two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_score_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_score_two_sided_pval(delta)
		},
		#' @description Computes the score confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_score_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_score_confidence_interval(alpha)
		},
		#' @description Computes the likelihood-ratio two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_lik_ratio_two_sided_pval(delta)
		},
		#' @description Computes the likelihood-ratio confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_lik_ratio_confidence_interval(alpha)
		},
		#' @description Computes the gradient two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_gradient_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_gradient_two_sided_pval(delta)
		},
		#' @description Computes the gradient confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_gradient_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_gradient_confidence_interval(alpha)
		},
		#' @description Computes a parametric likelihood-ratio bootstrap p-value.
		#' @param delta Null treatment effect.
		#' @param B Number of bootstrap replicates.
		#' @param show_progress Whether to show progress.
		#' @param min_number_usable_samples Minimum usable bootstrap samples.
		#' @param max_attempts_per_replicate Maximum attempts per replicate.
		compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_lik_ratio_bootstrap_two_sided_pval(
				delta = delta,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate
			)
		},
		#' @description Computes a parametric likelihood-ratio bootstrap CI.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of bootstrap replicates.
		#' @param show_progress Whether to show progress.
		#' @param min_number_usable_samples Minimum usable bootstrap samples.
		#' @param max_attempts_per_replicate Maximum attempts per replicate.
		#' @param root_tolerance Root tolerance.
		#' @param max_root_iterations Maximum root iterations.
		compute_lik_ratio_bootstrap_confidence_interval = function(alpha = 0.05, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L, root_tolerance = NULL, max_root_iterations = 8L){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_lik_ratio_bootstrap_confidence_interval(
				alpha = alpha,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate,
				root_tolerance = root_tolerance,
				max_root_iterations = max_root_iterations
			)
		}
	)

inference_count_likelihood_private = list(
		# --- Count-specific shared logic ---

		shared = function(estimate_only = FALSE){
			if (estimate_only &&
					!is.null(private$cached_values$beta_hat_T) &&
					is.finite(as.numeric(private$cached_values$beta_hat_T)[1L])) {
				return(invisible(NULL))
			}
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T) && is.finite(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			model_output = private$generate_mod(estimate_only = estimate_only)
			private$cached_mod = model_output

			if (is.null(model_output)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}

			# Count models usually return treatment effect at index 2 of the cond component
			# or in a specific field. We allow models to set beta_hat_T directly in output.
			private$cached_values$beta_hat_T = model_output$beta_hat_T %||% model_output$b[2]

			if (!is.null(model_output$params) || !is.null(model_output$b)) {
				private$set_fit_warm_start(
					as.numeric(model_output$params %||% model_output$b),
					type = if (!is.null(model_output$params)) "params" else "beta",
					fisher = model_output$fisher_information %||% model_output$XtWX,
					weights = model_output$w %||% model_output$mu,
					force_pd = TRUE
				)
			}

			if (estimate_only) return(invisible(NULL))

			ssq = model_output$ssq_b_j %||% model_output$ssq_b_2
			if (!is.null(ssq) && is.finite(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = model_output$df %||% Inf
		},

		get_standard_error = function(){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$shared(estimate_only = FALSE)
			# Try information-based SE first if supported
			if (isTRUE(private$supports_information_preference())) {
				se = tryCatch(private$compute_standard_error_from_information_matrix(), error = function(e) NA_real_)
				if (is.finite(se)) return(se)
			}
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df %||% Inf
		},
		get_backend_warm_start_args = function(expected_length, expected_fisher_dim = expected_length) {
			private$get_optimal_warm_start_config(expected_length, expected_fisher_dim)
		},

		# --- Likelihood test support ---

		supports_lik_ratio_param_bootstrap = function() TRUE,
		supports_likelihood_tests = function(){
			TRUE
		},

		get_likelihood_test_spec = function(){
			# This is still abstract, but we provide the structure
			NULL
		},

		compute_score_two_sided_pval_impl = function(delta){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		compute_score_confidence_interval_impl = function(alpha){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			private$invert_test_pval_confidence_interval(alpha)
		},
		compute_gradient_confidence_interval_impl = function(alpha){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			private$invert_test_pval_confidence_interval(alpha)
		},
		compute_lik_ratio_confidence_interval_impl = function(alpha){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			private$invert_test_pval_confidence_interval(alpha)
		},
		count_likelihood_block_asymp_unsupported = function(){
			private$jackknife_block_size_gt_one_unsupported(unit = "auto")
		},
		mark_count_likelihood_block_asymp_nonestimable = function(){
			if (!private$count_likelihood_block_asymp_unsupported()) return(FALSE)
			if (is.null(private$cached_values$beta_hat_T)) {
				try(private$shared(estimate_only = TRUE), silent = TRUE)
			}
			private$cache_nonestimable_se("count_likelihood_asymp_block_size_gt_one_not_supported")
			TRUE
		},
		count_likelihood_missing_ci = function(alpha = 0.05){
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	)

InferenceCountLikelihood = R6::R6Class("InferenceCountLikelihood",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = inference_count_likelihood_public,
	private = inference_count_likelihood_private
)

#' Count-Specific Likelihood Inference Without Parametric LR Bootstrap
#'
#' @name InferenceCountLikelihoodNoParamBootstrap
#' @description Internal base class for count likelihood families that should
#' not expose the parametric LR bootstrap API.
#'
#' @keywords internal
#' @noRd
InferenceCountLikelihoodNoParamBootstrap = R6::R6Class("InferenceCountLikelihoodNoParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = inference_count_likelihood_public,
	private = inference_count_likelihood_private
)
