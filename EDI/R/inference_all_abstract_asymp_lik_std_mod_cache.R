#' GLM and Kaplan-Meier Inference
#'
#' Abstract class providing MLE/KM-based inference methods for GLM and survival models.
#'
#' @keywords internal
inference_asymp_lik_std_mod_cache_public = list(
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
			if (private$testing_type == "wald") {
				private$shared(estimate_only = FALSE)
				if (is.finite(private$cached_values$s_beta_hat_T %||% NA_real_)) {
					return(private$compute_z_or_t_two_sided_pval_from_s_and_df(delta))
				}
			}
			super$compute_asymp_two_sided_pval(delta)
		}
	)
inference_asymp_lik_std_mod_cache_private = list(
		supports_likelihood_tests = function(){
			TRUE
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		generate_mod = function(estimate_only = FALSE) stop(class(self)[1], " must implement generate_mod()"),
		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},
		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},
		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			if (isTRUE(private$supports_information_preference())) {
				se = private$compute_standard_error_from_information_matrix()
				if (is.finite(se)) return(se)
			}
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		},
		compute_score_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		get_likelihood_test_spec = function(){
			NULL
		},
		make_warm_fit_null_wrapper = function(spec, cache_key){
			last_start = NULL
			last_delta = NULL
			fit_null_formals = tryCatch(names(formals(spec$fit_null)), error = function(e) character())
			accepts_start = "start" %in% fit_null_formals
			function(delta){
				warm_enabled = isTRUE(private$null_fit_warm_start_enabled)
				cache_state = if (warm_enabled) private$get_likelihood_null_warm_state(cache_key) else NULL
				start = if (warm_enabled) last_start else NULL
				if (warm_enabled && is.null(start) && !is.null(cache_state)) start = cache_state$start
				fit = tryCatch(
					if (accepts_start) spec$fit_null(delta, start = start) else spec$fit_null(delta),
					error = function(e) NULL
				)
				extract_start = spec$extract_start %||% function(fit_obj) NULL
				last_start <<- if (warm_enabled && accepts_start) tryCatch(extract_start(fit), error = function(e) NULL) else NULL
				last_delta <<- delta
				if (warm_enabled && accepts_start) {
					private$set_likelihood_null_warm_state(cache_key, delta = delta, start = last_start)
				}
				fit
			}
		},
		compute_likelihood_test_two_sided_pval = function(delta, testing_type){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) {
				stop(class(self)[1], " does not expose a likelihood-test specification.", call. = FALSE)
			}
			p_value = private$get_memoized_likelihood_test_pval(
				delta = delta,
				testing_type = testing_type,
				spec = spec,
				warm_cache_key = paste0("likelihood_test:", testing_type)
			)
			if (!is.finite(p_value) && !isTRUE(self$is_nonestimable("estimate"))) {
				private$cache_nonestimable_se(paste0(testing_type, "_test_unavailable"))
			}
			p_value
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			has_cached_se = !is.null(private$cached_values$s_beta_hat_T) &&
				length(private$cached_values$s_beta_hat_T) == 1L &&
				isTRUE(is.finite(private$cached_values$s_beta_hat_T))
			if (isTRUE(!is.null(private$cached_values$beta_hat_T) && (estimate_only || has_cached_se))) return(invisible(NULL))
			
			# Abstract function implemented by daughter classes.
			# Should return a list with 'b' (coeff vector), 'params' (full vector if diff from b),
			# 'fisher_information' (Hessian matrix), and either 'ssq_b_2' or 'ssq_b_j'
			# for the treatment-effect variance.
			model_output = private$generate_mod(estimate_only = estimate_only) 
			
			private$cached_mod = model_output
			if (is.null(model_output)) {
				private$cache_nonestimable_estimate("model_fit_unavailable")
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			beta_hat_T = as.numeric(model_output$beta_hat_T %||% model_output$b[2])[1L]
			if (!is.finite(beta_hat_T)) {
				private$cache_nonestimable_estimate("model_treatment_estimate_unavailable")
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T = beta_hat_T
			
			if (!is.null(model_output$b)) {
				private$set_fit_warm_start(
					as.numeric(model_output$params %||% model_output$b),
					type = if (!is.null(model_output$params)) "params" else "beta",
					fisher = model_output$fisher_information %||% model_output$XtWX,
					weights = model_output$w %||% model_output$mu,
					force_pd = TRUE
				)
			}
			if (estimate_only) return(invisible(NULL))
			ssq = model_output$ssq_b_2 %||% model_output$ssq_b_j
			ssq = if (length(ssq) >= 1L) as.numeric(ssq)[1L] else NA_real_
			private$cached_values$df = model_output$df %||% NA_real_
			if (is.finite(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
				private$clear_nonestimable_state()
			} else {
				private$cache_nonestimable_se("model_standard_error_unavailable")
			}
		},
		# Helper for subclasses to extract the policy-driven warm start arguments for C++ calls.
		get_backend_warm_start_args = function(expected_length, expected_fisher_dim = expected_length) {
			private$get_optimal_warm_start_config(expected_length, expected_fisher_dim)
		}
	)

InferenceAsympLikStdModCache = R6::R6Class("InferenceAsympLikStdModCache",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = inference_asymp_lik_std_mod_cache_public,
	private = inference_asymp_lik_std_mod_cache_private
)

#' GLM and Kaplan-Meier Inference Without Parametric LR Bootstrap
#'
#' Abstract class parallel to \code{InferenceAsympLikStdModCache} that retains
#' the same asymptotic and likelihood-test behavior while not exposing the
#' parametric LR bootstrap API.
#'
#' @keywords internal
#' @noRd
InferenceAsympLikStdModCacheNoParamBootstrap = R6::R6Class("InferenceAsympLikStdModCacheNoParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = inference_asymp_lik_std_mod_cache_public,
	private = inference_asymp_lik_std_mod_cache_private
)
