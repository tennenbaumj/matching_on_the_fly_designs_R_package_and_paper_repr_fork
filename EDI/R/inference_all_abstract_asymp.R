#' Asymptotic Inference
#'
#' Abstract class for asymptotic inference.
#'
#' @keywords internal
InferenceAsymp = R6::R6Class("InferenceAsymp",
	lock_objects = FALSE,
	inherit = InferenceJackknife,
	public = list(
		#' @description Computes an asymptotic confidence interval using the configured test.
		#'
		#' @param alpha  				Significance level 1 - \code{alpha}. Default 0.05.
		#'
		#' @return 	A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_wald_confidence_interval_impl(alpha)
		},
		#' @description Computes an asymptotic two-sided p-value using the configured test.
		#'
		#' @param delta  				Null treatment effect to test against. Default 0.
		#'
		#' @return 	The asymptotic p-value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$compute_wald_two_sided_pval_impl(delta)
		},
		#' @description Gets the asymptotic testing methods supported by this inference object.
		get_supported_testing_types = function(){
			private$get_supported_testing_types_impl()
		},
		#' @description Computes the Wald two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_wald_two_sided_pval = function(delta = 0){
			private$compute_wald_two_sided_pval_impl(delta)
		},
		#' @description Computes the Wald confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_wald_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_wald_confidence_interval_impl(alpha)
		},
		#' @description Abstract method to compute the treatment estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		#' @return 	A scalar treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			stop("Must be implemented by concrete class.")
		},
		#' @description Returns the model object from the last call that produced the treatment
		#' estimate and SE. Calls \code{compute_estimate()} first if needed.
		#'
		#' @return The cached model object (type depends on the concrete class).
		get_mod = function(){
			if (is.null(private$cached_mod)) self$compute_estimate()
			private$cached_mod
		},
			#' @description Prints a summary of the model from the last call that produced the
			#' treatment estimate and SE.
			get_summary = function(){
				if (is.null(private$cached_mod) || is.null(private$cached_values$summary_table)) {
					tryCatch(
						self$compute_estimate(estimate_only = FALSE),
						error = function(e) tryCatch(self$compute_estimate(), error = function(e2) NULL)
					)
				}
				mod = self$get_mod()
				if (is.null(mod)) {
					cat("No model available (call compute_estimate() first).\n")
					return(invisible(NULL))
				}
				fallback = private$cached_values$model_fit_fallback
				if (is.list(fallback) && isTRUE(fallback$used)) {
					cat("Model fit fallback: ", fallback$fitted_model, "\n", sep = "")
					if (!is.null(fallback$requested_model)) {
						cat("Requested model: ", fallback$requested_model, "\n", sep = "")
					}
					if (!is.null(fallback$reason)) {
						cat("Fallback reason: ", fallback$reason, "\n", sep = "")
					}
					if (length(fallback$omitted_conditional %||% character(0))) {
						cat("Omitted conditional covariates: ", paste(fallback$omitted_conditional, collapse = ", "), "\n", sep = "")
					}
					if (length(fallback$omitted_auxiliary %||% character(0))) {
						cat("Omitted auxiliary covariates: ", paste(fallback$omitted_auxiliary, collapse = ", "), "\n", sep = "")
					}
					cat("Omitted coefficient rows are shown as NA because those terms were not included in the fallback fit.\n")
				}
				if (identical(class(mod), "list")) {
					if (!is.null(private$cached_values$summary_table)) {
						print(private$cached_values$summary_table)
				} else {
					print(mod)
				}
			} else {
				print(summary(mod))
			}
			invisible(NULL)
		}
	),
	private = list(
		cached_mod = NULL,
		get_standard_error = function(){
			se = private$cached_values$s_beta_hat_T
			if (!is.null(se) && is.finite(se) && se > 0) return(se)
			stop(class(self)[1], " must implement get_standard_error() to support Wald-type inference.")
		},
		get_degrees_of_freedom = function(){
			df = private$cached_values$df
			if (!is.null(df)) return(df)
			NA_real_
		},
		get_supported_testing_types_impl = function(){
			"wald"
		},
		compute_wald_confidence_interval_impl = function(alpha){
			est = self$compute_estimate()
			se = private$get_standard_error()
			df = private$get_degrees_of_freedom()
			
			if (length(est) != 1L || length(se) != 1L) return(c(NA_real_, NA_real_))
			if (!is.finite(est) || !is.finite(se) || se <= 0) return(c(NA_real_, NA_real_))
			
			critical_val = if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else stats::qnorm(1 - alpha / 2)
			
			ci = c(est - critical_val * se, est + critical_val * se)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		compute_wald_two_sided_pval_impl = function(delta){
			est = self$compute_estimate()
			se = private$get_standard_error()
			df = private$get_degrees_of_freedom()
			
			if (length(est) != 1L || length(se) != 1L) return(NA_real_)
			if (!is.finite(est) || !is.finite(se) || se <= 0) return(NA_real_)
			
			t_stat = (est - delta) / se
			
			if (is.finite(df)) {
				2 * stats::pt(-abs(t_stat), df = df)
			} else {
				2 * stats::pnorm(-abs(t_stat))
			}
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},
		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},
		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},
		
		# Shared helpers for normal/t tests
		compute_z_or_t_ci_from_s_and_df = function(alpha){
			beta_hat_T = private$cached_values$beta_hat_T
			s_beta_hat_T = private$cached_values$s_beta_hat_T
			df = private$cached_values$df
			if (is.null(df)) df = NA_real_
			
			if (length(beta_hat_T) != 1L || length(s_beta_hat_T) != 1L) return(c(NA_real_, NA_real_))
			if (!is.finite(beta_hat_T) || !is.finite(s_beta_hat_T) || s_beta_hat_T <= 0) return(c(NA_real_, NA_real_))
			
			mult = if (!is.finite(df)) stats::qnorm(1 - alpha / 2) else stats::qt(1 - alpha / 2, df = df)
			ci = c(beta_hat_T - mult * s_beta_hat_T, beta_hat_T + mult * s_beta_hat_T)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		compute_z_or_t_two_sided_pval_from_s_and_df = function(delta){
			beta_hat_T = private$cached_values$beta_hat_T
			s_beta_hat_T = private$cached_values$s_beta_hat_T
			df = private$cached_values$df
			if (is.null(df)) df = NA_real_
			
			if (length(beta_hat_T) != 1L || length(s_beta_hat_T) != 1L) return(NA_real_)
			if (!is.finite(beta_hat_T) || !is.finite(s_beta_hat_T) || s_beta_hat_T <= 0) return(NA_real_)
			
			val = (beta_hat_T - delta) / s_beta_hat_T
			if (!is.finite(df)) 2 * stats::pnorm(-abs(val)) else 2 * stats::pt(-abs(val), df = df)
		},
		invert_ci_to_find_two_sided_pval_for_treatment_effect = function(delta = 0){
			# Use bisection to find alpha such that delta is on the boundary of the CI.
			# The p-value is the largest alpha such that delta is outside the (1-alpha) CI.
			
			f = function(alpha) {
				ci = self$compute_asymp_confidence_interval(alpha = alpha)
				if (any(!is.finite(ci))) return(NA_real_)
				# Distance to the nearest boundary. If delta is inside, return positive.
				# If delta is outside, return negative.
				# Actually, easier: return 0 if delta is on boundary.
				min(abs(ci - delta))
			}
			
			# More robust: check if delta is within the range of the estimate
			est = self$compute_estimate()
			if (!is.finite(est)) return(NA_real_)
			
			# We want to find alpha such that ci_boundary(alpha) == delta
			# This is monotonic in alpha.
			
			target_fn = function(alpha) {
				ci = self$compute_asymp_confidence_interval(alpha = alpha)
				if (est > delta) {
					# Null is below estimate, we care about the lower bound
					ci[1] - delta
				} else {
					# Null is above estimate, we care about the upper bound
					ci[2] - delta
				}
			}
			
			# p-value is usually between 0 and 1.
			lower = .Machine$double.eps
			upper = 1 - .Machine$double.eps
			
			# Check if target_fn has different signs at boundaries
			tl = tryCatch(target_fn(lower), error = function(e) NA_real_)
			tu = tryCatch(target_fn(upper), error = function(e) NA_real_)
			
			if (!is.finite(tl) || !is.finite(tu)) return(NA_real_)
			if (tl * tu > 0) {
				# If both are same sign, delta is either always outside or always inside
				if (abs(est - delta) < .Machine$double.eps) return(1)
				if (abs(tl) < abs(tu)) return(lower) else return(upper)
			}
			
			stats::uniroot(target_fn, lower = lower, upper = upper, tol = 1e-6)$root
		},
		# Complexity tiers for optimization paths:
		# "heavy": Expensive (quadrature, transcendental). Always use full warm start (Beta + Info).
		# "medium": Moderate (standard IRLS). Use Beta + Weights.
		# "light": Extremely fast. Use Beta-only to avoid R/C++ bridge overhead for matrices.
		get_complexity_tier = function() "light",
		
		# Helper to extract optimal warm start components based on complexity tier and problem scale.
		# Returns a list of arguments for the backend solver.
		get_optimal_warm_start_config = function(expected_length, expected_fisher_dim = expected_length) {
			tier = private$get_complexity_tier()
			p = expected_length
			
			start_beta = private$get_fit_warm_start_for_length("beta", p)
			start_params = private$get_fit_warm_start_for_length("params", p)
			coef_start = start_beta %||% start_params
			
			# Tier 3 (Light) or very high dimensionality: Beta-only is safest and fastest.
			# If p > 50, matrix copying overhead for Hessian usually negates its benefit.
			if (tier == "light" || p > 50) {
				return(list(
					start_beta = coef_start,
					warm_start_beta = coef_start,
					start_params = start_params
				))
			}

			if (tier == "medium" && identical(private$active_resampling_operation, "bayesian_boot")) {
				return(list(
					start_beta = coef_start,
					warm_start_beta = coef_start,
					start_params = start_params
				))
			}
			
			# Tier 2 (Medium): Beta + Weights
			if (tier == "medium") {
				return(list(
					start_beta = coef_start,
					warm_start_beta = coef_start,
					start_params = start_params,
					warm_start_weights = private$get_fit_warm_start_weights(private$n)
				))
			}
			
			# Tier 1 (Heavy): Full Warm Start (Beta + Info)
			return(list(
				start_beta = coef_start,
				warm_start_beta = coef_start,
				start_params = start_params,
				warm_start_weights = private$get_fit_warm_start_weights(private$n),
				warm_start_fisher_info = private$get_fit_warm_start_fisher(expected_fisher_dim)
			))
		}
	)
)
