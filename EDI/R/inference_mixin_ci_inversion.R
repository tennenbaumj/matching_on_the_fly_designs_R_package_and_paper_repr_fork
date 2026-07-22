#' Mixin for Likelihood-Test Confidence-Interval Inversion
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' bundling \code{InferenceAsympLik}'s three CI-inversion engines --
#' \code{invert_test_pval_confidence_interval()} (score / lik-ratio / Bartlett
#' pval inversion via root-finding on a supplied \code{pval_fn}),
#' \code{invert_gradient_ci_uniroot()} (gradient-test CI via the same
#' root-finding machinery, with an extra "reject at the estimate" short
#' circuit), and \code{invert_lik_ratio_ci_newton()} (likelihood-ratio CI via
#' the dedicated Newton-Raphson C++ inverter). All three seed from the Wald CI,
#' invert a per-testing-type p-value function, and then validate/clamp the
#' result against the same fallback-to-Wald-or-give-up policy -- that shared
#' policy was previously copy-pasted three times (differing only in the
#' nonestimable-reason string) and now lives once, in
#' \code{finalize_inverted_ci()}.
#'
#' Splice into \code{InferenceAsympLik} via
#' \code{private = c(InferenceMixinCIInversion$private, list(...))}.
#' All 13 direct \code{InferenceAsympLik} subclasses (and everything below
#' them) get these through ordinary R6 private-environment inheritance.
#'
#' @keywords internal
#' @noRd
InferenceMixinCIInversion = list(
	public = list(),
	private = list(
		# Shared "accept, fall back to Wald, or give up" policy for an inverted CI.
		# `ci_vals` is the raw two-element (possibly empty/non-finite) root-finder
		# output; `wald_ci` is the already-computed Wald CI used as a fallback when
		# the inverted CI is missing, out of order, doesn't bracket `est`, or is
		# absurdly wide. `unavailable_reason` is the nonestimable-SE reason cached
		# when neither the inverted CI nor the Wald fallback is usable.
		finalize_inverted_ci = function(ci_vals, alpha, est, wald_ci, unavailable_reason){
			ci = if (length(ci_vals) == 2L) sort(as.numeric(ci_vals[1:2]), na.last = TRUE) else c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			if (length(ci) < 2L || !all(is.finite(ci[1:2])) || ci[1L] > est || ci[2L] < est || any(abs(ci[1:2]) > private$likelihood_ci_max_abs)) {
				fallback = sort(as.numeric(wald_ci[1:2]), na.last = TRUE)
				if (length(fallback) >= 2L && all(is.finite(fallback)) && fallback[1L] <= est && fallback[2L] >= est && all(abs(fallback[1:2]) <= private$likelihood_ci_max_abs)) {
					ci = fallback
				} else {
					private$cache_nonestimable_se(unavailable_reason)
					ci = c(NA_real_, NA_real_)
				}
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			}
			ci
		},
		invert_test_pval_confidence_interval = function(alpha, testing_type = private$testing_type, bartlett_B = NULL){
			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			se = private$get_standard_error()
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = private$compute_wald_confidence_interval_impl(alpha)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else NA_real_
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else NA_real_

			testing_type = private$normalize_testing_type(testing_type)
			spec = private$get_likelihood_test_spec()
			bartlett_types = c("lik_ratio_bartlett_approx", "lik_ratio_bartlett_exact")
			pval_fn = if (!is.null(spec) && testing_type %in% c("score", "gradient", "lik_ratio", bartlett_types)) {
				function(delta) {
					private$get_memoized_likelihood_test_pval(
						delta = delta,
						testing_type = testing_type,
						spec = spec,
						warm_cache_key = if (identical(testing_type, "score")) {
							"likelihood_test:score"
						} else {
							paste0(testing_type, "_ci")
						},
						bartlett_B = bartlett_B
					)
				}
			} else {
				function(delta) self$compute_asymp_two_sided_pval(delta)
			}

			# A likelihood-backed CI cannot be inverted when the test is unavailable
			# at the fitted estimate.  In particular, some Cox score paths expose
			# an explicit non-estimable standard error and otherwise pass NA values
			# into the C++ root search, which may return an empty endpoint vector.
			p_est = tryCatch(pval_fn(est), error = function(e) NA_real_)
			if (!is.finite(p_est)) {
				private$cache_nonestimable_se(paste0(testing_type, "_test_unavailable"))
				return(c(NA_real_, NA_real_))
			}

			ci_vals = pval_invert_ci_cpp(
				pval_fn    = pval_fn,
				est        = est,
				alpha      = alpha,
				step       = step,
				lower_seed = lower_seed,
				upper_seed = upper_seed
			)

			private$finalize_inverted_ci(ci_vals, alpha, est, wald_ci, "test_inversion_confidence_interval_unavailable")
		},

		invert_gradient_ci_uniroot = function(alpha){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) return(private$invert_test_pval_confidence_interval(alpha))

			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			j = as.integer(spec$j)
			if (length(j) != 1L || !is.finite(j) || j < 1L) return(c(NA_real_, NA_real_))

			se = private$get_standard_error()
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = private$compute_wald_confidence_interval_impl(alpha)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else est - step
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else est + step

			pval_fn = function(delta){
				private$get_memoized_likelihood_test_pval(
					delta = delta,
					testing_type = "gradient",
					spec = spec,
					warm_cache_key = "gradient_ci"
				)
			}

			p_est = tryCatch(pval_fn(est), error = function(e) NA_real_)
			if (!is.finite(p_est)) {
				private$cache_nonestimable_se("gradient_test_unavailable")
				return(c(NA_real_, NA_real_))
			}
			if (p_est < alpha) {
				ci = c(est, est)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}

			ci_vals = tryCatch(
				pval_invert_ci_cpp(
					pval_fn    = pval_fn,
					est        = est,
					alpha      = alpha,
					step       = step,
					lower_seed = lower_seed,
					upper_seed = upper_seed
				),
				error = function(e) {
					private$cache_nonestimable_se("gradient_ci_inversion_failed")
					numeric(0)
				}
			)
			private$finalize_inverted_ci(ci_vals, alpha, est, wald_ci, "gradient_confidence_interval_unavailable")
		},

		invert_lik_ratio_ci_newton = function(alpha){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) return(private$invert_test_pval_confidence_interval(alpha))

			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			full_eval = tryCatch(private$get_memoized_likelihood_test_eval(
				delta = est,
				testing_type = "lik_ratio",
				spec = spec,
				warm_cache_key = "lik_ratio_ci",
				include_null_fit = FALSE,
				include_full_negloglik = TRUE
			), error = function(e) list(full_negloglik = NA_real_))
			full_negloglik = full_eval$full_negloglik %||% NA_real_
			if (!is.finite(full_negloglik)) return(private$invert_test_pval_confidence_interval(alpha))

			se = private$get_standard_error()
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = private$compute_wald_confidence_interval_impl(alpha)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else est - step
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else est + step

			j = as.integer(spec$j)
				fit_null_fn = function(delta){
					eval = tryCatch(
						private$get_memoized_likelihood_test_eval(
							delta = delta,
							testing_type = "lik_ratio",
							spec = spec,
							warm_cache_key = "lik_ratio_ci",
							include_score = FALSE,
							include_full_negloglik = FALSE,
							include_null_negloglik = FALSE
						),
						error = function(e) NULL
					)
					if (is.null(eval) || isTRUE(eval$invalid) || is.null(eval$null_fit)) return(NULL)
					eval$null_fit
				}
				neg_loglik_fn = function(fit){
					delta_fit = attr(fit, "edi_likelihood_test_delta", exact = TRUE)
					eval = if (is.finite(delta_fit %||% NA_real_)) {
						tryCatch(
							private$get_memoized_likelihood_test_eval(
								delta = delta_fit,
								testing_type = "lik_ratio",
								spec = spec,
								warm_cache_key = "lik_ratio_ci",
								include_score = FALSE,
								include_full_negloglik = FALSE,
								include_null_negloglik = TRUE
							),
							error = function(e) NULL
						)
					} else NULL
					if (!is.null(eval) && is.finite(eval$null_negloglik %||% NA_real_)) return(eval$null_negloglik)
					tryCatch(spec$neg_loglik(fit), error = function(e) NA_real_)
				}
				score_fn = function(fit){
					delta_fit = attr(fit, "edi_likelihood_test_delta", exact = TRUE)
					eval = if (is.finite(delta_fit %||% NA_real_)) {
						tryCatch(
							private$get_memoized_likelihood_test_eval(
								delta = delta_fit,
								testing_type = "lik_ratio",
								spec = spec,
								warm_cache_key = "lik_ratio_ci",
								include_score = TRUE,
								include_full_negloglik = FALSE,
								include_null_negloglik = FALSE
							),
							error = function(e) NULL
						)
					} else NULL
					if (!is.null(eval) && !is.null(eval$score)) return(eval$score)
					tryCatch(spec$score(fit), error = function(e) NULL)
				}

			ci_vals = tryCatch(lrt_ci_nr_cpp(
				fit_null_fn    = fit_null_fn,
				neg_loglik_fn  = neg_loglik_fn,
				score_fn       = score_fn,
				est            = est,
				full_negloglik = full_negloglik,
				alpha          = alpha,
				step           = step,
				lower_seed     = lower_seed,
				upper_seed     = upper_seed,
				j              = j
			), error = function(e) numeric(0))

			private$finalize_inverted_ci(ci_vals, alpha, est, wald_ci, "lik_ratio_confidence_interval_unavailable")
		}
	)
)
