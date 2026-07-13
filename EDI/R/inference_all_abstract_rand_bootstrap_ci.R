#' Bootstrap Randomization Confidence Intervals
#'
#' Abstract class implementing confidence intervals by inverting the bootstrap
#' randomization test of \code{InferenceRandBootstrap} over the null effect \code{delta}.
#'
#' @description
#' The CI is the set of \code{delta} values whose bootstrap randomization p-value exceeds
#' \code{alpha}. All p-value evaluations across candidate \code{delta} values share one set
#' of pre-generated draws (resampled row indices and fresh design assignments) — common
#' random numbers — so the p-value is a deterministic, near-monotone function of
#' \code{delta} and the bound search is stable. See \code{InferenceRandBootstrap} for the
#' statistical justification of the test being inverted; the resulting interval inherits
#' its unconditional, superpopulation interpretation and asymptotic validity.
#'
#' @keywords internal
InferenceRandBootstrapCI = R6::R6Class("InferenceRandBootstrapCI",
	lock_objects = FALSE,
	inherit = InferenceRandBootstrap,
	public = list(
		#' @description Computes a confidence interval by inverting the bootstrap randomization
		#'   test over the null effect \code{delta}. When the p-value does not drop below
		#'   \code{alpha / 2} anywhere within the search radius, a conservative bound is returned
		#'   at the search boundary rather than \code{NA}; each such event emits a
		#'   \code{message()} and increments the private field
		#'   \code{rand_bootstrap_ci_conservative_count}.
		#'
		#' @param alpha  				The confidence level 1 - \code{alpha}. Default 0.05.
		#' @param B  					Number of bootstrap randomization draws. Default 501.
		#' @param pval_epsilon  		Bisection tolerance (on both the \code{delta} bracket width
		#'   and the p-value span). Default 0.005.
		#' @param show_progress  		A flag indicating whether progress should be displayed.
		#' @param max_expansions 		Maximum number of bound-doubling expansions when the seed
		#'   interval does not bracket the target p-value. Default 7.
		#' @param bootstrap_type 		Optional bootstrap-resampling scheme; see
		#'   \code{approximate_bootstrap_distribution_beta_hat_T} for legal values. Default \code{NULL}.
		#' @param zero_one_logit_clamp The clamping amount for exact 0 and 1 values when logging.
		#' @return A bootstrap randomization confidence interval. The interval lives on the
		#'   response-transformation scale used by the test (identity for continuous, logit for
		#'   proportion, log for count and survival). Bounds may be conservative (wider than
		#'   necessary) when the p-value inversion cannot be completed within the search radius.
		compute_rand_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, pval_epsilon = 0.005, show_progress = TRUE, max_expansions = 7L, bootstrap_type = NULL, zero_one_logit_clamp = .Machine$double.eps){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Bootstrap randomization inference")
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
				assertCount(as.integer(max_expansions), positive = TRUE)
				assertLogical(show_progress)
				if (private$des_obj_priv_int$response_type == "incidence" && is.null(private$custom_randomization_statistic_function) && is.null(private[["compiled_cpp_stat_fn"]])) {
					stop("Bootstrap randomization confidence intervals are not supported for incidence.")
				}
			}
			# The Monte Carlo p-value has floor 2/B; if that floor is >= alpha/2 the inversion
			# can never bracket and both bounds degrade to the conservative search boundary.
			if (2 / as.integer(B) >= alpha / 2) {
				message(sprintf(
					"compute_rand_bootstrap_confidence_interval: B = %d is too small for alpha = %g (p-value floor 2/B = %.4g >= alpha/2 = %.4g); bounds will be conservative. Use B > %d.",
					as.integer(B), alpha, 2 / as.integer(B), alpha / 2, ceiling(4 / alpha)))
			}
			transform_arg = switch(
				private$des_obj_priv_int$response_type,
				continuous = "none",
				proportion = "logit",
				count = "log",
				survival = "log",
				"none"
			)
			missing_ci = function(reason){
				private$cache_nonestimable_se(reason)
				ci = c(NA_real_, NA_real_)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				ci
			}
			# Common random numbers: one set of (row indices, fresh design assignment) draws
			# shared across every delta evaluation, so the p-value varies smoothly in delta.
			if (!is.null(private$seed)) {
				had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				on.exit(
					if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv),
					add = TRUE
				)
				set.seed(private$seed)
			}
			draws = private$generate_rand_bootstrap_draws(B, bootstrap_type = bootstrap_type, materialize_w = TRUE)
			ci_pval_cache = new.env(parent = emptyenv())
			evaluate_pval = function(delta) {
				private$compute_rand_bootstrap_ci_pval_cached(
					delta = delta,
					B = B,
					transform_arg = transform_arg,
					draws = draws,
					resolution = pval_epsilon,
					ci_pval_cache = ci_pval_cache,
					zero_one_logit_clamp = zero_one_logit_clamp
				)
			}
			est = as.numeric(tryCatch(self$compute_estimate(), error = function(e) NA_real_))
			est = if (length(est) == 0L) NA_real_ else est[1]
			# Null distribution at delta = 0 doubles as the scale estimate for the bound search.
			t0s = self$approximate_rand_bootstrap_distribution_beta_hat_T(
				B = B,
				delta = 0,
				transform_responses = transform_arg,
				show_progress = FALSE,
				rand_bootstrap_draws = draws,
				zero_one_logit_clamp = zero_one_logit_clamp
			)
			t0s_finite = t0s[is.finite(t0s)]
			if (length(t0s_finite) < 2L) {
				return(missing_ci("rand_bootstrap_ci_too_few_finite_null_draws"))
			}
			if (!is.finite(est)) est = stats::median(t0s_finite)
			se_guess = stats::sd(t0s_finite)
			response_scale = stats::sd(private$y, na.rm = TRUE)
			if (!is.finite(response_scale) || response_scale <= 0) response_scale = stats::IQR(private$y, na.rm = TRUE) / 1.349
			if (!is.finite(response_scale) || response_scale <= 0) response_scale = 1
			if (!is.finite(se_guess) || se_guess <= 0) se_guess = response_scale / sqrt(max(1, private$n))
			max_radius = max(25 * se_guess, 6 * response_scale, 1)
			z_target = stats::qnorm(1 - alpha / 2)
			l = private$expand_rand_bootstrap_bound(est - z_target * se_guess, est, alpha / 2, TRUE, max_radius, as.integer(max_expansions), evaluate_pval)
			u = private$expand_rand_bootstrap_bound(est + z_target * se_guess, est, alpha / 2, FALSE, max_radius, as.integer(max_expansions), evaluate_pval)
			if (!all(is.finite(c(l, u)))) {
				return(missing_ci("rand_bootstrap_ci_search_bounds_failed"))
			}
			ci = c(
				private$invert_rand_bootstrap_test_bisection(l, est, alpha / 2, pval_epsilon, TRUE, show_progress, evaluate_pval),
				private$invert_rand_bootstrap_test_bisection(est, u, alpha / 2, pval_epsilon, FALSE, show_progress, evaluate_pval)
			)
			if (length(ci) < 2L || !all(is.finite(ci[1:2]))) {
				private$cache_nonestimable_se("rand_bootstrap_ci_bisection_failed")
			}
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		rand_bootstrap_ci_conservative_count = 0L,
		compute_rand_bootstrap_ci_pval_cached = function(delta, B, transform_arg, draws, resolution, ci_pval_cache, zero_one_logit_clamp = .Machine$double.eps){
			cache_key = private$normalize_delta_for_cache(delta, resolution)
			if (!is.null(ci_pval_cache[[cache_key]])) return(ci_pval_cache[[cache_key]])
			pval = tryCatch(
				as.numeric(self$compute_rand_bootstrap_two_sided_pval(
					B = B,
					delta = delta,
					transform_responses = transform_arg,
					na.rm = TRUE,
					show_progress = FALSE,
					rand_bootstrap_draws = draws,
					zero_one_logit_clamp = zero_one_logit_clamp
				)),
				error = function(e) NA_real_
			)
			pval = if (length(pval) == 0L) NA_real_ else pval[1]
			ci_pval_cache[[cache_key]] = pval
			pval
		},
		expand_rand_bootstrap_bound = function(bound, est, target_pval, lower, max_radius, max_expansions, evaluate_pval){
			if (!is.finite(bound) || !is.finite(est) || !is.finite(max_radius) || max_radius <= 0) return(NA_real_)
			bound = if (lower) max(bound, est - max_radius) else min(bound, est + max_radius)
			pval_bound = evaluate_pval(bound)
			if (is.finite(pval_bound) && pval_bound < target_pval) return(bound)
			step = abs(est - bound)
			if (!is.finite(step) || step <= 0) step = min(max_radius / 4, 1)
			for (iter in seq_len(max_expansions)) {
				step = min(step * 2, max_radius)
				candidate = if (lower) est - step else est + step
				pval_candidate = evaluate_pval(candidate)
				if (is.finite(pval_candidate) && pval_candidate < target_pval) return(candidate)
				if (step >= max_radius) break
			}
			# Conservative fallback: no bracket found within max_radius; return the boundary so
			# the bisection can detect and report a conservative bound.
			if (lower) est - max_radius else est + max_radius
		},
		invert_rand_bootstrap_test_bisection = function(l, u, pval_th, tol, lower, show_progress, evaluate_pval){
			pval_l = evaluate_pval(l)
			pval_u = evaluate_pval(u)
			if (length(pval_l) == 0) pval_l = NA_real_; if (length(pval_u) == 0) pval_u = NA_real_
			for (k in seq_len(30L)) {
				if (!is.na(pval_l) && !is.na(pval_u)) break
				if (is.na(pval_l)) { l = (l + u) / 2; pval_l = evaluate_pval(l) }
				if (is.na(pval_u)) { u = (l + u) / 2; pval_u = evaluate_pval(u) }
			}
			if (is.na(pval_l) || is.na(pval_u) || !all(is.finite(c(l, u)))) return(NA_real_)
			# Conservative bound: the p-value at the search boundary is still >= alpha/2, so the
			# true CI bound lies beyond the search radius. Return the boundary as a valid
			# (conservative) one-sided limit.
			if (lower && is.finite(pval_l) && pval_l >= pval_th) {
				message(sprintf(
					"Bootstrap randomization CI lower bound is conservative: p-value at search boundary delta=%.4g is %.4g >= %.4g. True CI lower bound may extend further left. (rand_bootstrap_ci_conservative_count++)",
					l, pval_l, pval_th))
				private[["rand_bootstrap_ci_conservative_count"]] = private[["rand_bootstrap_ci_conservative_count"]] + 1L
				return(l)
			}
			if (!lower && is.finite(pval_u) && pval_u >= pval_th) {
				message(sprintf(
					"Bootstrap randomization CI upper bound is conservative: p-value at search boundary delta=%.4g is %.4g >= %.4g. True CI upper bound may extend further right. (rand_bootstrap_ci_conservative_count++)",
					u, pval_u, pval_th))
				private[["rand_bootstrap_ci_conservative_count"]] = private[["rand_bootstrap_ci_conservative_count"]] + 1L
				return(u)
			}
			iter = 0; progress_label = if (lower) "BRT CI lower" else "BRT CI upper"
			repeat {
				pval_span = abs(pval_u - pval_l)
				if ((abs(u - l)) <= tol || pval_span <= tol) {
					if (isTRUE(show_progress)) cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g) done\n", progress_label, iter, pval_span, tol))
					return(if (lower) l else u)
				}
				m = (l + u) / 2.0; pval_m = evaluate_pval(m)
				if (is.na(pval_m)) { if (lower) { l = m; pval_l = 0 } else { u = m; pval_u = 0 }; iter = iter + 1; next }
				if (pval_m >= pval_th && lower) { u = m; pval_u = pval_m }
				else if (pval_m >= pval_th && !lower) { l = m; pval_l = pval_m }
				else if (lower) { l = m; pval_l = pval_m }
				else { u = m; pval_u = pval_m }
				iter = iter + 1
				if (isTRUE(show_progress)) cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g)", progress_label, iter, pval_span, tol))
			}
		}
	)
)
