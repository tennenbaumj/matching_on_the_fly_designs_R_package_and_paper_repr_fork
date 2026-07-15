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
#' \code{delta} and the bound search is stable. See \code{\link{InferenceRandBootstrap}}
#' for the statistical justification of the test being inverted; the resulting interval
#' inherits its unconditional, superpopulation interpretation and asymptotic validity.
#'
#' Users do not instantiate this class directly: every concrete inference class in the
#' package inherits from it, so \code{compute_rand_bootstrap_confidence_interval} is
#' available on any inference object.
#'
#' @examples
#' \dontrun{
#' seq_des = DesignSeqOneByOneKK14$new(n = 100, response_type = "continuous")
#' # ... run the experiment: add subjects and responses ...
#' seq_des_inf = InferenceAllSimpleMeanDiff$new(seq_des)
#' seq_des_inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 501)
#' }
#' @keywords internal
InferenceRandBootstrapCI = R6::R6Class("InferenceRandBootstrapCI",
	lock_objects = FALSE,
	inherit = InferenceRandBootstrap,
	public = list(
		#' @description Computes a confidence interval by inverting the bootstrap randomization
		#'   test over the null effect \code{delta}. For statistics that are affine in the
		#'   additive sharp-null shift (e.g. the simple mean difference and the OLS treatment
		#'   coefficient), the inversion is performed in closed form from the breakpoints of the
		#'   p-value step function — exact given the draws and requiring no bisection. Otherwise,
		#'   the generic bisection search is used. When the p-value does not drop below
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
		#' @param type  				CI type. \code{"percentile"} (default) inverts the raw BRT
		#'   p-value via bisection (or closed-form for affine statistics).
		#'   \code{"studentized"} inverts the signed-pivot BRT p-value
		#'   \eqn{p(\delta) = 2\min(P(z^0_b \ge z), P(z^0_b \le z))} where
		#'   \eqn{z^0_b = (t^0_b(\delta) - \delta)/\hat{s}^0_b}; yields asymmetric CI.
		#'   \code{"symmetric-percentile-t"} inverts the absolute-pivot version
		#'   \eqn{p(\delta) = P(|t^0_b - \delta|/\hat{s}^0_b \ge |t - \delta|/\hat{s})};
		#'   yields a CI symmetric around the observed estimate.
		#'   Both SE-based types pre-compute \eqn{\hat{s}^0_b} once at \eqn{\delta = 0} and
		#'   reuse it across bisection steps. Yield O(\eqn{n^{-1}}) coverage error versus
		#'   O(\eqn{n^{-1/2}}) for \code{"percentile"} when the pivot is asymptotically normal.
		#'   Require the class to expose a standard error; return \code{NA} bounds in harden
		#'   mode if unavailable.
		#'   \code{"smoothed"} adds per-draw kernel noise \eqn{\varepsilon_b \sim N(0, \hat{\sigma}/\sqrt{n})}
		#'   to the resampled responses before imposing the null shift, reducing discreteness.
		#'   Only meaningful for continuous responses.
		#' @return A bootstrap randomization confidence interval. The interval lives on the
		#'   response-transformation scale used by the test (identity for continuous, logit for
		#'   proportion, log for count and survival). Bounds may be conservative (wider than
		#'   necessary) when the p-value inversion cannot be completed within the search radius.
		compute_rand_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, pval_epsilon = 0.005, show_progress = TRUE, max_expansions = 7L, bootstrap_type = NULL, zero_one_logit_clamp = .Machine$double.eps, type = "percentile"){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Bootstrap randomization inference")
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
				assertCount(as.integer(max_expansions), positive = TRUE)
				assertLogical(show_progress)
				assertChoice(tolower(type), c("percentile", "studentized", "symmetric-percentile-t", "smoothed"))
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
			# Enable SMC early-stopping during the CI bisection: each evaluate_pval(delta) call
			# on the pre-materialized draws will stop once the Clopper-Pearson band clears alpha/2,
			# saving the bulk of the B worker fits for GLM-class CIs that cannot use a C++ kernel.
			{
				brt_mc_b = as.integer(B)
				brt_default_batch = min(brt_mc_b, max(25L, as.integer(ceiling(2 * sqrt(brt_mc_b)))))
				old_brt_mc_ctrl = private$brt_mc_control
				private$brt_mc_control = list(
					mc_enable = brt_mc_b >= 200L,
					mc_batch_size = brt_default_batch,
					mc_min_draws = min(brt_mc_b, max(100L, 2L * brt_default_batch)),
					mc_conf_level = 0.99,
					mc_stop_threshold = alpha / 2
				)
				on.exit({ private$brt_mc_control = old_brt_mc_ctrl }, add = TRUE)
			}
			# Closed-form shortcut: when the statistic is affine in delta under the additive
			# sharp-null shift (mean difference, OLS treatment coefficient), each null draw is
			# t0_b(delta) = A_b + delta * c_b, so the p-value is a step function of delta whose
			# breakpoints are (t - A_b) / c_b and the CI is read off exactly — no bisection.
			# Only valid for the raw (percentile) test — skipped for all other types.
			if (identical(tolower(type), "percentile") && identical(transform_arg, "none") &&
				is.null(private$custom_randomization_statistic_function) &&
				is.null(private[["compiled_cpp_stat_fn"]]) &&
				private$has_private_method("compute_rand_bootstrap_ci_affine_coefs")) {
				ci_closed = tryCatch({
					affine = private$compute_rand_bootstrap_ci_affine_coefs(draws)
					if (is.null(affine)) NULL else {
						t_obs = private$compute_treatment_estimate_during_randomization_inference()
						if (is.list(t_obs) && "b" %in% names(t_obs)) t_obs = t_obs$b[1]
						t_obs = as.numeric(t_obs)[1]
						if (!is.finite(t_obs)) NULL else {
							private$closed_form_ci_from_affine_null_draws(affine$A, affine$c, t_obs, alpha)
						}
					}
				}, error = function(e) NULL)
				if (!is.null(ci_closed)) {
					if (length(ci_closed) != 2L) {
						private$cache_nonestimable_se("rand_bootstrap_ci_closed_form_unavailable")
						ci_closed = c(NA_real_, NA_real_)
					}
					names(ci_closed) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
					return(ci_closed)
				}
			}
			type_lc = tolower(type)
			# Studentized / symmetric-percentile-t BRT CI: pivot each null draw by its per-draw SE.
			# SE0_b is pre-computed once at delta = 0 and held fixed across all bisection steps.
			# "studentized" uses signed pivots (asymmetric CI); "symmetric-percentile-t" uses |pivot|.
			if (type_lc %in% c("studentized", "symmetric-percentile-t")) {
				symmetric_pivot = identical(type_lc, "symmetric-percentile-t")
				se_obs = tryCatch(private$infer_original_se(), error = function(e) NA_real_)
				if (!is.finite(se_obs) || se_obs <= 0) {
					return(missing_ci("rand_bootstrap_ci_studentized_se_unavailable"))
				}
				y0_full_0 = private$y  # sharp null at delta = 0: y0 = y_obs
				brt_stats_0 = private$compute_brt_null_statistics_with_se(
					draws, 0, transform_arg, y0_full_0, zero_one_logit_clamp
				)
				t0s = brt_stats_0$t0
				se0_b = brt_stats_0$se0
				# Pre-populate the BRT distribution cache at delta = 0 from the SE pass so the
				# bisection can reuse these draws without a second round of sub-inference creation.
				{
					draws_id_stud = attr(draws, "draws_id")
					delta_key_0 = formatC(0.0, digits = 17L, format = "fg", flag = "#")
					ck0 = paste(as.integer(B), delta_key_0, transform_arg, if (!is.null(draws_id_stud)) draws_id_stud else "fresh", sep = "|")
					if (is.null(private$cached_values$rand_boot_distr_cache)) private$cached_values$rand_boot_distr_cache = list()
					private$cached_values$rand_boot_distr_cache[[ck0]] = t0s
				}
				est = as.numeric(tryCatch(self$compute_estimate(), error = function(e) NA_real_))
				est = if (length(est) == 0L) NA_real_ else est[1]
				t_obs_stud = est
				stud_pval_cache = new.env(parent = emptyenv())
				evaluate_pval = function(delta) {
					ck = private$normalize_delta_for_cache(delta, pval_epsilon)
					if (!is.null(stud_pval_cache[[ck]])) return(stud_pval_cache[[ck]])
					t0_b_d = self$approximate_rand_bootstrap_distribution_beta_hat_T(
						B = B, delta = delta, transform_responses = transform_arg,
						show_progress = FALSE, rand_bootstrap_draws = draws,
						zero_one_logit_clamp = zero_one_logit_clamp
					)
					pval = private$compute_two_sided_brt_pval_studentized(
						t_obs_stud, t0_b_d, se0_b, delta, se_obs, symmetric = symmetric_pivot
					)
					stud_pval_cache[[ck]] = pval
					pval
				}
			} else if (identical(type_lc, "smoothed")) {
				# Smoothed BRT CI: pre-generate one noise vector per draw (CRN across delta).
				# Tags draws with a "_smoothed" draws_id suffix so the distribution cache
				# is keyed separately from the same draws used without smoothing.
				n_int = as.integer(private$n)
				sd_y = stats::sd(private$y, na.rm = TRUE)
				if (!is.finite(sd_y) || sd_y <= 0) sd_y = 1.0
				bw = sd_y / sqrt(n_int)
				for (b in seq_along(draws)) draws[[b]][["smooth_noise"]] = stats::rnorm(n_int, 0, bw)
				attr(draws, "draws_id") = paste0(attr(draws, "draws_id"), "_smoothed")
				est = as.numeric(tryCatch(self$compute_estimate(), error = function(e) NA_real_))
				est = if (length(est) == 0L) NA_real_ else est[1]
				t0s = self$approximate_rand_bootstrap_distribution_beta_hat_T(
					B = B, delta = 0, transform_responses = transform_arg,
					show_progress = FALSE, rand_bootstrap_draws = draws,
					zero_one_logit_clamp = zero_one_logit_clamp
				)
				smooth_pval_cache = new.env(parent = emptyenv())
				t_obs_sm = est
				evaluate_pval = function(delta) {
					ck = private$normalize_delta_for_cache(delta, pval_epsilon)
					if (!is.null(smooth_pval_cache[[ck]])) return(smooth_pval_cache[[ck]])
					t0s_d = self$approximate_rand_bootstrap_distribution_beta_hat_T(
						B = B, delta = delta, transform_responses = transform_arg,
						show_progress = FALSE, rand_bootstrap_draws = draws,
						zero_one_logit_clamp = zero_one_logit_clamp
					)
					pval = private$compute_two_sided_randomization_pval_from_t0s(t0s_d, t_obs_sm)
					smooth_pval_cache[[ck]] = pval
					pval
				}
			} else {
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
			} # end percentile/studentized branch
			t0s_finite = t0s[is.finite(t0s)]
			if (length(t0s_finite) < 2L) {
				return(missing_ci("rand_bootstrap_ci_too_few_finite_null_draws"))
			}
			if (!is.finite(est)) {
				est = stats::median(t0s_finite)
				if (type_lc %in% c("studentized", "symmetric-percentile-t") && exists("t_obs_stud") && !is.finite(t_obs_stud)) t_obs_stud = est
				if (identical(type_lc, "smoothed") && exists("t_obs_sm") && !is.finite(t_obs_sm)) t_obs_sm = est
			}
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
			if (length(ci) != 2L || !all(is.finite(ci[1:2]))) {
				private$cache_nonestimable_se("rand_bootstrap_ci_bisection_failed")
				ci = c(NA_real_, NA_real_)
			}
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		rand_bootstrap_ci_conservative_count = 0L,
		# Exact inversion of the bootstrap randomization test when every null draw is affine in
		# delta: t0_b(delta) = A[b] + delta * c_slopes[b]. The two-sided Monte Carlo p-value
		# (with its 2/B floor) is piecewise constant with breakpoints at (t_obs - A_b) / c_b;
		# the acceptance region {delta : pval(delta) >= alpha} is delimited by breakpoints, so
		# the CI is found by evaluating the count-based p-value on each open interval between
		# sorted breakpoints. Endpoints are closed (a draw sitting exactly at t_obs counts in
		# both tails, so the p-value at a breakpoint is >= that of its neighboring intervals).
		# Returns c(lower, upper) or NULL when the inversion is not possible (too few finite
		# draws, an unbounded acceptance region, or a degenerate breakpoint set) — callers
		# fall back to the generic bisection.
		closed_form_ci_from_affine_null_draws = function(A, c_slopes, t_obs, alpha){
			# Match the bisection path's threshold convention (inherited from InferenceRandCI):
			# the two-sided p-value is inverted at alpha / 2.
			pval_th = alpha / 2
			ok = is.finite(A) & is.finite(c_slopes)
			A = as.numeric(A[ok]); cs = as.numeric(c_slopes[ok])
			B_valid = length(A)
			# p-value floor 2/B must be able to drop below the threshold for the inversion to bracket
			if (B_valid == 0L || 2 / B_valid >= pval_th) return(NULL)
			# Slope ~ 0 draws never flip; carry their constant tail contributions separately
			is_const = abs(cs) < 1e-12
			n_const_ge = sum(A[is_const] >= t_obs)
			n_const_le = sum(A[is_const] <= t_obs)
			A_v = A[!is_const]; c_v = cs[!is_const]
			if (length(A_v) == 0L) return(NULL)
			breakpoints = sort(unique((t_obs - A_v) / c_v))
			K = length(breakpoints)
			if (K < 2L) return(NULL)
			pval_at = function(d){
				t0 = A_v + d * c_v
				G = sum(t0 >= t_obs) + n_const_ge
				L = sum(t0 <= t_obs) + n_const_le
				min(1, max(2 / B_valid, 2 * min(G, L) / B_valid))
			}
			# Probe each open interval: outer probes detect an unbounded acceptance region
			span_pad = max(1, diff(range(breakpoints)))
			probes = c(
				breakpoints[1] - span_pad,
				(breakpoints[-1] + breakpoints[-K]) / 2,
				breakpoints[K] + span_pad
			)
			accepted = vapply(probes, pval_at, numeric(1)) >= pval_th
			if (accepted[1] || accepted[K + 1L]) return(NULL)
			interior = which(accepted)
			if (length(interior) == 0L) return(NULL)
			# probe j (2..K) covers the interval (breakpoints[j-1], breakpoints[j]); closed
			# breakpoint endpoints bound the acceptance region
			c(breakpoints[interior[1] - 1L], breakpoints[interior[length(interior)]])
		},
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
