#' Randomization-based Confidence Intervals
#'
#' Abstract class for randomization-based confidence interval inference.
#'
#' @keywords internal
InferenceRandCI = R6::R6Class("InferenceRandCI",
	lock_objects = FALSE,
	inherit = InferenceRand,
	public = list(
		#' @description Compute a randomization-based two-sided p-value for the treatment effect.
		#' @param r Number of randomization vectors.
		#' @param delta Null treatment effect value.
		#' @param transform_responses Response transformation to apply during the test.
		#' @param na.rm Whether to remove non-finite simulated statistics.
		#' @param show_progress Whether to show progress.
		#' @param permutations Optional pre-generated assignment draws.
		#' @param type Optional incidence-specific exact randomization type.
		#' @param args_for_type Optional arguments keyed by \code{type}.
		#' @param zero_one_logit_clamp The clamping amount for exact 0 and 1 values when logging
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @return A two-sided p-value.
		compute_rand_two_sided_pval = function(r = 501, delta = 0, transform_responses = "none", na.rm = TRUE, show_progress = TRUE, permutations = NULL, type = NULL, args_for_type = NULL, zero_one_logit_clamp = .Machine$double.eps){
			# message("In InferenceRandCI$compute_rand_two_sided_pval")
			# message("Class: ", class(self)[1])
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Randomization inference")
				assertLogical(na.rm)
			}
			if (should_run_asserts()) {
				if (private$des_obj_priv_int$response_type == "incidence" && is.null(private$custom_randomization_statistic_function)){
					if (!identical(transform_responses, "none")) {
						stop("transform_responses is not supported for incidence randomization inference.")
					}
					rand_type = if (is.null(type)) "Zhang" else type
					exact_args = private$normalize_exact_inference_args(
						rand_type,
						args_for_type = args_for_type
					)
					return(private$compute_exact_two_sided_pval_rand(rand_type, delta, exact_args))
				}
			}
			if (should_run_asserts()) {
				private$assert_no_incidence_only_randomization_args(private$des_obj_priv_int$response_type, type, args_for_type)
			}
			super$compute_rand_two_sided_pval(
				r = r,
				delta = delta,
				transform_responses = transform_responses,
				na.rm = na.rm,
				show_progress = show_progress,
				permutations = permutations,
				zero_one_logit_clamp = zero_one_logit_clamp
			)
		},
		#' @description Computes a randomization-based confidence interval.
		#' @param alpha  				Significance level.
		#' @param r  	Number of randomization vectors.
		#' @param pval_epsilon  		Bisection tolerance.
		#' @param show_progress  	Show progress.
		#' @param type Optional incidence-specific exact randomization type.
		#' @param args_for_type Optional arguments keyed by \code{type}.
		#' @param ci_search_control Optional control list for randomization-CI search. Supported
		#'   entries are \code{fallback}, \code{seed}, \code{max_radius_se_mult} (default 25),
		#'   \code{max_radius_scale_mult} (default 6), \code{max_expansions} (default 7),
		#'   \code{seed_boot_B}, Monte Carlo settings \code{mc_enable}, \code{mc_batch_size},
		#'   \code{mc_min_draws}, and \code{mc_conf_level}, midpoint-cache settings
		#'   \code{pval_cache_enable} and \code{pval_cache_resolution}, and model-fit
		#'   reuse settings \code{fit_warm_start_enable} and
		#'   \code{fit_reuse_factorizations}. Set \code{mc_enable = FALSE} to force
		#'   full enumeration of all requested randomization draws.
		#'   The search radius is \code{max(max_radius_se_mult * se_guess, max_radius_scale_mult * sd(y))}.
		#'   When the randomization p-value does not drop below \code{alpha/2} anywhere within
		#'   the search radius (e.g. the test has low power or the design has few unique
		#'   permutations), a \emph{conservative} CI bound is returned at the search boundary
		#'   rather than \code{NA}. This guarantees a valid (though possibly wide) interval.
		#'   Each such event emits a \code{message()} and increments the private field
		#'   \code{rand_ci_conservative_count} for monitoring.
		#' @return Randomization CI. Bounds may be conservative (wider than necessary) when the
		#'   p-value inversion cannot be completed within the search radius; see
		#'   \code{ci_search_control} for details.
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, type = NULL, args_for_type = NULL, ci_search_control = NULL){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Randomization inference")
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(r, positive = TRUE); assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
				assertLogical(show_progress); show_progress = isTRUE(show_progress) && self$num_cores == 1
			}
			ci_search_control = private$normalize_randomization_ci_search_control(ci_search_control, r, pval_epsilon)
			ci_search_control$mc_stop_threshold = alpha / 2
			dispatch_cores = private$effective_parallel_cores("rand_ci", self$num_cores)
			if (dispatch_cores != self$num_cores) {
				serial_inf = self$duplicate(verbose = private$verbose)
				serial_inf$num_cores = dispatch_cores
				return(serial_inf$compute_rand_confidence_interval(
					alpha = alpha,
					r = r,
					pval_epsilon = pval_epsilon,
					show_progress = show_progress,
					type = type,
					args_for_type = args_for_type,
					ci_search_control = ci_search_control
				))
			}
			resp_type = private$des_obj_priv_int$response_type
			if (resp_type == "incidence" && is.null(private$custom_randomization_statistic_function)){
				rand_type = if (is.null(type)) "Zhang" else type
				exact_args = private$normalize_exact_inference_args(
					rand_type,
					args_for_type = args_for_type,
					pval_epsilon = pval_epsilon
				)
				return(private$compute_exact_confidence_interval_rand(rand_type, alpha, exact_args))
			}
			if (should_run_asserts()) {
				private$assert_no_incidence_only_randomization_args(resp_type, type, args_for_type)
			}
			is_glm = inherits(self, "InferenceAsympLikStdModCache") ||
			         inherits(self, "InferenceAsympLikStdModCacheNoParamBootstrap") ||
			         inherits(self, "InferenceCountLikelihoodNoParamBootstrap") ||
			         isTRUE(private$kk_gee_engine) ||
			         isTRUE(private$kk_glmm_engine) ||
			         isTRUE(private$kk_passthrough) ||
			         inherits(self, "InferencePropZeroOneInflatedBetaRegr") ||
			         inherits(self, "InferencePropGCompAbstract") ||
			         inherits(self, "InferenceCountZeroAugmentedPoissonAbstract") ||
			         inherits(self, "InferenceCountHurdleNegBin")
			temp_inf = if (resp_type %in% c("count", "proportion", "survival")) self$duplicate() else self
			transform_arg = "none"
			
			if (resp_type == "count"){
				transform_arg = if (is_glm) "log" else "already_transformed"
				if (!is_glm) temp_inf$.__enclos_env__$private$y = log1p(temp_inf$.__enclos_env__$private$y)
			} else if (resp_type == "proportion"){
				transform_arg = if (is_glm) "logit" else "already_transformed"
				y_clamped = pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, temp_inf$.__enclos_env__$private$y))
				temp_inf$.__enclos_env__$private$y = if (is_glm) y_clamped else logit(y_clamped)
			} else if (resp_type == "survival"){
				transform_arg = if (is_glm) "log" else "already_transformed"
				if (!is_glm) temp_inf$.__enclos_env__$private$y = log(pmax(.Machine$double.eps, temp_inf$.__enclos_env__$private$y))
			}
			if (resp_type %in% c("count", "proportion", "survival")) temp_inf$.__enclos_env__$private$cached_values = list()
			old_mc_control = temp_inf$.__enclos_env__$private$randomization_mc_control
			temp_inf$.__enclos_env__$private$randomization_mc_control = ci_search_control
			on.exit({ temp_inf$.__enclos_env__$private$randomization_mc_control = old_mc_control }, add = TRUE)
			perms = temp_inf$.__enclos_env__$private$generate_permutations(r)
			ci_pval_cache = if (isTRUE(ci_search_control$pval_cache_enable)) new.env(parent = emptyenv()) else NULL
			bounds = private$build_randomization_ci_search_bounds(temp_inf, r, alpha, transform_arg, perms, ci_search_control, ci_pval_cache)
			if (!all(is.finite(c(bounds$l, bounds$u)))) {
				fallback_ci = as.numeric(bounds$fallback_ci)
				fallback_mode = ci_search_control$fallback
				if (identical(fallback_mode, "error")) {
					stop("Randomization CI search failed to bracket the target p-value within the configured search radius.")
				}
				if (identical(fallback_mode, "na") || length(fallback_ci) < 2L || !all(is.finite(fallback_ci[1:2]))) {
					private$cache_nonestimable_se("rand_ci_search_bounds_failed")
					ci = c(NA_real_, NA_real_)
				} else {
					ci = sort(fallback_ci[1:2])
				}
				if (length(ci) != 2L) ci = c(NA_real_, NA_real_)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}
			# Run the lower and upper bounds sequentially and reserve all available
			# cores for the inner p-value computations inside each bisection step.
			ci = c(
				temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
					r, bounds$l, bounds$est, alpha / 2, pval_epsilon, transform_arg, TRUE, show_progress, perms, ci_search_control, ci_pval_cache
				),
				temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
					r, bounds$est, bounds$u, alpha / 2, pval_epsilon, transform_arg, FALSE, show_progress, perms, ci_search_control, ci_pval_cache
				)
			)
			if (length(ci) != 2L || !all(is.finite(ci[1:2]))) {
				private$cache_nonestimable_se("rand_ci_bisection_failed")
				ci = c(NA_real_, NA_real_)
			}
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		assert_no_incidence_only_randomization_args = function(resp_type, type, args_for_type){
			if (should_run_asserts()) {
				if (!is.null(type)) {
					stop("Randomization type dispatch is only supported for incidence outcomes.")
				}
				if (!is.null(args_for_type)) {
					stop("args_for_type is only used for incidence randomization inference.")
				}
			}
			invisible(NULL)
		},
		compute_randomization_ci_pval_cached = function(inf_obj, r, delta, transform_responses, permutations, ci_search_control, ci_pval_cache){
			cache_enabled = isTRUE(ci_search_control$pval_cache_enable) && !is.null(ci_pval_cache)
			if (cache_enabled) {
				cache_key = private$normalize_delta_for_cache(delta, ci_search_control$pval_cache_resolution)
				if (!is.null(ci_pval_cache[[cache_key]])) return(ci_pval_cache[[cache_key]])
			}
			pval = tryCatch(
				as.numeric(inf_obj$compute_rand_two_sided_pval(
					r,
					delta = delta,
					transform_responses = transform_responses,
					na.rm = FALSE,
					show_progress = FALSE,
					permutations = permutations
				)),
				error = function(e) NA_real_
			)
			if (length(pval) == 0L) pval = NA_real_ else pval = pval[1]
			if (cache_enabled) ci_pval_cache[[cache_key]] = pval
			pval
		},
		normalize_randomization_ci_search_control = function(ci_search_control, r, pval_epsilon){
			default_mc_batch = min(as.integer(r), max(25L, as.integer(ceiling(2 * sqrt(as.integer(r))))))
			defaults = list(
				fallback = "fallback",
				seed = "asymp_then_boot",
				max_radius_se_mult = 25,
				max_radius_scale_mult = 6,
				max_expansions = 7L,
				seed_boot_B = max(51L, min(as.integer(r), 201L)),
				pval_cache_enable = TRUE,
				pval_cache_resolution = pval_epsilon,
				mc_enable = as.integer(r) >= 200L,
				mc_batch_size = default_mc_batch,
				mc_min_draws = min(as.integer(r), max(100L, 2L * default_mc_batch)),
				mc_conf_level = 0.99,
				fit_warm_start_enable = TRUE,
				fit_reuse_factorizations = TRUE
			)
			if (should_run_asserts()) {
				assertList(ci_search_control, null.ok = TRUE)
			}
			ctrl = utils::modifyList(defaults, if (is.null(ci_search_control)) list() else ci_search_control)
			if (should_run_asserts()) {
				assertChoice(ctrl$fallback, c("fallback", "na", "error"))
				assertChoice(ctrl$seed, c("asymp_then_boot", "boot_then_asymp", "asymp_only", "boot_only", "none"))
			}
			if (should_run_asserts()) {
				assertNumber(ctrl$max_radius_se_mult, lower = 0, finite = TRUE)
				assertNumber(ctrl$max_radius_scale_mult, lower = 0, finite = TRUE)
			}
			if (should_run_asserts()) {
				assertCount(as.integer(ctrl$max_expansions), positive = TRUE)
				assertCount(as.integer(ctrl$seed_boot_B), positive = TRUE)
				assertFlag(ctrl$pval_cache_enable)
			}
			if (should_run_asserts()) {
				assertNumber(ctrl$pval_cache_resolution, lower = .Machine$double.eps, finite = TRUE)
			}
			if (should_run_asserts()) {
				assertFlag(ctrl$mc_enable)
				assertCount(as.integer(ctrl$mc_batch_size), positive = TRUE)
				assertCount(as.integer(ctrl$mc_min_draws), positive = TRUE)
			}
			if (should_run_asserts()) {
				assertNumber(ctrl$mc_conf_level, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (should_run_asserts()) {
				assertFlag(ctrl$fit_warm_start_enable)
				assertFlag(ctrl$fit_reuse_factorizations)
			}
			ctrl$max_expansions = as.integer(ctrl$max_expansions)
			ctrl$seed_boot_B = as.integer(ctrl$seed_boot_B)
			ctrl$pval_cache_resolution = as.numeric(ctrl$pval_cache_resolution)
			ctrl$mc_batch_size = min(as.integer(r), as.integer(ctrl$mc_batch_size))
			ctrl$mc_min_draws = min(as.integer(r), max(as.integer(ctrl$mc_batch_size), as.integer(ctrl$mc_min_draws)))
			ctrl
		},
		normalize_exact_inference_args = function(type, args_for_type = NULL, pval_epsilon = NULL){
			zhang_normalize_exact_inference_args(type, args_for_type = args_for_type, pval_epsilon = pval_epsilon)
		},
		assert_exact_inference_params = function(type, args_for_type){
			zhang_assert_exact_inference_params(self, type, args_for_type)
		},
		compute_exact_confidence_interval_rand = function(type, alpha, args_for_type){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Randomization inference")
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Zhang = zhang_ci_exact_combined(
					self,
					alpha = alpha,
					pval_epsilon = args_for_type[[type]]$pval_epsilon,
					combination_method = args_for_type[[type]]$combination_method
				)
			)
		},
		compute_exact_two_sided_pval_rand = function(type, delta, args_for_type){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Randomization inference")
				assertNumeric(delta)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Zhang = zhang_pval_exact_combined(
					self,
					delta_0 = delta,
					combination_method = args_for_type[[type]]$combination_method
				)
			)
		},
		build_randomization_ci_search_bounds = function(inf_obj, r, alpha, transform_arg, permutations, ci_search_control, ci_pval_cache){
			normalize_ci = function(ci){
				ci = as.numeric(ci)
				if (length(ci) < 2L || !all(is.finite(ci[1:2]))) return(c(NA_real_, NA_real_))
				sort(ci[1:2])
			}
			obj_private = inf_obj$.__enclos_env__$private
			if (transform_arg == "none" && (is.null(obj_private$cached_values$t0s_rand) || length(obj_private$cached_values$t0s_rand) < r)) {
				private$compute_randomization_ci_pval_cached(inf_obj, r, 0, transform_arg, permutations, ci_search_control, ci_pval_cache)
			}
			est = as.numeric(inf_obj$compute_estimate())
			if (length(est) == 0L || !is.finite(est[1])) est = NA_real_ else est = est[1]
			seed_cis = private$get_randomization_ci_seed_candidates(inf_obj, alpha)
			wald_ci = normalize_ci(seed_cis$wald_ci)
			asym_ci = normalize_ci(seed_cis$asym_ci)
			if (!is.finite(est) && all(is.finite(wald_ci))) est = mean(wald_ci)
			if (!is.finite(est) && all(is.finite(asym_ci))) est = mean(asym_ci)
			response_scale = stats::sd(obj_private$y, na.rm = TRUE)
			if (!is.finite(response_scale) || response_scale <= 0) response_scale = stats::IQR(obj_private$y, na.rm = TRUE) / 1.349
			if (!is.finite(response_scale) || response_scale <= 0) response_scale = 1
			if (!is.finite(est)) est = stats::median(obj_private$y, na.rm = TRUE)
			if (!is.finite(est)) est = 0
			se_guess = NA_real_
			if (all(is.finite(wald_ci))) se_guess = max(abs(wald_ci - est)) / max(stats::qnorm(1 - alpha), .Machine$double.eps)
			if (!is.finite(se_guess) && all(is.finite(asym_ci))) se_guess = max(abs(asym_ci - est)) / max(stats::qnorm(1 - alpha), .Machine$double.eps)
			if (!is.finite(se_guess) || se_guess <= 0) se_guess = response_scale / sqrt(max(1, obj_private$n))
			if (!is.finite(se_guess) || se_guess <= 0) se_guess = response_scale
			default_seed_ci = sort(c(est - 2 * se_guess, est + 2 * se_guess))
			fallback_ci = if (all(is.finite(wald_ci))) wald_ci else if (all(is.finite(asym_ci))) asym_ci else c(NA_real_, NA_real_)
			seed_ci = switch(ci_search_control$seed,
				asymp_then_boot = if (all(is.finite(wald_ci))) wald_ci else if (all(is.finite(asym_ci))) asym_ci else default_seed_ci,
				boot_then_asymp = if (all(is.finite(wald_ci))) wald_ci else if (all(is.finite(asym_ci))) asym_ci else default_seed_ci,
				asymp_only = if (all(is.finite(wald_ci))) wald_ci else if (all(is.finite(asym_ci))) asym_ci else c(NA_real_, NA_real_),
				boot_only = if (all(is.finite(wald_ci))) wald_ci else if (all(is.finite(asym_ci))) asym_ci else c(NA_real_, NA_real_),
				none = default_seed_ci
			)
			max_radius = max(ci_search_control$max_radius_se_mult * se_guess, ci_search_control$max_radius_scale_mult * response_scale, 1)
			min_radius = min(max(se_guess, 10 * .Machine$double.eps), max_radius)
			l = max(seed_ci[1], est - max_radius)
			u = min(seed_ci[2], est + max_radius)
			if (!is.finite(l) || l >= est) l = est - min_radius
			if (!is.finite(u) || u <= est) u = est + min_radius
			l = private$expand_bound(inf_obj, l, est, r, transform_arg, permutations, alpha / 2, TRUE, max_radius, ci_search_control$max_expansions, ci_search_control, ci_pval_cache)
			u = private$expand_bound(inf_obj, u, est, r, transform_arg, permutations, alpha / 2, FALSE, max_radius, ci_search_control$max_expansions, ci_search_control, ci_pval_cache)
			list(est = est, l = l, u = u, fallback_ci = fallback_ci)
		},
		get_randomization_ci_seed_candidates = function(inf_obj, alpha){
			normalize_ci = function(ci){
				ci = as.numeric(ci)
				if (length(ci) < 2L || !all(is.finite(ci[1:2]))) return(c(NA_real_, NA_real_))
				sort(ci[1:2])
			}
			wald_ci = c(NA_real_, NA_real_)
			asym_ci = c(NA_real_, NA_real_)
			if (is(inf_obj, "InferenceAsymp")) {
				asym_ci = normalize_ci(tryCatch(inf_obj$compute_asymp_confidence_interval(alpha = alpha * 2), error = function(e) c(NA_real_, NA_real_)))
				supported = tryCatch(inf_obj$get_supported_testing_types(), error = function(e) character())
				if ("wald" %in% supported) {
					old_testing_type = tryCatch(inf_obj$get_testing_type(), error = function(e) NULL)
					if (!is.null(old_testing_type)) {
						on.exit(tryCatch(inf_obj$set_testing_type(old_testing_type), error = function(e) NULL), add = TRUE)
					}
					wald_ci = normalize_ci(tryCatch({
						inf_obj$set_testing_type("wald")
						inf_obj$compute_asymp_confidence_interval(alpha = alpha * 2)
					}, error = function(e) c(NA_real_, NA_real_)))
				}
			}
			list(wald_ci = wald_ci, asym_ci = asym_ci)
		},
		expand_bound = function(inf_obj, bound, est, r, transform_arg, permutations, target_pval, lower, max_radius, max_expansions, ci_search_control, ci_pval_cache){
			evaluate_pval = function(delta) {
				private$compute_randomization_ci_pval_cached(inf_obj, r, delta, transform_arg, permutations, ci_search_control, ci_pval_cache)
			}
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
			# Conservative fallback: no bracket found within max_radius.
			# Return the search boundary so the bisection can detect and report it.
			if (lower) est - max_radius else est + max_radius
		},
			compute_ci_by_inverting_the_randomization_test_iteratively = function(r, l, u, pval_th, tol, transform_responses, lower, show_progress = TRUE, permutations = NULL, ci_search_control = NULL, ci_pval_cache = NULL){
			evaluate_pval = function(delta) {
				private$compute_randomization_ci_pval_cached(self, r, delta, transform_responses, permutations, ci_search_control, ci_pval_cache)
			}
			pval_l = evaluate_pval(l)
			pval_u = evaluate_pval(u)
			if (length(pval_l) == 0) pval_l = NA_real_; if (length(pval_u) == 0) pval_u = NA_real_
			for (k in seq_len(30L)) {
				if (!is.na(pval_l) && !is.na(pval_u)) break
				if (is.na(pval_l)) { l = (l + u) / 2; pval_l = evaluate_pval(l) }
				if (is.na(pval_u)) { u = (l + u) / 2; pval_u = evaluate_pval(u) }
			}
			if (is.na(pval_l) || is.na(pval_u) || !all(is.finite(c(l, u)))) return(NA_real_)
			# Conservative bound: p-value at the search boundary still >= alpha/2,
			# meaning the true CI bound lies beyond the search radius.
			# Return the boundary as a valid (conservative) one-sided limit.
			if (lower && is.finite(pval_l) && pval_l >= pval_th) {
				message(sprintf(
					"Randomization CI lower bound is conservative: p-value at search boundary delta=%.4g is %.4g >= %.4g. True CI lower bound may extend further left. (rand_ci_conservative_count++)",
					l, pval_l, pval_th))
				private[["rand_ci_conservative_count"]] = (if (is.null(private[["rand_ci_conservative_count"]])) 0L else private[["rand_ci_conservative_count"]]) + 1L
				return(l)
			}
			if (!lower && is.finite(pval_u) && pval_u >= pval_th) {
				message(sprintf(
					"Randomization CI upper bound is conservative: p-value at search boundary delta=%.4g is %.4g >= %.4g. True CI upper bound may extend further right. (rand_ci_conservative_count++)",
					u, pval_u, pval_th))
				private[["rand_ci_conservative_count"]] = (if (is.null(private[["rand_ci_conservative_count"]])) 0L else private[["rand_ci_conservative_count"]]) + 1L
				return(u)
			}
			iter = 0; progress_label = if (lower) "CI lower" else "CI upper"
			repeat {
				pval_span = abs(pval_u - pval_l)
				if ((abs(u - l)) <= tol || pval_span <= tol) {
					if (isTRUE(show_progress)) cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g) done\n", progress_label, iter, pval_span, tol))
					return(if(lower) l else u)
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
