#' Robust-Regression IVWC Compound Inference for KK Designs
#'
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using robust linear regression (`MASS::rlm`) for the
#' matched-pair and reservoir components separately.
#'
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferenceContinKKRobustRegrIVWC = R6::R6Class("InferenceContinKKRobustRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompoundNoParamBootstrap,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param method  		Robust-regression fitting method for `MASS::rlm`; one of `"M"` or `"MM"`.
		#' @param maxit  		Maximum number of robust-regression iterations. If `NULL`, a
		#'   data-adaptive default is chosen at fit time.
		#' @param acc  			Convergence tolerance for `MASS::rlm`. If `NULL`, a
		#'   data-adaptive default is chosen at fit time.
		#' @param start_with_ols  Whether to compute an OLS warm start and pass it to `MASS::rlm`
		#'   when the fit method honors `init`. This affects the `M` path only; `MM`
		#'   uses its own LQS-based start. Default `TRUE`.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart starting values for the optimizer.
		initialize = function(des_obj, model_formula = NULL, method = "MM", maxit = NULL, acc = NULL, start_with_ols = TRUE, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertFormula(model_formula, null.ok = TRUE)
				assertChoice(method, c("M", "MM"))
				if (!is.null(maxit)) assertCount(maxit, positive = TRUE)
				if (!is.null(acc)) assertNumeric(acc, lower = .Machine$double.xmin, upper = 1)
				assertFlag(start_with_ols)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			private$rlm_method = method
			private$rlm_maxit = maxit
			private$rlm_acc = acc
			private$rlm_start_with_ols = start_with_ols
		},
		#' @description Returns the estimated treatment effect.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the approximate confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the approximate p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Duplicate
		#' @param verbose Logical. A flag indicating whether messages should be displayed.
		#' @param make_fork_cluster Logical. Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			i = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			i
		}
	),
	private = list(
		rlm_method = NULL,
		rlm_maxit = NULL,
		rlm_acc = NULL,
		rlm_start_with_ols = TRUE,
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			preserve = if (is.null(permutations$m_mat)) c("kk_robust_ivwc_matched_reduced_design", "kk_robust_ivwc_reservoir_reduced_design") else character()
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, preserve_cache_keys = preserve, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		rlm_force_M = FALSE,
				
		reduce_design_matrix_once = function(X, j_treat, cache_key){
			cached = private$cached_values[[cache_key]]
			if (!is.null(cached)) return(cached)
			qr_X = qr(X)
			if (qr_X$rank < ncol(X)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_treat %in% keep)) keep[qr_X$rank] = j_treat
				keep = sort(unique(keep))
				X = X[, keep, drop = FALSE]
				j_treat = which(keep == j_treat)
			}
			cached = list(X = X, j_treat = j_treat)
			private$cached_values[[cache_key]] = cached
			cached
		},
		resolve_rlm_control = function(X){
			p = ncol(X)
			maxit = private$rlm_maxit
			if (is.null(maxit)) {
				maxit = if (isTRUE(private$rlm_start_with_ols)) {
					if (p <= 3L) 10L else if (p <= 10L) 15L else 20L
				} else {
					if (p <= 3L) 15L else if (p <= 10L) 20L else 25L
				}
			}
			acc = private$rlm_acc
			if (is.null(acc)) {
				acc = if (isTRUE(private$rlm_start_with_ols)) 1e-4 else 1e-3
			}
			list(maxit = as.integer(maxit), acc = as.numeric(acc))
		},
		is_rlm_nonconvergence_warning = function(w){
			msg = conditionMessage(w)
			grepl("'rlm' failed to converge", msg, fixed = TRUE) ||
				grepl("alternation limit reached", msg, fixed = TRUE)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (m > 0){
				private$robust_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m  = private$cached_values$ssq_beta_T_matched
			m_ok   = !is.null(beta_m) && is.finite(beta_m) && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0
			if (nRT > 0 && nRC > 0){
				private$robust_for_reservoir(estimate_only = estimate_only)
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r  = private$cached_values$ssq_beta_T_reservoir
			r_ok   = !is.null(beta_r) && is.finite(beta_r) && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
			}
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		# estimate_only = TRUE forces "M" (fast, no LQS phase) for the point estimate.
		fit_rlm_with_treatment = function(X, y, j_treat, estimate_only = FALSE){
			if (nrow(X) <= ncol(X)) return(NULL)
			ctrl = private$resolve_rlm_control(X)
			run_rlm = function(method, init = NULL){
				nonconverged = FALSE
				tryCatch({
					mod = withCallingHandlers({
						if (identical(method, "M")) {
							args = list(
								x = X,
								y = y,
								method = "M",
								psi = MASS::psi.huber,
								maxit = ctrl$maxit,
								acc = ctrl$acc
							)
							if (!is.null(init)) args$init = init
							do.call(MASS::rlm, args)
						} else {
							do.call(MASS::rlm, list(
								x = X,
								y = y,
								method = method,
								maxit = ctrl$maxit,
								acc = ctrl$acc
							))
						}
					}, warning = function(w){
						if (private$is_rlm_nonconvergence_warning(w)) {
							nonconverged <<- TRUE
							invokeRestart("muffleWarning")
						}
					})
					list(mod = mod, nonconverged = nonconverged)
				}, error = function(e) e)
			}
			# When only the point estimate is needed, skip the expensive MM/LQS phase.
			method_to_try = if (estimate_only || isTRUE(private$rlm_force_M)) "M" else private$rlm_method
			start_coef = NULL
			if (isTRUE(private$rlm_start_with_ols)) {
				start_coef = tryCatch(as.numeric(stats::coef(stats::lm.fit(x = X, y = y))), error = function(e) NULL)
				if (!is.null(start_coef) && (length(start_coef) != ncol(X) || any(!is.finite(start_coef)))) {
					start_coef = NULL
				}
				if (!is.null(start_coef)) start_coef = list(coef = unname(start_coef))
			}
			if (identical(method_to_try, "M") && is.null(start_coef)) {
				start_coef = rep(0, ncol(X))
			}
			rlm_attempt = run_rlm(method_to_try, init = start_coef)
			if (!estimate_only && identical(method_to_try, "MM") &&
					(inherits(rlm_attempt, "error") || isTRUE(rlm_attempt$nonconverged))){
				msg = if (inherits(rlm_attempt, "error") && length(rlm_attempt$message) > 0L) rlm_attempt$message else ""
				if (isTRUE(rlm_attempt$nonconverged) || grepl("'lqs' failed", msg, fixed = TRUE) || grepl("singular", msg, ignore.case = TRUE)) {
					private$rlm_force_M = TRUE
					rlm_attempt = run_rlm("M", init = start_coef)
				}
			}
			if (inherits(rlm_attempt, "error") || isTRUE(rlm_attempt$nonconverged)) return(NULL)
			mod = rlm_attempt$mod
			if (is.null(mod)) return(NULL)
			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if (!estimate_only && (is.null(coef_table) || nrow(coef_table) < j_treat) && identical(method_to_try, "MM")){
				private$rlm_force_M = TRUE
				rlm_attempt = run_rlm("M", init = start_coef)
				if (inherits(rlm_attempt, "error") || isTRUE(rlm_attempt$nonconverged)) return(NULL)
				mod = rlm_attempt$mod
				if (is.null(mod)) return(NULL)
				coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			}
			if (is.null(coef_table) || nrow(coef_table) < j_treat) return(NULL)
			beta = as.numeric(coef_table[j_treat, "Value"])
			se   = as.numeric(coef_table[j_treat, "Std. Error"])
			if (!estimate_only && (!is.finite(beta) || !is.finite(se) || se <= 0) && identical(method_to_try, "MM")){
				private$rlm_force_M = TRUE
				rlm_attempt = run_rlm("M", init = start_coef)
				if (inherits(rlm_attempt, "error") || isTRUE(rlm_attempt$nonconverged)) return(NULL)
				mod = rlm_attempt$mod
				if (is.null(mod)) return(NULL)
				coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
				if (is.null(coef_table) || nrow(coef_table) < j_treat) return(NULL)
				beta = as.numeric(coef_table[j_treat, "Value"])
				se   = as.numeric(coef_table[j_treat, "Std. Error"])
			}
			if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)
			list(beta = beta, ssq = se^2)
		},
		robust_for_matched_pairs = function(estimate_only = FALSE){
			yd = private$cached_values$KKstats$y_matched_diffs
			m  = length(yd)
			if (ncol(as.matrix(private$X)) > 0){
				Xd = as.matrix(private$cached_values$KKstats$X_matched_diffs)
				X = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				reduced = private$reduce_design_matrix_once(
					X,
					1L,
					cache_key = "kk_robust_ivwc_matched_reduced_design"
				)
				X = reduced$X
			} else {
				X = matrix(1, nrow = m, ncol = 1L)
			}
			fit = private$fit_rlm_with_treatment(X, yd, 1L, estimate_only = estimate_only)
			if (is.null(fit)) {
				private$cached_values$beta_T_matched     = if (m >= 1) mean(yd) else NA_real_
				private$cached_values$ssq_beta_T_matched = if (m >= 2) var(yd) / m else NA_real_
			} else {
				private$cached_values$beta_T_matched     = fit$beta
				private$cached_values$ssq_beta_T_matched = fit$ssq
			}
		},
		robust_for_reservoir = function(estimate_only = FALSE){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			X_r = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L
			if (ncol(as.matrix(private$X)) > 0){
				X_full = cbind(1, w_r, X_r)
				reduced = private$reduce_design_matrix_once(
					X_full,
					j_treat,
					cache_key = "kk_robust_ivwc_reservoir_reduced_design"
				)
				X_full = reduced$X
				j_treat = reduced$j_treat
			} else {
				X_full = cbind(1, w_r)
			}
			fit = private$fit_rlm_with_treatment(X_full, y_r, j_treat, estimate_only = estimate_only)
			if (is.null(fit)) {
				private$cached_values$beta_T_reservoir     = NA_real_
				private$cached_values$ssq_beta_T_reservoir = NA_real_
			} else {
				private$cached_values$beta_T_reservoir     = fit$beta
				private$cached_values$ssq_beta_T_reservoir = fit$ssq
			}
		}
	)
)
#' Robust-Regression IVWC Compound Inference for KK Designs
#'
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using robust regression for matched-pair differences
#' and reservoir outcomes, with treatment and, optionally, all recorded covariates
#' as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinKKRobustRegrIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
