#' Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a single stacked robust regression over matched-pair differences and reservoir
#' observations for KK matching-on-the-fly designs with continuous responses.
#'
#' @export
InferenceContinKKRobustRegrOneLik = R6::R6Class("InferenceContinKKRobustRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompoundNoParamBootstrap,
	public = list(
				
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param method  		Robust-regression fitting method for `MASS::rlm`; one of `"M"` or `"MM"`.
		#' @param maxit  		Maximum number of robust-regression iterations. If `NULL`, a
		#'   data-adaptive default is chosen at fit time.
		#' @param acc  			Convergence tolerance for `MASS::rlm`. If `NULL`, a
		#'   data-adaptive default is chosen at fit time.
		#' @param start_with_ols  Whether to compute an OLS warm start and pass it to `MASS::rlm`
		#'   when the fit method honors `init`. This affects the `M` path only; `MM`
		#'   uses its own LQS-based start. Default `TRUE`.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, method = "MM", maxit = NULL, acc = NULL, start_with_ols = TRUE, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertChoice(method, c("M", "MM"))
				if (!is.null(maxit)) assertCount(maxit, positive = TRUE)
				if (!is.null(acc)) assertNumeric(acc, lower = .Machine$double.xmin, upper = 1)
				assertFlag(start_with_ols)
				assertFormula(model_formula, null.ok = TRUE)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			
			private$rlm_method = method
			private$rlm_maxit = maxit
			private$rlm_acc = acc
			private$rlm_start_with_ols = start_with_ols
		},
		#' @description Returns the combined robust-regression estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations and use
		#'   the faster \code{"M"} robust estimator regardless of the \code{method}
		#'   argument passed at construction time.
		compute_estimate = function(estimate_only = FALSE){
			private$fit_combined(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the combined robust-regression estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				self$compute_estimate(estimate_only = estimate_only)
				beta_hat_T = as.numeric(private$cached_values$beta_hat_T)[1L]
				if (is.finite(beta_hat_T)) {
					return(private$cached_values$beta_hat_T)
				}
			}
			result = private$fit_weighted_combined(row_weights, estimate_only = estimate_only)
			if (is.list(result)) {
				private$cached_values$beta_hat_T = result$beta
				private$cached_values$s_beta_hat_T = result$se
			} else {
				private$cached_values$beta_hat_T = result
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$beta_hat_T
		},
		#' @description Computes the approximate confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 -
		#'   \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$fit_combined()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the approximate p-value.
		#' @param delta The null difference to test against. For any treatment effect at all this
		#'   is set to zero (the default).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$fit_combined()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes the Wald confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 -
		#'   \code{alpha}. The default is 0.05.
		compute_wald_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$fit_combined()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the Wald two-sided p-value.
		#' @param delta The null difference to test against. Default 0.
		compute_wald_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$fit_combined()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Duplicate
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
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
			preserve = if (is.null(permutations$m_mat)) c("kk_robust_combined_reduced_design") else character()
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp, preserve_cache_keys = preserve)
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
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		get_standard_error = function(){
			private$fit_combined(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T %||% NA_real_
		},
		get_degrees_of_freedom = function(){
			private$fit_combined(estimate_only = FALSE)
			private$cached_values$df %||% NA_real_
		},
		# estimate_only = TRUE forces "M" (fast, no LQS phase) and skips summary().
		fit_rlm = function(X, y, j_treat, estimate_only = FALSE){
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
			# MASS::rlm only honors `init` for the M path; MM uses its own LQS start.
			# If no OLS warm start is available, explicitly start the M fit at zero.
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
			# Fast path: extract coefficient directly without calling summary().
			if (estimate_only) {
				coef_vec = tryCatch(as.numeric(stats::coef(mod)), error = function(e) NULL)
				if (is.null(coef_vec) || length(coef_vec) < j_treat) return(NULL)
				beta = coef_vec[j_treat]
				if (!is.finite(beta)) return(NULL)
				return(list(beta = beta, se = NA_real_, mod = NULL))
			}
			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if ((is.null(coef_table) || nrow(coef_table) < j_treat) && identical(method_to_try, "MM")){
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
			if ((!is.finite(beta) || !is.finite(se) || se <= 0) && identical(method_to_try, "MM")){
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
			list(beta = beta, se = se, mod = mod)
		},
		fit_combined = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			nR  = nRT + nRC
			if (m > 0 && nRT > 0 && nRC > 0){
				if (ncol(as.matrix(private$X)) > 0){
					Xd_full = as.matrix(KKstats$X_matched_diffs_full)
					p = ncol(Xd_full)
				} else {
					Xd_full = matrix(nrow = m, ncol = 0L)
					p = 0L
				}
				X_comb = rbind(
					if (p > 0L) cbind(0, 1, Xd_full) else cbind(numeric(m), rep(1.0, m)),
					if (p > 0L) cbind(rep(1, nR), KKstats$w_reservoir, as.matrix(KKstats$X_reservoir)) else cbind(rep(1, nR), KKstats$w_reservoir)
				)
				y_comb = c(KKstats$y_matched_diffs, KKstats$y_reservoir)
				j_treat = 2L
			} else if (m > 0){
				Xd = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_matched_diffs) else matrix(nrow = m, ncol = 0L)
				X_comb = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				y_comb = KKstats$y_matched_diffs
				j_treat = 1L
			} else if (nRT > 0 && nRC > 0){
				X_r = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_reservoir) else matrix(nrow = nR, ncol = 0L)
				X_comb = if (ncol(X_r) > 0L) cbind(rep(1, nR), KKstats$w_reservoir, X_r) else cbind(rep(1, nR), KKstats$w_reservoir)
				y_comb = KKstats$y_reservoir
				j_treat = 2L
			} else {
				private$cache_nonestimable_estimate("no_usable_matched_or_reservoir_data")
				return(invisible(NULL))
			}
			reduced = private$reduce_design_matrix_once(
				X_comb,
				j_treat,
				cache_key = "kk_robust_combined_reduced_design"
			)
			X_comb = reduced$X
			j_treat = reduced$j_treat
			fit = private$fit_rlm(X_comb, y_comb, j_treat, estimate_only = estimate_only)
			if (is.null(fit)){
				private$cache_nonestimable_estimate("rlm_fit_unavailable")
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T   = fit$beta
			if (!estimate_only) {
				private$cached_values$s_beta_hat_T = fit$se
				if (!is.null(fit$mod)) {
					private$cached_values$full_coefficients = stats::coef(fit$mod)
					private$cached_values$full_vcov = tryCatch(stats::vcov(fit$mod), error = function(e) NULL)
				}
			}
			invisible(NULL)
		},
		fit_weighted_combined = function(row_weights, estimate_only = TRUE){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}
			if (is.null(KKstats)) return(NA_real_)
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			nR  = nRT + nRC
			kk_w = kk_pair_and_reservoir_bootstrap_weights(private, row_weights)
			if (m > 0 && nRT > 0 && nRC > 0){
				if (ncol(as.matrix(private$X)) > 0){
					Xd_full = as.matrix(KKstats$X_matched_diffs_full)
					p = ncol(Xd_full)
				} else {
					Xd_full = matrix(nrow = m, ncol = 0L)
					p = 0L
				}
				X_comb = rbind(
					if (p > 0L) cbind(0, 1, Xd_full) else cbind(numeric(m), rep(1.0, m)),
					if (p > 0L) cbind(rep(1, nR), KKstats$w_reservoir, as.matrix(KKstats$X_reservoir)) else cbind(rep(1, nR), KKstats$w_reservoir)
				)
				y_comb = c(KKstats$y_matched_diffs, KKstats$y_reservoir)
				w_comb = c(kk_w$pair_weights, kk_w$reservoir_weights)
				j_treat = 2L
			} else if (m > 0){
				Xd = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_matched_diffs) else matrix(nrow = m, ncol = 0L)
				X_comb = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				y_comb = KKstats$y_matched_diffs
				w_comb = kk_w$pair_weights
				j_treat = 1L
			} else if (nRT > 0 && nRC > 0){
				X_r = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_reservoir) else matrix(nrow = nR, ncol = 0L)
				X_comb = if (ncol(X_r) > 0L) cbind(rep(1, nR), KKstats$w_reservoir, X_r) else cbind(rep(1, nR), KKstats$w_reservoir)
				y_comb = KKstats$y_reservoir
				w_comb = kk_w$reservoir_weights
				j_treat = 2L
			} else {
				return(NA_real_)
			}
			ok = is.finite(w_comb) & w_comb > 0 & is.finite(y_comb)
			if (!any(ok)) return(NA_real_)
			X_comb = X_comb[ok, , drop = FALSE]
			y_comb = y_comb[ok]
			w_comb = w_comb[ok]
			reduced = private$reduce_design_matrix_once(X_comb, j_treat, cache_key = "kk_robust_combined_reduced_design_weighted")
			X_comb = reduced$X
			j_treat = reduced$j_treat
			fit = tryCatch(
				suppressWarnings(MASS::rlm(x = X_comb, y = y_comb, weights = w_comb, method = "M", maxit = 20L, acc = 1e-4)),
				error = function(e) NULL
			)
			coef_vec = tryCatch(as.numeric(stats::coef(fit)), error = function(e) NULL)
			if (!is.null(coef_vec) && length(coef_vec) >= j_treat && is.finite(coef_vec[j_treat])) {
				beta = coef_vec[j_treat]
				if (estimate_only) return(beta)
				se = tryCatch({
					ct = summary(fit)$coefficients
					if (!is.null(ct) && nrow(ct) >= j_treat) {
						se_val = as.numeric(ct[j_treat, "Std. Error"])
						if (is.finite(se_val) && se_val > 0) se_val else NA_real_
					} else NA_real_
				}, error = function(e) NA_real_)
				return(list(beta = beta, se = se))
			}
			lm_fit = tryCatch(stats::lm.wfit(x = X_comb, y = y_comb, w = w_comb), error = function(e) NULL)
			coef_lm = if (is.null(lm_fit)) NULL else as.numeric(lm_fit$coefficients)
			if (is.null(coef_lm) || length(coef_lm) < j_treat || !is.finite(coef_lm[j_treat])) {
				return(if (estimate_only) NA_real_ else list(beta = NA_real_, se = NA_real_))
			}
			beta = coef_lm[j_treat]
			if (estimate_only) return(beta)
			df = nrow(X_comb) - ncol(X_comb)
			se = NA_real_
			if (df > 0) {
				sw = sqrt(w_comb)
				post_fit = tryCatch(
					ols_hc2_post_fit_cpp(X_comb * sw, as.numeric(y_comb * sw), coef_lm, j_treat),
					error = function(e) NULL
				)
				if (!is.null(post_fit)) {
					se_val = as.numeric(post_fit$std_err[j_treat])
					if (is.finite(se_val) && se_val > 0) se = se_val
				}
			}
			list(beta = beta, se = se)
		}
	)
)
#' Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a single stacked robust regression over matched-pair differences and
#' reservoir observations for KK matching-on-the-fly designs with continuous
#' responses, using the treatment indicator and, optionally, all recorded
#' covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinKKRobustRegrOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
