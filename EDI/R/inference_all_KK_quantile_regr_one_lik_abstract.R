#' Abstract Quantile Regression Combined-Likelihood Compound Estimator for KK Designs
#'
#' Fits a single joint quantile regression over all KK design data by stacking
#' matched-pair differences and reservoir observations into one design matrix.
#'
#' Column layout of X_stack: [beta_0 | beta_T | beta_xs (p cols)]
#' Pair rows:      [0 | 1   | Xd_k] -> Q_tau(yd_k) = beta_T + Xd_k' beta_xs
#' Reservoir rows: [1 | w_i | X_i]  -> Q_tau(y_i)  = beta_0 + w_i*beta_T + X_i'*beta_xs
#' Fitting a single rq() on the stacked dataset minimises the combined check-function loss.
#'
#' Special cases:
#' Pairs only:     beta_0 column is all-zero and dropped; layout [beta_T | beta_xs].
#' Reservoir only: standard quantile regression layout [beta_0 | beta_T | beta_xs].
#'
#' Standard errors use Powell's "nid" sandwich estimator, falling back to "iid".
#'
#' @keywords internal
InferenceAbstractKKQuantileRegrOneLik = R6::R6Class("InferenceAbstractKKQuantileRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractQuantileRandCI,
	public = list(
		#' @description Initialize KK quantile-regression combined-likelihood inference.
		#' @param des_obj A completed KK design object.
		#' @param tau The quantile level for regression, strictly between 0 and 1. The default
		#'   \code{tau = 0.5} estimates the median treatment effect. Values of exactly 0 or 1
		#'   are excluded because quantile regression is undefined at the boundary (the check
		#'   function \eqn{\rho_\tau(u) = u(\tau - \mathbf{1}_{u < 0})} degenerates there);
		#'   the bound is enforced as \code{(.Machine$double.eps, 1 - .Machine$double.eps)},
		#'   i.e. the smallest representable positive number away from 0 and 1.
		#' @param transform_y_fn Optional response transformation.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @return A new inference object.
		initialize = function(des_obj, model_formula = NULL, tau = 0.5, transform_y_fn = identity,  verbose = FALSE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertNumeric(tau, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			}
			if (should_run_asserts()) {
				if (!check_package_installed("quantreg")) {
					stop("Package 'quantreg' is required. Please install it with install.packages(\"quantreg\").")
				}
			}
			private$tau = tau
			private$transform_y_fn_list = list(fn = transform_y_fn)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (private$is_KK){
				private$m = des_obj$.__enclos_env__$private$m
				private$compute_basic_match_data()
			}
		},
		#' @description Compute the quantile-regression treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute the treatment estimate with bootstrap weights.
		#' @param subject_or_block_weights Numeric vector. Row weights for bootstrap.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		#' @return The treatment estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				self$compute_estimate(estimate_only = estimate_only)
				beta_hat_T = as.numeric(private$cached_values$beta_hat_T)[1L]
				if (is.finite(beta_hat_T)) {
					return(private$cached_values$beta_hat_T)
				}
			}
			result = private$compute_weighted_combined_estimate(row_weights, estimate_only = estimate_only)
			if (is.list(result)) {
				private$cached_values$beta_hat_T = result$beta
				private$cached_values$s_beta_hat_T = result$se
			} else {
				private$cached_values$beta_hat_T = result
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$beta_hat_T
		},
		#' @description Compute an asymptotic confidence interval.
		#' @param alpha Significance level.
		#' @return A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute an asymptotic two-sided p-value for the treatment effect.
		#' @param delta Null treatment effect value.
		#' @return A two-sided p-value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		tau = NULL,
		transform_y_fn_list = NULL,  # list(fn = ...) wrapping avoids R6 treating function as a locked method
		m = NULL,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		compute_basic_match_data = function(){
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv     = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},
		assert_finite_se = function(){
		},
		get_standard_error = function(){
			private$shared_combined_likelihood(estimate_only = FALSE)
			se = private$cached_values$s_beta_hat_T
			if (is.null(se) || length(se) == 0L) {
				return(NA_real_)
			}
			as.numeric(se)[1L]
		},
		# Fit the combined check-function loss over matched-pair differences and
		# reservoir observations with SHARED covariate effects beta_xs.
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			tau = private$tau
			fn  = private$transform_y_fn_list$fn
			has_reservoir = nRT > 0 && nRC > 0
			y_stack  = NULL
			X_stack  = NULL
			j_beta_T = 2L
			if (m > 0){
				yd = fn(KKstats$yTs_matched) - fn(KKstats$yCs_matched)
				if (has_reservoir){
					# Combined case: pair rows and reservoir rows must share the same
					# covariate columns. Use the full-width pair differences from the
					# shared C++ preprocessing rather than the reduced pair-only matrix.
					Xd = as.matrix(KKstats$X_matched_diffs_full)
					p   = ncol(Xd)
					y_r = fn(KKstats$y_reservoir)
					w_r = KKstats$w_reservoir
					X_r = as.matrix(KKstats$X_reservoir)
					X_pairs = if (p > 0) cbind(0, 1, Xd) else cbind(numeric(m), rep(1.0, m))
					X_res   = if (p > 0) cbind(1, w_r, X_r) else cbind(1, w_r)
					X_stack  = rbind(X_pairs, X_res)
					y_stack  = c(yd, y_r)
					j_beta_T = 2L
				} else {
					# Pairs only: drop all-zero beta_0 column; intercept = beta_T.
					Xd = as.matrix(KKstats$X_matched_diffs)
					p  = ncol(Xd)
					X_stack  = if (p > 0) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1)
					y_stack  = yd
					j_beta_T = 1L
				}
			} else if (has_reservoir){
				y_r = fn(KKstats$y_reservoir)
				w_r = KKstats$w_reservoir
				X_r = as.matrix(KKstats$X_reservoir)
				p   = ncol(X_r)
				X_stack  = if (p > 0) cbind(1, w_r, X_r) else cbind(1, w_r)
				y_stack  = y_r
				j_beta_T = 2L
			}
			if (is.null(X_stack)){
				private$cache_nonestimable_estimate("no_usable_matched_or_reservoir_data")
				return(invisible(NULL))
			}
			# QR-reduce to full rank, preserving beta_T column
			reduced = qr_reduce_preserve_cols_cpp(X_stack, j_beta_T)
			X_stack  = reduced$X_reduced
			j_beta_T = match(j_beta_T, reduced$keep)
			n_total  = nrow(X_stack)
			n_params = ncol(X_stack)
			if (n_total <= n_params){
				private$cache_nonestimable_estimate("insufficient_data_for_quantile_regr")
				return(invisible(NULL))
			}
			cn = paste0("x", seq_len(n_params))
			cn[j_beta_T] = "trt__"
			colnames(X_stack) = cn
			dat = as.data.frame(X_stack)
			dat$y_stack__ = y_stack
			fit = tryCatch(
				suppressWarnings(quantreg::rq(y_stack__ ~ . - 1, tau = tau, data = dat)),
				error = function(e) NULL
			)
			if (is.null(fit)){
				private$cache_nonestimable_estimate("quantile_regr_fit_unavailable")
				return(invisible(NULL))
			}
			beta = tryCatch(coef(fit)[["trt__"]], error = function(e) NA_real_)
			if (!is.finite(beta)) {
				private$cache_nonestimable_estimate("quantile_regr_nonfinite_coef")
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T   = beta
			if (!estimate_only) {
				se = private$extract_se_from_rq(fit, "trt__")
				private$cached_values$s_beta_hat_T = if (!is.na(se)) se else NA_real_
			}
			invisible(NULL)
		},
		# Helper: extract SE from rq fit by coefficient name, trying "nid" then "iid".
		# SEs above 1e6 are treated as invalid (the "nid" sparsity estimator can return
		# astronomically large but finite values when the density at the quantile is near
		# zero, which bypasses the usual !is.finite() || <= 0 guard in callers).
		extract_se_from_rq = function(fit, coef_name){
			.extract_se_from_rq_fit(fit, coef_name)
		},
		compute_weighted_combined_estimate = function(row_weights, estimate_only = TRUE){
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(NA_real_)
			kk_w = kk_pair_and_reservoir_bootstrap_weights(private, row_weights)
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			tau = private$tau
			fn  = private$transform_y_fn_list$fn
			has_reservoir = nRT > 0 && nRC > 0
			y_stack  = NULL
			X_stack  = NULL
			w_stack  = NULL
			j_beta_T = 2L
			if (m > 0){
				yd = fn(KKstats$yTs_matched) - fn(KKstats$yCs_matched)
				if (has_reservoir){
					Xd = as.matrix(KKstats$X_matched_diffs_full)
					p = ncol(Xd)
					y_r = fn(KKstats$y_reservoir)
					w_r = KKstats$w_reservoir
					X_r = as.matrix(KKstats$X_reservoir)
					X_pairs = if (p > 0) cbind(0, 1, Xd) else cbind(numeric(m), rep(1.0, m))
					X_res = if (p > 0) cbind(1, w_r, X_r) else cbind(1, w_r)
					X_stack = rbind(X_pairs, X_res)
					y_stack = c(yd, y_r)
					w_stack = c(kk_w$pair_weights, kk_w$reservoir_weights)
					j_beta_T = 2L
				} else {
					Xd = as.matrix(KKstats$X_matched_diffs)
					p = ncol(Xd)
					X_stack = if (p > 0) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1)
					y_stack = yd
					w_stack = kk_w$pair_weights
					j_beta_T = 1L
				}
			} else if (has_reservoir){
				y_r = fn(KKstats$y_reservoir)
				w_r = KKstats$w_reservoir
				X_r = as.matrix(KKstats$X_reservoir)
				p = ncol(X_r)
				X_stack = if (p > 0) cbind(1, w_r, X_r) else cbind(1, w_r)
				y_stack = y_r
				w_stack = kk_w$reservoir_weights
				j_beta_T = 2L
			}
			if (is.null(X_stack) || is.null(w_stack)) return(NA_real_)
			ok = is.finite(w_stack) & w_stack > 0 & is.finite(y_stack)
			if (!any(ok)) return(NA_real_)
			X_stack = X_stack[ok, , drop = FALSE]
			y_stack = y_stack[ok]
			w_stack = w_stack[ok]
			reduced = qr_reduce_preserve_cols_cpp(X_stack, j_beta_T)
			X_stack = reduced$X_reduced
			j_beta_T = match(j_beta_T, reduced$keep)
			if (nrow(X_stack) <= ncol(X_stack)) return(NA_real_)
			cn = paste0("x", seq_len(ncol(X_stack)))
			cn[j_beta_T] = "trt__"
			colnames(X_stack) = cn
			dat = as.data.frame(X_stack)
			dat$y_stack__ = y_stack
			fit = tryCatch(
				suppressWarnings(quantreg::rq(y_stack__ ~ . - 1, tau = tau, data = dat, weights = w_stack)),
				error = function(e) NULL
			)
			if (is.null(fit)) return(if (estimate_only) NA_real_ else list(beta = NA_real_, se = NA_real_))
			beta = tryCatch(coef(fit)[["trt__"]], error = function(e) NA_real_)
			beta = if (is.finite(beta)) as.numeric(beta) else NA_real_
			if (estimate_only) return(beta)
			se = tryCatch(private$extract_se_from_rq(fit, "trt__"), error = function(e) NA_real_)
			list(beta = beta, se = if (!is.na(se) && se > 0) se else NA_real_)
		}
	)
)
