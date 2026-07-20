#' Stratified Cox / Standard Cox Compound Inference for KK Designs
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses. For matched pairs, it uses stratified Cox proportional hazards
#' regression (each pair is a stratum). For reservoir subjects, it uses standard Cox
#' regression. The two estimates (both log-hazard ratios) are combined via a
#' variance-weighted linear combination.
#'
#' Under \code{harden = TRUE}, multivariate fits preserve the treatment column and
#' progressively retry reduced covariate sets after QR-based rank reduction and
#' correlation-based pruning. Extreme finite coefficients / standard errors are
#' rejected and treated as non-estimable.
#'
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferenceSurvivalKKStratCoxPHIVWC = R6::R6Class("InferenceSurvivalKKStratCoxPHIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the estimated treatment effect (log-hazard ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the asymptotic confidence interval.
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
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
		#' @description Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				if (should_run_asserts()) {
					stop("Testing non-zero delta is not yet implemented for this class.")
				}
				NA_real_
			}
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		max_abs_reasonable_coef = 1e4,
		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		cox_design_candidates = function(w, X){
			X_full = matrix(w, ncol = 1)
			colnames(X_full) = "w"
			X_covs = as.matrix(X)
			if (ncol(X_covs) > 0L){
				X_full = cbind(X_full, X_covs)
			}
			if (!private$harden || ncol(X_full) <= 1L){
				return(list(X_full))
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit) X_fit,
				fit_ok = function(mod, X_fit, keep) TRUE
			)
			candidates = list(attempt$X)
			keys = paste(colnames(candidates[[1L]]), collapse = "|")
			thresholds = c(0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			X_cov_orig = X_full[, -1, drop = FALSE]
			for (thresh in thresholds){
				X_cov = drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M
				X_try = matrix(w, ncol = 1)
				colnames(X_try) = "w"
				if (ncol(X_cov) > 0){
					X_try = cbind(X_try, X_cov)
				}
				X_try_df = as.data.frame(X_try)
				key = paste(colnames(X_try_df), collapse = "|")
				if (key %in% keys) next
				keys = c(keys, key)
				attempt_try = private$fit_with_hardened_qr_column_dropping(
					X_full = X_try,
					required_cols = 1L,
					fit_fun = function(X_fit) X_fit,
					fit_ok = function(mod, X_fit, keep) TRUE
				)
				candidates[[length(candidates) + 1L]] = attempt_try$X
			}
			candidates
		},
		rcpp_cox_fit_is_usable = function(fit, estimate_only = FALSE){
			if (is.null(fit) || !isTRUE(fit$converged)) return(FALSE)
			coef1 = fit$coefficients[1L]
			if (!is.finite(coef1) || abs(coef1) > private$max_abs_reasonable_coef) return(FALSE)
			if (estimate_only) return(TRUE)
			se = sqrt(fit$vcov[1L, 1L])
			is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(invisible(NULL))
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (sum(private$dead) == 0L){
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_no_events")
				return(invisible(NULL))
			}
			# ── Matched pairs ───────────────────────────────────────────────
			beta_m = NA_real_; ssq_m = NA_real_
			if (m > 0){
				y_m = KKstats$y_matched_long
				dead_m = KKstats$dead_matched_long
				w_m = KKstats$w_matched_long
				strata_m = KKstats$m_matched_long
				X_cov_m = KKstats$X_matched_long
				X_m = cbind(w = w_m, X_cov_m)
				
				res = tryCatch(
					fast_stratified_coxph_regression_cpp(X_m, y_m, dead_m, as.integer(strata_m), estimate_only = estimate_only),
					error = function(e) NULL
				)
				if (private$rcpp_cox_fit_is_usable(res, estimate_only = estimate_only)){
					beta_m = res$coefficients[1L]
					if (!estimate_only) ssq_m = res$vcov[1L, 1L]
				}
			}
			# ── Reservoir ───────────────────────────────────────────────────
			beta_r = NA_real_; ssq_r = NA_real_
			res = NULL
			if (nRT > 0 && nRC > 0){
				y_r = KKstats$y_reservoir
				dead_r = KKstats$dead_reservoir
				w_r = KKstats$w_reservoir
				X_cov_r = as.matrix(KKstats$X_reservoir)
				candidates = private$cox_design_candidates(w_r, X_cov_r)
				for (X_candidate in candidates){
					res_try = tryCatch(
						fast_coxph_regression_cpp(X_candidate, y_r, dead_r, estimate_only = estimate_only),
						error = function(e) NULL
					)
					if (private$rcpp_cox_fit_is_usable(res_try, estimate_only = estimate_only)){
						res = res_try
						break
					}
				}
			}
			if (is.null(res)){
				X_mat = matrix(private$w[private$m == 0], ncol = 1L)
				colnames(X_mat) = "w"
				res = tryCatch(
					fast_coxph_regression_cpp(X_mat, y_r, dead_r, estimate_only = estimate_only),
					error = function(e) NULL
				)
			}
			if (!is.null(res)){
				beta_r = res$coefficients[1L]
				if (!estimate_only) ssq_r = res$vcov[1L, 1L]
			}
			# ── Combine ─────────────────────────────────────────────────────
			m_ok = is.finite(beta_m) && (estimate_only || (is.finite(ssq_m) && ssq_m > 0))
			r_ok = is.finite(beta_r) && (estimate_only || (is.finite(ssq_r) && ssq_r > 0))
			if (m_ok && r_ok){
				if (estimate_only){
					# Use simple sample size weighting as fallback if SEs not requested
					w_star = m / (m + (nRT+nRC)/2)
					private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
				} else {
					w_star = ssq_r / (ssq_r + ssq_m)
					private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
					private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
				}
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_both_failed")
			}
		}
	)))
)

#' Stratified Cox Combined-Likelihood Compound Inference for KK Designs
#'
#' @export
InferenceSurvivalKKStratCoxPHOneLik = R6::R6Class("InferenceSurvivalKKStratCoxPHOneLik",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed KK survival design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only If \code{TRUE}, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the one-likelihood stratified Cox treatment
		#'   estimate under Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$beta_hat_T = beta_hat_T
					private$cached_values$s_beta_hat_T = NA_real_
					return(private$cached_values$beta_hat_T)
				}
			}
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			strata = m_vec
			res_idx = which(strata == 0L)
			if (length(res_idx) > 0L) {
				max_m = max(strata)
				strata[res_idx] = max_m + seq_along(res_idx)
			}
			X_cov = private$get_X()
			X_fit = if (ncol(as.matrix(X_cov)) > 0) cbind(treatment = private$w, X_cov) else matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			fit = weighted_cox_bootstrap_surrogate_fit(
				private$y, private$dead, X_fit, row_weights,
				strata = strata,
				warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(X_fit)) %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit))
			)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an asymptotic confidence interval for the treatment effect.
		#' @param alpha Confidence level.
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
		#' @description Returns a 2-sided p-value for H0: beta_T = delta.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		max_abs_reasonable_coef = 1e4,
		best_X_colnames = NULL,
		optimization_alg = "lbfgs",
		design_matrix_candidates = function(){
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			X_covs = private$get_X()
			if (ncol(as.matrix(X_covs)) > 0L){
				X_full = cbind(X_full, as.matrix(X_covs))
			}
			X_full
		},
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			strata_ids = m_vec
			res_idx = which(strata_ids == 0L)
			if (length(res_idx) > 0L){
				max_m = if (any(strata_ids > 0L)) max(strata_ids) else 0L
				strata_ids[res_idx] = max_m + 1L
			}
			X_full = private$design_matrix_candidates()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					fast_stratified_coxph_regression_cpp(
						X = as.matrix(X_fit),
						y = private$y,
						dead = private$dead,
						strata = as.integer(strata_ids),
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					)
				},
				fit_ok = function(res, X_fit, keep){
					if (is.null(res) || !isTRUE(res$converged)) return(FALSE)
					beta = res$coefficients[1L]
					if (!is.finite(beta) || abs(beta) > private$max_abs_reasonable_coef) return(FALSE)
					if (estimate_only) return(TRUE)
					se = tryCatch(sqrt(res$vcov[1L, 1L]), error = function(e) NA_real_)
					is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
				}
			)
			res = attempt$fit
			if (!is.null(res)){
				private$best_X_colnames = setdiff(colnames(attempt$X), "w")
				private$cached_mod = res
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					y = private$y,
					dead = private$dead,
					strata = strata_ids,
					j_treat = 1L
				)
				private$cached_values$beta_hat_T = as.numeric(res$coefficients[1L])
				if (!estimate_only){
					se = sqrt(res$vcov[1L, 1L])
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
				return(invisible(NULL))
			}
			private$cache_nonestimable_estimate("kk_strat_cox_combined_fit_failed")
			invisible(NULL)
		},
		supports_likelihood_tests = function() TRUE,
		get_likelihood_test_spec = function(){
			private$shared_combined_likelihood(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			dead = as.numeric(ctx$dead)
			strata = as.integer(ctx$strata)
			j_treat = as.integer(ctx$j_treat %||% 1L)
			list(
				X = X_fit, y = y, dead = dead, strata = strata, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm = start %||% private$get_fit_warm_start_for_length("params", ncol(X_fit))
					fast_stratified_coxph_regression_cpp(
						X = X_fit, y = y, dead = dead, strata = strata,
						warm_start_beta = warm,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg,
						fixed_idx = j_treat, fixed_values = delta
					)
				},
				extract_start = function(fit) as.numeric(fit$coefficients %||% fit$b),
				score = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					get_stratified_coxph_score_cpp(X_fit, y, dead, strata, beta)
				},
				observed_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					-get_stratified_coxph_hessian_cpp(X_fit, y, dead, strata, beta)
				},
				fisher_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					-get_stratified_coxph_hessian_cpp(X_fit, y, dead, strata, beta)
				},
				information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					-get_stratified_coxph_hessian_cpp(X_fit, y, dead, strata, beta)
				},
				neg_loglik = function(fit) as.numeric(fit$neg_ll %||% fit$neg_loglik)
			)
		},
		get_standard_error = function(){
			private$shared_combined_likelihood()
			as.numeric(private$cached_values$s_beta_hat_T)
		},
		get_degrees_of_freedom = function() Inf,
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$coefficients %||% null_fit$b)
			if (!all(is.finite(b_null))) return(NULL)
			X_fit = spec$X
			y_obs = as.numeric(spec$y)
			dead_obs = as.numeric(spec$dead)
			strata = as.integer(spec$strata)
			j = spec$j
			sim = .cox_simulate_stratified(y_obs, dead_obs, X_fit, b_null, strata)
			y_sim = sim$y_sim; dead_sim = sim$dead_sim
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			full_res = tryCatch(
				fast_stratified_coxph_regression_cpp(
					X = X_fit, y = y_sim, dead = dead_sim, strata = strata,
					estimate_only = FALSE,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$coefficients[j])) return(NULL)
			full_fit_boot = list(b = as.numeric(full_res$coefficients), neg_loglik = as.numeric(full_res$neg_ll))
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = tryCatch(
						fast_stratified_coxph_regression_cpp(
							X = X_fit, y = y_sim, dead = dead_sim, strata = strata,
							warm_start_beta = start %||% as.numeric(full_res$coefficients),
							fixed_idx = j, fixed_values = d,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), neg_loglik = as.numeric(res$neg_ll))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik %||% fit$neg_ll)
			)
		}
	)))
)
