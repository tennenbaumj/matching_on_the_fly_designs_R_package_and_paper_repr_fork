#' Abstract class for Weibull Frailty / Standard Weibull Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Weibull AFT GLMM for matched pairs.
#' The matched-pair component uses a shared log-normal random intercept per pair,
#' fitted by the package's native Rcpp likelihood optimizer. The reservoir component uses standard Weibull AFT regression.
#' The two treatment-effect estimates are combined by inverse-variance weighting.
#'
#' @details
#' This compound estimator accounts for the dependence within matched pairs by
#' modeling it as a shared frailty.
#'
#' \strong{Univariate} (\code{ncol(as.matrix(private$X)) == 0}): uses the native
#' Rcpp Weibull frailty likelihood with \code{formula = survival::Surv(y, dead) ~ w}
#' and a pair-level random intercept.
#'
#' \strong{Multivariate} (\code{ncol(as.matrix(private$X)) > 0}): fits the same
#' native Rcpp likelihood with covariate adjustment, dropping rank-deficient
#' columns when needed.
#'
#' @keywords internal
InferenceAbstractKKWeibullFrailtyIVWC = R6::R6Class("InferenceAbstractKKWeibullFrailtyIVWC",
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
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the estimated treatment effect (log-time ratio).
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
		#' @param delta                                   The null difference to test against. Default is 0.
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
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		},
		#' @description Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			inf_obj = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			inf_obj
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		optimization_alg = NULL,
		best_par = NULL,
		best_X_colnames = NULL,
		any_censoring = NULL,
		m = NULL,
		cached_mod = NULL,
		best_X_colnames_matched = NULL,
		best_X_colnames_reservoir = NULL,
		max_abs_reasonable_coef = 1e4,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		supports_likelihood_tests = function() FALSE,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Use private fields (survive duplicate()) instead of cached_values
			if (is.null(private$best_X_colnames_matched) && is.null(private$best_X_colnames_reservoir)){
				private$shared()
			}
			if (is.null(private$best_X_colnames_matched) && is.null(private$best_X_colnames_reservoir)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(NA_real_)
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			X_data = private$get_X()
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			# Matched pairs component
			beta_m = NA_real_
			ssq_m = NA_real_
			if (m > 0 && !is.null(private$best_X_colnames_matched)){
				i_matched = which(m_vec > 0L)
				X_cov = X_data[i_matched, intersect(private$best_X_colnames_matched, colnames(X_data)), drop = FALSE]
				X_m = cbind(w = private$w[i_matched], X_cov)
				# Fixed-VC fast path
				if (!is.null(private$cached_vc_params_matched) && all(is.finite(private$cached_vc_params_matched))) {
					p_m = ncol(X_m)
					group_m = m_vec[i_matched]
					fit_fast_m = tryCatch(
						fast_weibull_frailty_cpp(
							X = as.matrix(X_m),
							y = private$y[i_matched], dead = private$dead[i_matched],
							group_id = as.integer(group_m),
							estimate_only = TRUE,
							optimization_alg = private$optimization_alg,
							fixed_idx    = as.integer(c(p_m + 1L, p_m + 2L)),
							fixed_values = as.numeric(private$cached_vc_params_matched)
						),
						error = function(e) NULL
					)
					if (!is.null(fit_fast_m) && isTRUE(fit_fast_m$converged) && length(fit_fast_m$b) >= 1L && is.finite(fit_fast_m$b[1L]))
						beta_m = as.numeric(fit_fast_m$b[1L])
				}
				if (!is.finite(beta_m)) {
					fit_m = .fit_weibull_frailty(
						y = private$y[i_matched],
						dead = private$dead[i_matched],
						X = X_m,
						pair_id = m_vec[i_matched],
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					)
					if (!is.null(fit_m) && is.finite(fit_m$beta)){
						beta_m = fit_m$beta
						if (!estimate_only && !is.null(fit_m$ssq) && is.finite(fit_m$ssq) && fit_m$ssq > 0)
							ssq_m = fit_m$ssq
					}
				}
			}
			# Reservoir component
			beta_r = NA_real_
			ssq_r = NA_real_
			if (nRT > 0 && nRC > 0 && !is.null(private$best_X_colnames_reservoir)){
				i_reservoir = which(m_vec == 0L)
				X_cov_r = X_data[i_reservoir, intersect(private$best_X_colnames_reservoir, colnames(X_data)), drop = FALSE]
				X_r = cbind(w = private$w[i_reservoir], X_cov_r)
				# Fixed-VC fast path (Weibull regression, X must include intercept)
				if (!is.null(private$cached_vc_params_reservoir) && is.finite(private$cached_vc_params_reservoir[1L])) {
					X_r_int = cbind("(Intercept)" = 1, as.matrix(X_r))
					p_r = ncol(X_r_int)
					fit_fast_r = tryCatch(
						fast_weibull_regression_cpp(
							y    = as.numeric(private$y[i_reservoir]),
							dead = as.numeric(private$dead[i_reservoir]),
							X    = X_r_int,
							estimate_only = TRUE,
							fixed_idx    = as.integer(p_r),
							fixed_values = as.numeric(private$cached_vc_params_reservoir[1L])
						),
						error = function(e) NULL
					)
					if (!is.null(fit_fast_r) && isTRUE(fit_fast_r$converged) && length(fit_fast_r$b) >= 2L && is.finite(fit_fast_r$b[2L]))
						beta_r = as.numeric(fit_fast_r$b[2L])
				}
				if (!is.finite(beta_r)) {
					fit_r = .fit_standard_weibull_aft_from_matrix(
						y = private$y[i_reservoir],
						dead = private$dead[i_reservoir],
						X = X_r,
						estimate_only = estimate_only
					)
					if (!is.null(fit_r) && is.finite(fit_r$beta)){
						beta_r = fit_r$beta
						if (!estimate_only && !is.null(fit_r$ssq) && is.finite(fit_r$ssq) && fit_r$ssq > 0)
							ssq_r = fit_r$ssq
					}
				}
			}
			# Pooling
			m_ok = is.finite(beta_m) && (estimate_only || is.finite(ssq_m))
			r_ok = is.finite(beta_r) && (estimate_only || is.finite(ssq_r))
			if (m_ok && r_ok){
				if (estimate_only) {
					ssq_m_orig = private$cached_values$ssq_beta_T_matched
					ssq_r_orig = private$cached_values$ssq_beta_T_reservoir
					if (!is.null(ssq_m_orig) && !is.null(ssq_r_orig) && is.finite(ssq_m_orig) && is.finite(ssq_r_orig)){
						w_star = ssq_r_orig / (ssq_r_orig + ssq_m_orig)
						return(w_star * beta_m + (1 - w_star) * beta_r)
					}
					return(0.5 * beta_m + 0.5 * beta_r)
				}
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		},
		best_X_colnames_matched = NULL,
		best_X_colnames_reservoir = NULL,
		cached_vc_params_matched = NULL,
		cached_vc_params_reservoir = NULL,
		max_abs_reasonable_coef = 1e4,
		frailty_for_matched_pairs = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			i_matched = which(m_vec > 0L)
			if (length(i_matched) == 0L) return(invisible(NULL))
			X_full = if (ncol(as.matrix(private$X)) == 0L){
				matrix(private$w[i_matched], ncol = 1L, dimnames = list(NULL, "w"))
			} else {
				cbind(w = private$w[i_matched], private$get_X()[i_matched, , drop = FALSE])
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				fit_fun = function(X_fit, keep){
					res = .fit_weibull_frailty(
						y = private$y[i_matched],
						dead = private$dead[i_matched],
						X = X_fit,
						pair_id = m_vec[i_matched],
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					)
					res
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || !is.finite(mod$beta)) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq) && mod$ssq > 0
				}
			)
			
			if (!is.null(attempt$fit)){
				private$cached_values$beta_T_matched = attempt$fit$beta
				private$cached_values$ssq_beta_T_matched = attempt$fit$ssq
				best_cols = setdiff(colnames(attempt$X), "w")
				private$cached_values$best_X_colnames_matched = best_cols
				private$best_X_colnames_matched = best_cols
				if (!is.null(attempt$fit$log_sigma_eps) && !is.null(attempt$fit$log_sigma_u))
					private$cached_vc_params_matched = as.numeric(c(attempt$fit$log_sigma_eps, attempt$fit$log_sigma_u))
			}
		},
		weibull_for_reservoir = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			i_reservoir = which(m_vec == 0L)
			if (length(i_reservoir) == 0L) return(invisible(NULL))
			X_full = if (ncol(as.matrix(private$X)) == 0L){
				matrix(private$w[i_reservoir], ncol = 1L, dimnames = list(NULL, "w"))
			} else {
				cbind(w = private$w[i_reservoir], private$get_X()[i_reservoir, , drop = FALSE])
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				fit_fun = function(X_fit, keep){
					res = .fit_standard_weibull_aft_from_matrix(
						y = private$y[i_reservoir],
						dead = private$dead[i_reservoir],
						X = X_fit,
						estimate_only = estimate_only
					)
					res
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || !is.finite(mod$beta)) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq) && mod$ssq > 0
				}
			)
			
			if (!is.null(attempt$fit)){
				private$cached_values$beta_T_reservoir = attempt$fit$beta
				private$cached_values$ssq_beta_T_reservoir = attempt$fit$ssq
				best_cols_r = setdiff(colnames(attempt$X), "w")
				private$cached_values$best_X_colnames_reservoir = best_cols_r
				private$best_X_colnames_reservoir = best_cols_r
				X_r_int = cbind("(Intercept)" = 1, as.matrix(attempt$X))
				res_log_s = tryCatch(
					fast_weibull_regression_cpp(
						y    = as.numeric(private$y[i_reservoir]),
						dead = as.numeric(private$dead[i_reservoir]),
						X    = X_r_int,
						estimate_only = TRUE
					),
					error = function(e) NULL
				)
				if (!is.null(res_log_s) && isTRUE(res_log_s$converged))
					private$cached_vc_params_reservoir = as.numeric(res_log_s$log_sigma)
			}
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(invisible(NULL))
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (m > 0){
				private$frailty_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
			       (!estimate_only && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0 || estimate_only)
			if (nRT > 0 && nRC > 0){
				private$weibull_for_reservoir(estimate_only = estimate_only)
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) &&
			       (!estimate_only && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0 || estimate_only)
			if (m_ok && r_ok){
				if (estimate_only){
					if (!is.null(ssq_m) && !is.null(ssq_r) && is.finite(ssq_m) && is.finite(ssq_r)){
						w_star = ssq_r / (ssq_r + ssq_m)
						private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
					} else {
						private$cached_values$beta_hat_T = 0.5 * beta_m + 0.5 * beta_r
					}
					return(invisible(NULL))
				}
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)))
)
#' Abstract class for Weibull Frailty Combined-Likelihood Inference
#'
#' @keywords internal
InferenceAbstractKKWeibullFrailtyOneLik = R6::R6Class("InferenceAbstractKKWeibullFrailtyOneLik",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed KK survival design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Logical. If \code{TRUE}, use the internal Rcpp backend.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg Optimization algorithm to use.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			private$use_rcpp = use_rcpp
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only If \code{TRUE}, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the one-likelihood Weibull-frailty treatment
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
			X_cov = private$get_X()
			X_fit = if (ncol(as.matrix(X_cov)) > 0) cbind(treatment = private$w, X_cov) else matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			cluster_ids = m_vec
			res_idx = which(cluster_ids == 0L)
			if (length(res_idx) > 0L) {
				max_m = max(cluster_ids)
				cluster_ids[res_idx] = max_m + seq_along(res_idx)
			}
			fit = weighted_weibull_bootstrap_surrogate_fit(private$y, private$dead, X_fit, row_weights, cluster = cluster_ids)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an asymptotic confidence interval for the treatment effect.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Returns a 2-sided p-value for H0: beta_T = delta.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
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
		cached_mod = NULL,
		best_X_colnames = NULL,
		optimization_alg = "lbfgs",
		best_par = NULL,
		any_censoring = NULL,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		use_rcpp = TRUE,
		cached_vc_params = NULL,
		max_abs_reasonable_coef = 1e4,
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			if (sum(private$dead) == 0L){
				private$cache_nonestimable_estimate("kk_weibull_frailty_no_events")
				return(invisible(NULL))
			}
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			res_idx = which(group_id == 0L)
			if (length(res_idx) > 0L){
				max_m = if (any(group_id > 0L)) max(group_id) else 0L
				group_id[res_idx] = max_m + seq_along(res_idx)
			}
			X_full = if (ncol(as.matrix(private$X)) == 0L){
				matrix(private$w, ncol = 1L, dimnames = list(NULL, "w"))
			} else {
				cbind(w = private$w, private$get_X())
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					fast_weibull_frailty_cpp(
						X = as.matrix(X_fit),
						y = private$y,
						dead = private$dead,
						group_id = as.integer(group_id),
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					)
				},
				fit_ok = function(res, X_fit, keep){
					if (is.null(res) || !isTRUE(res$converged)) return(FALSE)
					beta = as.numeric(res$b[1L])
					if (!is.finite(beta) || abs(beta) > private$max_abs_reasonable_coef) return(FALSE)
					if (estimate_only) return(TRUE)
					se = tryCatch(sqrt(as.numeric(res$ssq_b_T)), error = function(e) NA_real_)
					is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
				}
			)
			res = attempt$fit
			if (!is.null(res)){
				private$cached_values$best_X_colnames = setdiff(colnames(attempt$X), "w")
				private$cached_mod = res
				if (!is.null(res$log_sigma_eps) && !is.null(res$log_sigma_u))
					private$cached_vc_params = as.numeric(c(res$log_sigma_eps, res$log_sigma_u))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					y = private$y,
					dead = private$dead,
					group_id = group_id,
					j_treat = 1L
				)
				private$cached_values$beta_hat_T = as.numeric(res$b[1L])
				if (!estimate_only){
					se = sqrt(as.numeric(res$ssq_b_T))
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
				return(invisible(NULL))
			}
			private$cache_nonestimable_estimate("kk_weibull_frailty_combined_fit_failed")
			invisible(NULL)
		},
		supports_likelihood_tests = function() FALSE,
		get_likelihood_test_spec = function(){
			private$shared_combined_likelihood(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			if (!is.finite(as.numeric(private$cached_mod$neg_loglik %||% NA_real_))) return(NULL)
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			dead = as.numeric(ctx$dead)
			group_id = as.integer(ctx$group_id)
			j_treat = as.integer(ctx$j_treat %||% 1L)
			p = ncol(X_fit)
			list(
				X = X_fit, y = y, dead = dead, j = j_treat,
				group_id = group_id,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm = start %||% private$get_fit_warm_start_for_length("params", p + 2L)
					fast_weibull_frailty_cpp(
						X = X_fit, y = y, dead = dead,
						group_id = group_id,
						warm_start_params = warm,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg,
						fixed_idx = j_treat, fixed_values = delta
					)
				},
				extract_start = function(fit){
					c(as.numeric(fit$b), as.numeric(fit$log_sigma_eps), as.numeric(fit$log_sigma_u))
				},
				score = function(fit) rep(NA_real_, p + 2L),
				observed_information = function(fit) matrix(NA_real_, p + 2L, p + 2L),
				fisher_information = function(fit) matrix(NA_real_, p + 2L, p + 2L),
				information = function(fit) matrix(NA_real_, p + 2L, p + 2L),
				neg_loglik = function(fit) as.numeric(fit$neg_loglik %||% fit$neg_ll)
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
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		simulate_under_lik_null = function(spec, delta, null_fit){
			p = ncol(spec$X)
			b_null = as.numeric(null_fit$b)
			if (length(b_null) < p || !all(is.finite(b_null))) return(NULL)
			log_sigma_eps = as.numeric(null_fit$log_sigma_eps)
			log_sigma_u = as.numeric(null_fit$log_sigma_u)
			if (!is.finite(log_sigma_eps) || !is.finite(log_sigma_u)) return(NULL)
			sigma_eps = exp(log_sigma_eps)
			sigma_u = exp(log_sigma_u)
			X = spec$X
			group_id = spec$group_id
			n = nrow(X)
			K = max(group_id)
			u_g = rnorm(K, 0, sigma_u)
			mu = as.numeric(X %*% b_null) + u_g[group_id]
			T_sim = rweibull(n, shape = 1 / sigma_eps, scale = exp(mu))
			if (!all(is.finite(T_sim)) || any(T_sim <= 0)) return(NULL)
			y_obs = as.numeric(spec$y)
			dead_obs = as.numeric(spec$dead)
			C_i = ifelse(dead_obs == 0, y_obs, Inf)
			y_sim = pmin(T_sim, C_i)
			dead_sim = as.numeric(T_sim <= C_i)
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			j = spec$j
			full_res = tryCatch(
				fast_weibull_frailty_cpp(
					X = X, y = y_sim, dead = dead_sim,
					group_id = group_id,
					estimate_only = FALSE,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || length(full_res$b) < j || !is.finite(full_res$b[j]) || !is.finite(as.numeric(full_res$neg_loglik %||% NA_real_))) return(NULL)
			list(
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					warm = start %||% c(as.numeric(full_res$b), as.numeric(full_res$log_sigma_eps), as.numeric(full_res$log_sigma_u))
					tryCatch(
						fast_weibull_frailty_cpp(
							X = X, y = y_sim, dead = dead_sim,
							group_id = group_id,
							warm_start_params = warm,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j, fixed_values = d
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik %||% fit$neg_ll)
			)
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			private$compute_basic_match_data()
			if (is.null(private$cached_values$best_X_colnames)){
				private$shared_combined_likelihood(estimate_only = TRUE)
			}
			if (is.null(private$cached_values$best_X_colnames)) return(NA_real_)
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			res_idx = which(group_id == 0L)
			if (length(res_idx) > 0L){
				max_m = if (any(group_id > 0L)) max(group_id) else 0L
				group_id[res_idx] = max_m + seq_along(res_idx)
			}
			X_data = private$get_X()
			X_full = matrix(private$w, ncol = 1L, dimnames = list(NULL, "w"))
			X_covs = X_data[, intersect(private$cached_values$best_X_colnames, colnames(X_data)), drop = FALSE]
			if (ncol(X_covs) > 0L) X_full = cbind(X_full, X_covs)
			p_ncol = ncol(X_full)
			if (!is.null(private$cached_vc_params) && all(is.finite(private$cached_vc_params))) {
				fit_fast = tryCatch(
					fast_weibull_frailty_cpp(
						X = as.matrix(X_full),
						y = private$y, dead = private$dead,
						group_id = as.integer(group_id),
						estimate_only = TRUE,
						optimization_alg = private$optimization_alg,
						fixed_idx    = as.integer(c(p_ncol + 1L, p_ncol + 2L)),
						fixed_values = as.numeric(private$cached_vc_params)
					),
					error = function(e) NULL
				)
				if (!is.null(fit_fast) && isTRUE(fit_fast$converged) && length(fit_fast$b) >= 1L && is.finite(fit_fast$b[1L]))
					return(as.numeric(fit_fast$b[1L]))
			}
			fit = tryCatch(
				fast_weibull_frailty_cpp(
					X = as.matrix(X_full),
					y = private$y, dead = private$dead,
					group_id = as.integer(group_id),
					estimate_only = estimate_only,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged) || length(fit$b) < 1L || !is.finite(fit$b[1L])) return(NA_real_)
			as.numeric(fit$b[1L])
		}
	)))
)
#' Weibull Frailty IVWC Inference for KK Designs
#' @export
InferenceSurvivalKKWeibullFrailtyIVWC = R6::R6Class("InferenceSurvivalKKWeibullFrailtyIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKWeibullFrailtyIVWC,
	public = list(
		#' @description Initialize the IVWC Weibull-frailty inference object.
		#' @param des_obj A completed KK survival design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg Optimization algorithm to use.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, optimization_alg = NULL){
			self$set_optimization_alg(optimization_alg)
			super$initialize(des_obj = des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)
#' Weibull Frailty Combined-Likelihood Inference for KK Designs
#' @export
InferenceSurvivalKKWeibullFrailtyOneLik = R6::R6Class("InferenceSurvivalKKWeibullFrailtyOneLik",
	inherit = InferenceAbstractKKWeibullFrailtyOneLik,
	public = list(
		#' @description Initialize the one-likelihood Weibull-frailty inference object.
		#' @param des_obj A completed KK survival design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Logical. If \code{TRUE}, use the internal Rcpp backend.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg Optimization algorithm to use.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL){
			self$set_optimization_alg(optimization_alg)
			super$initialize(des_obj = des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose)
		}
	)
)
