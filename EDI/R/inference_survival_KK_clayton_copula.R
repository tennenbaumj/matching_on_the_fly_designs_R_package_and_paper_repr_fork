#' Clayton Copula / Standard Weibull Compound Inference for KK Designs
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Clayton copula with Weibull AFT margins for matched
#' pairs and a standard Weibull AFT model for the reservoir. The two treatment-effect
#' estimates (on the log-time ratio scale) are combined by inverse-variance weighting.
#'
#' @export
InferenceSurvivalKKClaytonCopulaIVWC = R6::R6Class("InferenceSurvivalKKClaytonCopulaIVWC",
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
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the estimated treatment effect (log-time ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the IVWC Clayton-copula treatment estimate under
		#'   Bayesian-bootstrap weights.
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
			fit = weighted_weibull_bootstrap_surrogate_fit(private$y, private$dead, X_fit, row_weights)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
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
		optimization_alg = "lbfgs",
		best_par = NULL,
		best_X_colnames = NULL,
		cached_mod = NULL,
		best_X_colnames_matched = NULL,
		best_X_colnames_reservoir = NULL,
		cached_vc_params_matched = NULL,
		cached_vc_params_reservoir = NULL,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		max_abs_reasonable_coef = 1e4,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			private$compute_basic_match_data()
			if (is.null(private$best_X_colnames_matched) && is.null(private$best_X_colnames_reservoir)){
				private$shared()
			}
			if (is.null(private$best_X_colnames_matched) && is.null(private$best_X_colnames_reservoir)){
				return(NA_real_)
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
			# Matched component
			beta_m = NA_real_
			if (m > 0 && !is.null(private$best_X_colnames_matched)) {
				i_matched = which(m_vec > 0L)
				cov_cols_m = setdiff(private$best_X_colnames_matched, "w")
				X_cand_m = cbind(w = private$w[i_matched], X_data[i_matched, intersect(cov_cols_m, colnames(X_data)), drop = FALSE])
				if (!is.null(private$cached_vc_params_matched) && all(is.finite(private$cached_vc_params_matched))) {
					X_m_int = cbind("(Intercept)" = 1, X_cand_m)
					p_m = ncol(X_m_int)
					pair_idx_m_r = .complete_pair_index_matrix(m_vec[i_matched])
					pair_idx_m_0 = if (nrow(pair_idx_m_r) > 0L) pair_idx_m_r - 1L else matrix(0L, 0L, 2L)
					warm_m = private$get_fit_warm_start_for_length("params", p_m + 2L)
					if (is.null(warm_m) || length(warm_m) != p_m + 2L)
						warm_m = c(rep(0, p_m), private$cached_vc_params_matched)
					fit_fast_m = tryCatch(
						fast_clayton_weibull_aft_optim_cpp(
							X = as.matrix(X_m_int),
							y = as.numeric(private$y[i_matched]),
							dead = as.numeric(private$dead[i_matched]),
							pair_idx = pair_idx_m_0,
							singleton_rows = integer(0),
							warm_start_params = warm_m,
							estimate_only = TRUE,
							optimization_alg = private$optimization_alg,
							fixed_idx = as.integer(c(p_m + 1L, p_m + 2L)),
							fixed_values = as.numeric(private$cached_vc_params_matched)
						), error = function(e) NULL
					)
					if (!is.null(fit_fast_m) && isTRUE(fit_fast_m$converged) && length(fit_fast_m$params) >= 2L)
						beta_m = as.numeric(fit_fast_m$params[2L])
				}
				if (!is.finite(beta_m)) {
					fit_m = .fit_clayton_weibull_aft(
						y = private$y[i_matched], dead = private$dead[i_matched],
						X = X_cand_m, pair_id = m_vec[i_matched], estimate_only = TRUE,
						optimization_alg = private$optimization_alg
					)
					if (!is.null(fit_m) && is.finite(fit_m$beta)) beta_m = fit_m$beta
				}
			}
			# Reservoir component
			beta_r = NA_real_
			if (nRT > 0 && nRC > 0 && !is.null(private$best_X_colnames_reservoir)) {
				i_reservoir = which(m_vec == 0L)
				y_r = private$y[i_reservoir]
				w_r = private$w[i_reservoir]
				dead_r = private$dead[i_reservoir]
				cov_cols_r = setdiff(private$best_X_colnames_reservoir, "w")
				X_cov_r = X_data[i_reservoir, intersect(cov_cols_r, colnames(X_data)), drop = FALSE]
				X_r = cbind(w = w_r, X_cov_r)
				if (!is.null(private$cached_vc_params_reservoir) && is.finite(private$cached_vc_params_reservoir[1L])) {
					X_r_int = cbind("(Intercept)" = 1, as.matrix(X_r))
					p_r = ncol(X_r_int)
					fit_fast_r = tryCatch(
						fast_weibull_regression_cpp(
							y = as.numeric(y_r), dead = as.numeric(dead_r), X = X_r_int,
							estimate_only = TRUE,
							fixed_idx = as.integer(p_r),
							fixed_values = as.numeric(private$cached_vc_params_reservoir[1L])
						), error = function(e) NULL
					)
					if (!is.null(fit_fast_r) && isTRUE(fit_fast_r$converged) && length(fit_fast_r$b) >= 2L)
						beta_r = as.numeric(fit_fast_r$b[2L])
				}
				if (!is.finite(beta_r)) {
					fit_r = .fit_standard_weibull_aft_from_matrix(
						y = y_r, dead = dead_r, X = X_r, estimate_only = TRUE
					)
					if (!is.null(fit_r) && is.finite(fit_r$beta)) beta_r = fit_r$beta
				}
			}
			# Pooling
			m_ok = is.finite(beta_m)
			r_ok = is.finite(beta_r)
			if (m_ok && r_ok) {
				ssq_m_orig = private$cached_values$ssq_beta_T_matched
				ssq_r_orig = private$cached_values$ssq_beta_T_reservoir
				if (!is.null(ssq_m_orig) && !is.null(ssq_r_orig) && is.finite(ssq_m_orig) && is.finite(ssq_r_orig)) {
					w_star = ssq_r_orig / (ssq_r_orig + ssq_m_orig)
					return(w_star * beta_m + (1 - w_star) * beta_r)
				}
				return(0.5 * beta_m + 0.5 * beta_r)
			} else if (m_ok) {
				return(beta_m)
			} else if (r_ok) {
				return(beta_r)
			}
			NA_real_
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		filtered_covariate_candidates = function(X = as.matrix(private$X)){
			if (ncol(X) == 0L) return(list(matrix(nrow = nrow(X), ncol = 0L)))
			
			# Ensure no linearly dependent columns first
			res = drop_linearly_dependent_cols(X)
			X_reduced = res$M
			if (ncol(X_reduced) == 0L) return(list(matrix(nrow = nrow(X), ncol = 0L)))
			
			# Generate candidates by dropping highly correlated columns at various thresholds
			thresholds = c(0.99, 0.9, 0.7)
			candidates = list(X_reduced)
			for (thresh in thresholds) {
				X_try = drop_highly_correlated_cols(X_reduced, threshold = thresh)$M
				# Simple way to avoid duplicates in the list
				is_new = TRUE
				for (c in candidates) {
					if (ncol(c) == ncol(X_try) && all(colnames(c) == colnames(X_try))) {
						is_new = FALSE
						break
					}
				}
				if (is_new) candidates[[length(candidates) + 1L]] = X_try
			}
			candidates
		},
		design_matrix_candidates = function(){
			if (!is.null(private$cached_values$clayton_design_candidates)){
				return(private$cached_values$clayton_design_candidates)
			}
			candidates = list()
			# Case 1: no covariates
			if (ncol(as.matrix(private$X)) == 0L){
				Xcand = matrix(private$w, ncol = 1)
				colnames(Xcand) = "w"
				candidates = list(Xcand)
			} else {
				cov_candidates = private$filtered_covariate_candidates()
				for (X in cov_candidates){
					M = matrix(private$w, ncol = 1)
					colnames(M) = "w"
					if (ncol(X) > 0){
						M = cbind(M, X)
					}
					# Ensure we drop any additional linear dependencies
					qr_res = qr(M)
					if (qr_res$rank < ncol(M)){
						keep = qr_res$pivot[seq_len(qr_res$rank)]
						if (!(1L %in% keep)) keep = c(1L, keep) # Keep treatment
						keep = sort(unique(keep))
						M = M[, keep, drop = FALSE]
					}
					candidates[[length(candidates) + 1L]] = M
				}
			}
			private$cached_values$clayton_design_candidates = candidates
			candidates
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (m > 0){
				private$clayton_copula_for_matched_pairs(estimate_only = estimate_only)
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
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
		},
		clayton_copula_for_matched_pairs = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			i_matched = which(m_vec > 0L)
			if (length(i_matched) == 0L) return(invisible(NULL))
			X_data = private$get_X()
			X_matched = X_data[i_matched, , drop = FALSE]
			X_fit = cbind(w = private$w[i_matched], X_matched)
			n_params = ncol(X_fit) + 2L
			# Fit using optimized helper
			fit = .fit_clayton_weibull_aft(
				y = private$y[i_matched],
				dead = private$dead[i_matched],
				X = X_fit,
				pair_id = m_vec[i_matched],
				estimate_only = estimate_only,
				optimization_alg = private$optimization_alg,
				starts = if (!is.null(private$get_fit_warm_start_for_length("params", n_params))) 
				             list(private$get_fit_warm_start_for_length("params", n_params)) else NULL,
				warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
			)
			
			if (!is.null(fit) && is.finite(fit$beta)){
				private$cached_values$beta_T_matched = fit$beta
				private$cached_values$ssq_beta_T_matched = if (estimate_only) NA_real_ else fit$ssq
				private$best_X_colnames_matched = colnames(X_fit)
				if (!is.null(fit$best_par) && length(fit$best_par) >= 3L)
					private$cached_vc_params_matched = as.numeric(tail(fit$best_par, 2L))
			}
		},
		weibull_for_reservoir = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			y_r    = KKstats$y_reservoir
			w_r    = KKstats$w_reservoir
			m_vec_safe = private$m
			if (is.null(m_vec_safe)) m_vec_safe = rep(0L, private$n)
			m_vec_safe[is.na(m_vec_safe)] = 0L
			dead_r = private$dead[m_vec_safe == 0]
			X_r    = as.matrix(KKstats$X_reservoir)
			# Candidate reduced design matrices for the reservoir
			candidates = list()
			if (ncol(X_r) == 0L){
				Xcand = matrix(w_r, ncol = 1)
				colnames(Xcand) = "w"
				candidates = list(Xcand)
			} else {
				cov_candidates = private$filtered_covariate_candidates(X_r)
				for (X in cov_candidates){
					M = matrix(w_r, ncol = 1)
					colnames(M) = "w"
					if (ncol(X) > 0){
						M = cbind(M, X)
					}
					qr_res = qr(M)
					if (qr_res$rank < ncol(M)){
						keep = qr_res$pivot[seq_len(qr_res$rank)]
						if (!(1L %in% keep)) keep = c(1L, keep)
						keep = sort(unique(keep))
						M = M[, keep, drop = FALSE]
					}
					candidates[[length(candidates) + 1L]] = M
				}
			}
			for (Xcand in candidates){
				n_params = ncol(Xcand) + 1L
				fit = .fit_standard_weibull_aft_from_matrix(
					y = y_r,
					dead = dead_r,
					X = Xcand,
					estimate_only = estimate_only,
					starts = if (!is.null(private$get_fit_warm_start_for_length("params", n_params)))
					             list(private$get_fit_warm_start_for_length("params", n_params)) else NULL,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
				)
				if (!is.null(fit) && is.finite(fit$beta) && (estimate_only || (is.finite(fit$ssq) && fit$ssq > 0))){
					private$cached_values$beta_T_reservoir = fit$beta
					private$cached_values$ssq_beta_T_reservoir = fit$ssq
					private$best_X_colnames_reservoir = colnames(Xcand)
					X_r_int = cbind("(Intercept)" = 1, as.matrix(Xcand))
					res_log_s = tryCatch(
						fast_weibull_regression_cpp(
							y    = as.numeric(y_r),
							dead = as.numeric(dead_r),
							X    = X_r_int,
							estimate_only = TRUE
						),
						error = function(e) NULL
					)
					if (!is.null(res_log_s) && isTRUE(res_log_s$converged))
						private$cached_vc_params_reservoir = as.numeric(res_log_s$log_sigma)
					return(invisible(NULL))
				}
			}
		}
	)))
)
#' Clayton Copula Combined-Likelihood Inference for KK Designs
#'
#' @export
InferenceSurvivalKKClaytonCopulaOneLik = R6::R6Class("InferenceSurvivalKKClaytonCopulaOneLik",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
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
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
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
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. Default is 0.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
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
		optimization_alg = "lbfgs",
		best_par = NULL,
		best_X_colnames = NULL,
		cached_mod = NULL,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			private$compute_basic_match_data()
			# Fixed-VC fast path: private$best_par and private$best_X_colnames survive duplicate()
			if (!is.null(private$best_par) && length(private$best_par) >= 3L &&
			    !is.null(private$best_X_colnames) && all(is.finite(tail(private$best_par, 2L)))) {
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				X_data = private$get_X()
				cov_cols = setdiff(private$best_X_colnames, "w")
				X_cand = matrix(private$w, ncol = 1L, dimnames = list(NULL, "w"))
				if (length(cov_cols) > 0L)
					X_cand = cbind(X_cand, X_data[, intersect(cov_cols, colnames(X_data)), drop = FALSE])
				X_full = cbind("(Intercept)" = 1, X_cand)
				p_full = ncol(X_full)
				vc_vals = as.numeric(tail(private$best_par, 2L))
				pair_idx_r = .complete_pair_index_matrix(m_vec)
				singleton_rows_r = setdiff(seq_len(private$n), sort(unique(as.vector(pair_idx_r)))) - 1L
				pair_idx_0 = if (nrow(pair_idx_r) > 0L) pair_idx_r - 1L else matrix(0L, 0L, 2L)
				warm = private$get_fit_warm_start_for_length("params", p_full + 2L)
				if (is.null(warm) || length(warm) != p_full + 2L)
					warm = c(rep(0, p_full), vc_vals)
				fit_fast = tryCatch(
					fast_clayton_weibull_aft_optim_cpp(
						X               = as.matrix(X_full),
						y               = as.numeric(private$y),
						dead            = as.numeric(private$dead),
						pair_idx        = pair_idx_0,
						singleton_rows  = as.integer(singleton_rows_r),
						warm_start_params = warm,
						estimate_only   = TRUE,
						optimization_alg = private$optimization_alg,
						fixed_idx       = as.integer(c(p_full + 1L, p_full + 2L)),
						fixed_values    = vc_vals
					),
					error = function(e) NULL
				)
				if (!is.null(fit_fast) && isTRUE(fit_fast$converged)) {
					b_w = as.numeric(fit_fast$params[2L])
					if (is.finite(b_w) && abs(b_w) <= private$max_abs_reasonable_coef)
						return(b_w)
				}
			}
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		get_standard_error = function(){
			private$shared()
			as.numeric(private$cached_values$s_beta_hat_T)
		},
		get_degrees_of_freedom = function() Inf,
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X = ctx$X
			y = as.numeric(ctx$y)
			dead = as.numeric(ctx$dead)
			pair_idx = ctx$pair_idx
			singleton_rows = ctx$singleton_rows
			j_treat = as.integer(ctx$j_treat %||% 1L)
			start_len = length(ctx$start %||% numeric(0))
			list(
				X = X,
				y = y,
				dead = dead,
				pair_idx = pair_idx,
				singleton_rows = singleton_rows,
				j = j_treat,
				full_fit = private$cached_mod$best_fit %||% private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm_start_params = start %||% private$get_fit_warm_start_for_length("params", start_len)
					if (length(warm_start_params) == 0L) warm_start_params = ctx$start
					warm_fisher = private$get_fit_warm_start_fisher(start_len)
					fast_clayton_weibull_aft_optim_cpp(
						X = X,
						y = y,
						dead = dead,
						pair_idx = pair_idx,
						singleton_rows = singleton_rows,
						warm_start_params = warm_start_params,
						warm_start_fisher_info = warm_fisher,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg,
						fixed_idx = j_treat,
						fixed_values = delta
					)
				},
				extract_start = function(fit){
					as.numeric(fit$params %||% fit$par %||% fit$best_par)
				},
				score = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.numeric(fit$score %||% get_clayton_weibull_aft_score_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				observed_information = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.matrix(fit$observed_information %||% fit$information %||% -get_clayton_weibull_aft_hessian_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				fisher_information = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.matrix(fit$information %||% fit$observed_information %||% -get_clayton_weibull_aft_hessian_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				information = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.matrix(fit$information %||% fit$observed_information %||% -get_clayton_weibull_aft_hessian_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll %||% fit$value)
				}
			)
		},
		filtered_covariate_candidates = function(){
			X = as.matrix(private$X)
			if (ncol(X) == 0L) return(list(matrix(nrow = private$n, ncol = 0L)))
			
			# Ensure no linearly dependent columns first
			res = drop_linearly_dependent_cols(X)
			X_reduced = res$M
			if (ncol(X_reduced) == 0L) return(list(matrix(nrow = private$n, ncol = 0L)))
			
			# Generate candidates by dropping highly correlated columns at various thresholds
			thresholds = c(0.99, 0.9, 0.7)
			candidates = list(X_reduced)
			for (thresh in thresholds) {
				X_try = drop_highly_correlated_cols(X_reduced, threshold = thresh)$M
				# Simple way to avoid duplicates in the list
				is_new = TRUE
				for (c in candidates) {
					if (ncol(c) == ncol(X_try) && all(colnames(c) == colnames(X_try))) {
						is_new = FALSE
						break
					}
				}
				if (is_new) candidates[[length(candidates) + 1L]] = X_try
			}
			candidates
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			# Optimization candidates
			w_only = matrix(private$w, ncol = 1); colnames(w_only) = "w"
			if (ncol(as.matrix(private$X)) == 0L){
				candidates = list(w_only)
			} else {
				cov_candidates = private$filtered_covariate_candidates()
				candidates = list()
				for (X_cov in cov_candidates){
					if (ncol(X_cov) == 0L){
						M = w_only
					} else {
						M = matrix(private$w, ncol = 1)
						colnames(M) = "w"
						M = cbind(M, X_cov)
						qr_M = qr(M)
						if (qr_M$rank < ncol(M)){
							keep = qr_M$pivot[seq_len(qr_M$rank)]
							if (!(1L %in% keep)) keep = c(1L, keep)
							keep = sort(unique(keep))
							M = M[, keep, drop = FALSE]
						}
						colnames(M)[1] = "w"
					}
					candidates[[length(candidates) + 1L]] = M
				}
				candidates[[length(candidates) + 1L]] = w_only
			}
			for (X in candidates){
				n_params = ncol(X) + 2L
				fit = .fit_clayton_weibull_aft(
					y = private$y,
					dead = private$dead,
					X = X,
					pair_id = m_vec,
					include_singletons = TRUE,
					estimate_only = estimate_only,
					starts = if (!is.null(private$get_fit_warm_start_for_length("params", n_params))) 
					             list(private$get_fit_warm_start_for_length("params", n_params)) else NULL,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
				)
				if (!is.null(fit) && is.finite(fit$beta) && (isTRUE(estimate_only) || (is.finite(fit$ssq) && fit$ssq > 0))){
					private$cached_values$beta_hat_T = fit$beta
					private$cached_values$s_beta_hat_T = if (is.finite(fit$ssq) && fit$ssq > 0) sqrt(fit$ssq) else NA_real_
					private$cached_values$theta = fit$theta
					private$best_par = fit$best_par
					private$best_X_colnames = colnames(X)
					private$cached_mod = fit
					private$set_fit_warm_start(fit$best_par, "params", fisher = fit$fisher_information)
					pair_idx = .complete_pair_index_matrix(m_vec) - 1L
					singleton_rows = if (length(pair_idx) > 0L) setdiff(seq_len(length(private$y)), sort(unique(as.vector(pair_idx + 1L)))) - 1L else seq_len(length(private$y)) - 1L
					private$cached_values$likelihood_test_context = list(
						X = X,
						y = as.numeric(private$y),
						dead = as.numeric(private$dead),
						pair_idx = pair_idx,
						singleton_rows = singleton_rows,
						start = fit$best_par,
						j_treat = 1L
					)
					return(invisible(NULL))
				}
			}
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			par = as.numeric(null_fit$par %||% null_fit$params %||% null_fit$b)
			if (!all(is.finite(par))) return(NULL)
			p = ncol(spec$X)
			if (length(par) < p + 2L) return(NULL)
			b_null = par[seq_len(p)]
			log_sigma = par[p + 1L]
			log_theta = par[p + 2L]
			sigma = exp(log_sigma)
			theta = exp(log_theta)
			if (!is.finite(sigma) || sigma <= 0 || !is.finite(theta) || theta <= 0) return(NULL)
			X = spec$X
			y_obs = as.numeric(spec$y)
			dead_obs = as.numeric(spec$dead)
			pair_idx = spec$pair_idx      # 0-indexed pairs matrix (n_pairs x 2)
			singleton_rows = spec$singleton_rows  # 0-indexed singleton row indices
			n = nrow(X)
			mu = as.numeric(X %*% b_null)
			T_sim = rep(NA_real_, n)
			# Simulate paired subjects via Clayton survival copula
			if (length(pair_idx) > 0L){
				pair_mat = as.matrix(pair_idx) + 1L  # convert to 1-indexed
				for (k in seq_len(nrow(pair_mat))){
					i1 = pair_mat[k, 1L]; i2 = pair_mat[k, 2L]
					U1 = runif(1L)
					W  = runif(1L)
					inner = (W * U1^(theta + 1L))^(-theta / (theta + 1L)) - U1^(-theta) + 1
					inner = max(inner, .Machine$double.xmin)
					U2 = inner^(-1 / theta)
					U2 = min(max(U2, .Machine$double.xmin), 1 - .Machine$double.xmin)
					T_sim[i1] = exp(mu[i1]) * (-log(max(U1, .Machine$double.xmin)))^sigma
					T_sim[i2] = exp(mu[i2]) * (-log(max(U2, .Machine$double.xmin)))^sigma
				}
			}
			# Simulate singleton subjects via marginal Weibull
			if (length(singleton_rows) > 0L){
				sg_idx = singleton_rows + 1L  # convert to 1-indexed
				for (i in sg_idx){
					U = runif(1L)
					T_sim[i] = exp(mu[i]) * (-log(max(U, .Machine$double.xmin)))^sigma
				}
			}
			if (!all(is.finite(T_sim)) || any(T_sim <= 0)) return(NULL)
			C_i = ifelse(dead_obs == 0, y_obs, Inf)
			y_sim = pmin(T_sim, C_i)
			dead_sim = as.numeric(T_sim <= C_i)
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			j = spec$j
			warm_full = c(as.numeric(par))
			full_res = tryCatch(
				fast_clayton_weibull_aft_optim_cpp(
					X = X, y = y_sim, dead = dead_sim,
					pair_idx = pair_idx, singleton_rows = singleton_rows,
					warm_start_params = warm_full,
					estimate_only = FALSE,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged)) return(NULL)
			full_par = as.numeric(full_res$par %||% full_res$params %||% full_res$b)
			if (length(full_par) < j || !is.finite(full_par[j])) return(NULL)
			list(
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					warm = start %||% full_par
					tryCatch(
						fast_clayton_weibull_aft_optim_cpp(
							X = X, y = y_sim, dead = dead_sim,
							pair_idx = pair_idx, singleton_rows = singleton_rows,
							warm_start_params = warm,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j, fixed_values = d
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik %||% fit$neg_ll %||% fit$value)
			)
		}
	)))
)
