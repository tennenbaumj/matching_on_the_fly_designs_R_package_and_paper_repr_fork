#' Conditional Logistic Combined-Likelihood Inference for KK Designs with Binary Responses
#'
#' Fits a single joint likelihood over all KK design data for incidence responses.
#' The matched-pair component uses the conditional logistic likelihood, and the
#' reservoir component uses the standard Bernoulli log-likelihood.
#'
#' @export
InferenceIncidKKCondLogitOneLik = R6::R6Class("InferenceIncidKKCondLogitOneLik",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart optimizer start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the combined conditional-logistic estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Numeric vector. Row weights for bootstrap.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
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
			private$cached_values$beta_hat_T = private$compute_weighted_combined_estimate(row_weights)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B Integer. Number of bootstrap samples (default 501).
		#' @param show_progress Logical. Whether to show a progress bar.
		#' @param debug Logical. Whether to return diagnostics.
		#' @param bootstrap_type Character. Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		cached_mod = NULL,
		max_abs_reasonable_coef = 1e4,
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$cached_values$likelihood_test_context = NULL
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(invisible(NULL))

			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			p             = ncol(as.matrix(private$X))
			has_reservoir = nRT > 0 && nRC > 0

			X_comb   = NULL
			y_comb   = NULL
			j_beta_T = 2L

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			if (m > 0){
				i_matched = which(m_vec > 0L)
				y_m      = private$y[i_matched]
				w_m      = private$w[i_matched]
				strata_m = m_vec[i_matched]
				X_mat    = if (p > 0L) as.matrix(private$get_X()[i_matched, , drop = FALSE]) else matrix(nrow = length(y_m), ncol = 0L)

				if (has_reservoir){
					y_r    = KKstats$y_reservoir
					w_r    = KKstats$w_reservoir
					X_r    = if (p > 0L) as.matrix(KKstats$X_reservoir) else matrix(nrow = length(y_r), ncol = 0L)
					design = build_matching_combined_clogit_design_cpp(
						as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m),
						as.double(y_r), as.double(w_r), X_r
					)
					X_comb   = design$X_comb
					y_comb   = design$y_comb
					j_beta_T = 2L
				} else {
					res = collect_discordant_pairs_cpp(
						as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m)
					)
					if (res$nd > 0){
						X_comb   = if (p > 0L) cbind(res$t_diffs, res$X_diffs) else matrix(res$t_diffs, ncol = 1L)
						y_comb   = res$y_01
						j_beta_T = 1L
					}
				}
			} else if (has_reservoir){
				y_r    = KKstats$y_reservoir
				w_r    = KKstats$w_reservoir
				X_comb = if (p > 0L) cbind(1, w_r, as.matrix(KKstats$X_reservoir)) else cbind(1, w_r)
				y_comb = y_r
			}

			if (is.null(X_comb)){
				private$cache_nonestimable_estimate("kk_clogit_combined_no_informative_data")
				return(invisible(NULL))
			}

			colnames(X_comb) = paste0("x", seq_len(ncol(X_comb)))
			colnames(X_comb)[j_beta_T] = "beta_T"

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_comb,
				required_cols = j_beta_T,
				fit_fun = function(X_fit){
					j_beta_fit = match("beta_T", colnames(X_fit))
					tryCatch(
						if (estimate_only) {
							fast_logistic_regression(
								X = X_fit,
								y = y_comb,
								optimization_alg = private$optimization_alg,
								warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
								warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
							)
						} else {
							fast_logistic_regression_with_var(
								X = X_fit,
								y = y_comb,
								j = j_beta_fit,
								optimization_alg = private$optimization_alg,
								warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
								warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
							)
						},
						error = function(e) NULL
					)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod)) return(FALSE)
					j_beta_fit = match("beta_T", colnames(X_fit))
					if (!is.finite(j_beta_fit) || is.na(j_beta_fit)) return(FALSE)
					beta = suppressWarnings(as.numeric(mod$b[j_beta_fit]))
					if (!is.finite(beta) || abs(beta) > private$max_abs_reasonable_coef) return(FALSE)
					if (estimate_only) return(TRUE)
					ssq = suppressWarnings(as.numeric(mod$ssq_b_j))
					is.finite(ssq) && ssq > 0 && sqrt(ssq) <= private$max_abs_reasonable_coef
				}
			)
			mod = attempt$fit
			j_beta_T = if (!is.null(attempt$X)) match("beta_T", colnames(attempt$X)) else NA_integer_
			if (is.null(mod) || !is.finite(j_beta_T) || is.na(j_beta_T) || !is.finite(mod$b[j_beta_T])){
				private$cache_nonestimable_estimate("kk_clogit_combined_fit_failed")
				return(invisible(NULL))
			}
			if (max(abs(mod$b), na.rm = TRUE) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_clogit_combined_extreme_coefficients")
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(mod$b[j_beta_T])
			private$cached_mod = mod
			private$set_fit_warm_start(as.numeric(mod$b), "beta")
			private$cached_values$likelihood_test_context = list(
				X = attempt$X,
				y = y_comb,
				j_treat = j_beta_T
			)
			if (!estimate_only) {
				se = sqrt(mod$ssq_b_j)
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se <= private$max_abs_reasonable_coef) se else NA_real_
				if (!is.finite(private$cached_values$s_beta_hat_T)){
					private$cache_nonestimable_se("kk_clogit_combined_standard_error_unavailable")
					return(invisible(NULL))
				}
			}
			private$clear_nonestimable_state()
			private$cached_values$df = NA_real_
			invisible(NULL)
		},
		get_standard_error = function(){
			private$shared_combined_likelihood(estimate_only = FALSE)
			se = private$cached_values$s_beta_hat_T
			if (is.null(se) || length(se) == 0L) {
				return(NA_real_)
			}
			as.numeric(se)[1L]
		},
		supports_likelihood_tests = function() TRUE,
		get_likelihood_test_spec = function(){
			private$shared_combined_likelihood(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta){
					fast_logistic_regression_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						fixed_idx = j_treat,
						fixed_values = delta
					)
				},
				score = function(fit){
					get_logistic_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_logistic_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_logistic_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				information = function(fit){
					-get_logistic_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					log_denom = ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta)))
					-sum(y * eta - log_denom)
				}
			)
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		compute_weighted_combined_estimate = function(row_weights){
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(NA_real_)
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			p             = ncol(as.matrix(private$X))
			has_reservoir = nRT > 0 && nRC > 0
			kk_w = kk_pair_and_reservoir_bootstrap_weights(private, row_weights)
			X_comb = NULL
			y_comb = NULL
			w_comb = NULL
			j_beta_T = 2L
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			if (m > 0){
				i_matched = which(m_vec > 0L)
				y_m = private$y[i_matched]
				w_m = private$w[i_matched]
				strata_m = m_vec[i_matched]
				X_mat = if (p > 0L) as.matrix(private$get_X()[i_matched, , drop = FALSE]) else matrix(nrow = length(y_m), ncol = 0L)
				if (has_reservoir){
					y_r = KKstats$y_reservoir
					w_r = KKstats$w_reservoir
					X_r = if (p > 0L) as.matrix(KKstats$X_reservoir) else matrix(nrow = length(y_r), ncol = 0L)
					design = build_matching_combined_clogit_design_cpp(
						as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m),
						as.double(y_r), as.double(w_r), X_r
					)
					X_comb = design$X_comb
					y_comb = design$y_comb
					w_comb = c(kk_w$pair_weights, kk_w$reservoir_weights)
					j_beta_T = 2L
				} else {
					res = collect_discordant_pairs_cpp(as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m))
					if (res$nd > 0){
						X_comb = if (p > 0L) cbind(res$t_diffs, res$X_diffs) else matrix(res$t_diffs, ncol = 1L)
						y_comb = res$y_01
						w_comb = kk_w$pair_weights[seq_len(res$nd)]
						j_beta_T = 1L
					}
				}
			} else if (has_reservoir){
				y_r = KKstats$y_reservoir
				w_r = KKstats$w_reservoir
				X_comb = if (p > 0L) cbind(1, w_r, as.matrix(KKstats$X_reservoir)) else cbind(1, w_r)
				y_comb = y_r
				w_comb = kk_w$reservoir_weights
			}
			if (is.null(X_comb) || is.null(w_comb)) return(NA_real_)
			colnames(X_comb) = paste0("x", seq_len(ncol(X_comb)))
			colnames(X_comb)[j_beta_T] = "beta_T"
			ok = is.finite(w_comb) & w_comb > 0 & is.finite(y_comb)
			if (!any(ok)) return(NA_real_)
			X_comb = X_comb[ok, , drop = FALSE]
			y_comb = y_comb[ok]
			w_comb = w_comb[ok]
			mod = tryCatch(
				fast_logistic_regression_weighted_cpp(
					X = X_comb,
					y = y_comb,
					weights = w_comb,
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_comb)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_comb)),
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			j_beta_fit = match("beta_T", colnames(X_comb))
			if (is.null(mod) || is.na(j_beta_fit) || length(mod$b) < j_beta_fit || !is.finite(mod$b[j_beta_fit])) return(NA_real_)
			as.numeric(mod$b[j_beta_fit])
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			X = spec$X
			j = spec$j
			n = nrow(X)
			b_null = as.numeric(null_fit$b)
			pi = plogis(as.numeric(X %*% b_null))
			y_sim = as.integer(rbinom(n, 1L, pi))
			full_res = tryCatch(
				fast_logistic_regression_cpp(
					X = X, y = y_sim,
					optimization_alg = private$optimization_alg %||% "lbfgs"
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !is.finite(full_res$b[j])) return(NULL)
			list(
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					tryCatch(
						fast_logistic_regression_with_var_cpp(
							X = X, y = y_sim, j = j,
							warm_start_beta = start %||% as.numeric(full_res$b),
							fixed_idx = j, fixed_values = d
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit){
					eta = as.numeric(X %*% as.numeric(fit$b))
					log_denom = ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta)))
					-sum(y_sim * eta - log_denom)
				}
			)
		}
	)))
)
#' Conditional Logistic IVWC Inference for KK Designs with Binary Responses
#'
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferenceIncidKKCondLogitIVWC = R6::R6Class("InferenceIncidKKCondLogitIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
		},
			#' @description Compute the treatment effect estimate.
			#' @param estimate_only Logical. If TRUE, skip variance component calculations.
			compute_estimate = function(estimate_only = FALSE){
				private$shared(estimate_only = estimate_only)
				private$cached_values$beta_hat_T
			},
			#' @description Computes an approximate confidence interval.
			#' @param alpha Numeric. Significance level (default 0.05).
			compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B Integer. Number of bootstrap samples (default 501).
		#' @param show_progress Logical. Whether to show a progress bar.
			#' @param debug Logical. Whether to return diagnostics.
			#' @param bootstrap_type Character. Optional resampling scheme.
			#' @return A numeric vector of bootstrap estimates.
			approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
				eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
			}
		))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			private$compute_basic_match_data()
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(invisible(NULL))
			X_covars = private$X
			
			# --- Matched pairs: Conditional Logistic ---
			if (KKstats$m > 0){
				private$clogit_for_matched_pairs(KKstats, X_covars)
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0
			# --- Reservoir: Logistic Regression ---
			if (KKstats$nRT > 0 && KKstats$nRC > 0){
				private$logistic_for_reservoir(KKstats, X_covars)
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0
			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = Inf
		},
		clogit_for_matched_pairs = function(KKstats, X_covars){
			yTs = KKstats$yTs_matched
			yCs = KKstats$yCs_matched
			y_m_all = yTs - yCs # 1 if (1,0), -1 if (0,1)
			i_m_disc = which(abs(y_m_all) == 1)
			if (length(i_m_disc) == 0L) {
				private$cached_values$beta_T_matched = NA_real_
				private$cached_values$ssq_beta_T_matched = NA_real_
				return(invisible(NULL))
			}
			# X_m should be the difference X_T - X_C. 
			# KKstats$X_matched_diffs is already X_T - X_C.
			# For treatment effect, the difference is 1 - 0 = 1.
			X_m = cbind(treatment = 1, KKstats$X_matched_diffs[i_m_disc, , drop = FALSE])
			y_m = (y_m_all[i_m_disc] + 1) / 2 # 1 if (1,0), 0 if (0,1)
			
			# Conditional logistic fit via internal Rcpp
			fit = tryCatch(fast_logistic_regression_with_var_cpp(X_m, y_m, j = 1L), error = function(e) NULL)
			if (is.null(fit) || !isTRUE(fit$converged)){
				private$cached_values$beta_T_matched = NA_real_
				private$cached_values$ssq_beta_T_matched = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_T_matched = as.numeric(fit$b[1])
			private$cached_values$ssq_beta_T_matched = as.numeric(fit$ssq_b_j)
		},
		logistic_for_reservoir = function(KKstats, X_covars){
			y_r = KKstats$y_reservoir
			w_r = KKstats$w_reservoir
			if (is.null(X_covars) || ncol(X_covars) == 0){
				X_r = cbind(`(Intercept)` = 1, treatment = w_r)
			} else {
				X_r = cbind(`(Intercept)` = 1, treatment = w_r, KKstats$X_reservoir)
			}
			# Logistic fit via internal Rcpp
			fit = tryCatch(fast_logistic_regression_with_var_cpp(X_r, y_r, j = 2L), error = function(e) NULL)
			if (is.null(fit) || !isTRUE(fit$converged)){
				private$cached_values$beta_T_reservoir = NA_real_
				private$cached_values$ssq_beta_T_reservoir = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_T_reservoir = as.numeric(fit$b[2])
			private$cached_values$ssq_beta_T_reservoir = as.numeric(fit$ssq_b_j)
		}
	)))
)
