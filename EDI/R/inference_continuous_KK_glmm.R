#' Linear Mixed Model Inference for KK Designs with Continuous Response
#'
#' Fits a linear mixed model for continuous responses under a KK
#' matching-on-the-fly design. The matched-pair strata enter as a subject-level
#' random intercept \code{(1 | group_id)}, accounting for within-pair correlation.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an
#' internal Rcpp/L-BFGS routine that requires no external packages. Set
#' \code{use_rcpp = FALSE} to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinKKGLMM = R6::R6Class("InferenceContinKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = utils::modifyList(as.list(InferenceMixinKKGLMMShared$public), list(
		#' @description Initialize a KK GLMM inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the optimised Rcpp
		#'   Gaussian LMM implementation (no external package required). If \code{FALSE},
		#'   use \pkg{glmmTMB}.
		#' @param use_gls_fast_path Logical. If \code{TRUE} (default), use a fast GLS
		#'   estimator (no optimisation) for \code{estimate_only} calls during
		#'   randomisation inference once variance components are cached from a prior
		#'   full fit.  Statistically exact: fixing VC at the null-fit MLE and
		#'   permuting only the treatment assignment gives a valid permutation test by
		#'   exchangeability.  Set \code{FALSE} to always run full L-BFGS.
		#' @param use_gls_fast_path_bootstrap Logical. If \code{TRUE}, also use the
		#'   fast GLS estimator for non-studentised bootstrap draws
		#'   (\code{estimate_only = TRUE} weighted calls).  Asymptotically valid by
		#'   the plug-in principle (VC orthogonal to beta_T in the Fisher information),
		#'   but not exact in finite samples.  Default \code{FALSE}.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		#' @param optimization_alg The optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, use_gls_fast_path = TRUE, use_gls_fast_path_bootstrap = FALSE, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
				assertFlag(use_gls_fast_path)
				assertFlag(use_gls_fast_path_bootstrap)
			}
			# If using Rcpp, skip glmmTMB package check in the parent initialize.
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_glmm_shared(des_obj)
			private$use_rcpp = use_rcpp
			private$use_gls_fast_path = use_gls_fast_path
			private$use_gls_fast_path_bootstrap = use_gls_fast_path_bootstrap
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
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
					private$cached_values$df = Inf
					private$cached_values$summary_table = NULL
					return(private$cached_values$beta_hat_T)
				}
			}
			if (isTRUE(private$use_rcpp)) {
				rcpp_result = private$weighted_rcpp_estimate(row_weights, estimate_only = estimate_only)
				if (!is.null(rcpp_result)) {
					if (is.list(rcpp_result)) {
						private$cached_values$beta_hat_T = rcpp_result$beta
						private$cached_values$s_beta_hat_T = rcpp_result$se
					} else {
						private$cached_values$beta_hat_T = as.numeric(rcpp_result)[1L]
						private$cached_values$s_beta_hat_T = NA_real_
					}
					private$cached_values$df = Inf
					private$cached_values$summary_table = NULL
					return(private$cached_values$beta_hat_T)
				}
			}
			result = private$compute_weighted_glmm_bootstrap_estimate(row_weights, estimate_only = estimate_only)
			if (is.list(result)) {
				private$cached_values$beta_hat_T = result$beta
				private$cached_values$s_beta_hat_T = result$se
			} else {
				private$cached_values$beta_hat_T = as.numeric(result)[1L]
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = Inf
			private$cached_values$summary_table = NULL
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKGLMMShared$private), list(
		use_rcpp = TRUE,
		use_gls_fast_path = TRUE,
		use_gls_fast_path_bootstrap = FALSE,
		cached_vc_params = NULL,   # c(log_sigma_e, log_sigma_b) from last full Rcpp fit
		glmm_response_type = function() "continuous",
		glmm_family        = function() stats::gaussian(link = "identity"),
		get_complexity_tier = function() "heavy",
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		# ── Dispatch ─────────────────────────────────────────────────────────
		shared = function(estimate_only = FALSE){
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				private$shared_glmm_tmb(estimate_only)
			}
		},
		# ── Rcpp Gaussian LMM path ────────────────────────────────────────────
		shared_rcpp = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			private$clear_nonestimable_state()
			private$cached_mod = NULL
			private$cached_values$likelihood_test_context = NULL
			# ── Group IDs (same logic as fit_glmm_on_data) ───────────────────
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			# ── Design matrix: (intercept, w, covariates) ────────────────────
			X_fit = private$create_design_matrix()
			# create_design_matrix uses "treatment"; Rcpp path expects "w" in the search for j_T_r
			if ("treatment" %in% colnames(X_fit))
				colnames(X_fit)[colnames(X_fit) == "treatment"] = "w"
			X_fit = as.matrix(X_fit)
			# Treatment column index (1-based in R, 0-based in C++ is handled via ssq_b_T)
			j_T_r = which(colnames(X_fit) == "w")
			if (length(j_T_r) == 0L) j_T_r = 2L   # fallback: second column
			
			start_len = ncol(X_fit) + 2L # beta + log_sigma + log_tau

			# GLS fast path: skip L-BFGS for estimate_only calls when VC is cached.
			# Exact for permutation tests: fixing VC at null-fit MLE and permuting w
			# gives a valid test statistic by exchangeability. Skip when sigma_b≈0:
			# L-BFGS terminates in 0 iters (already at boundary) so GLS overhead
			# exceeds the saving — about 2x slower than the full REML path.
			if (estimate_only && isTRUE(private$use_gls_fast_path) && !is.null(private$cached_vc_params) &&
					private$cached_vc_params[2L] > -8) {
				gls_b = tryCatch(
					EDI:::fast_gaussian_lmm_gls_cpp(
						X            = X_fit,
						y            = as.numeric(private$y),
						group_id     = as.integer(group_id),
						log_sigma_e  = private$cached_vc_params[1L],
						log_sigma_b  = private$cached_vc_params[2L]
					),
					error = function(e) NULL
				)
				if (!is.null(gls_b) && all(is.finite(gls_b))) {
					beta_hat_T_gls = as.numeric(gls_b[j_T_r])
					if (is.finite(beta_hat_T_gls) && abs(beta_hat_T_gls) <= private$max_abs_reasonable_coef) {
						private$cached_values$beta_hat_T = beta_hat_T_gls
						private$cached_values$df = Inf
						return(invisible(NULL))
					}
				}
			}

			ws_args = private$get_optimal_warm_start_config(start_len)

			fit = tryCatch(
				fast_gaussian_lmm_cpp(
					X             = X_fit,
					y             = as.numeric(private$y),
					group_id      = as.integer(group_id),
					warm_start_params = ws_args$start_params,
					estimate_only = estimate_only,
					maxit         = 300L,
					eps_g         = 1e-6,
					optimization_alg = private$optimization_alg,
					warm_start_fisher_info = ws_args$warm_start_fisher_info
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				return(private$shared_glmm_tmb(estimate_only = estimate_only))
			}
			beta_hat_T = as.numeric(fit$b[j_T_r])
			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				return(private$shared_glmm_tmb(estimate_only = estimate_only))
			}
			private$cached_mod = fit
			p_ncol = ncol(X_fit)
			private$cached_vc_params = as.numeric(fit$b[c(p_ncol + 1L, p_ncol + 2L)])
			full_params = as.numeric(c(fit$b, fit$log_sigma, fit$log_tau))
			private$set_fit_warm_start(full_params, "params", fisher = fit$fisher_information)
			private$cached_values$likelihood_test_context = list(
				X = X_fit,
				y = as.numeric(private$y),
				group_id = as.integer(group_id),
				j_treat = j_T_r,
				start = as.numeric(fit$b)
			)
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		}
		,
		get_likelihood_test_spec = function(){
			if (!isTRUE(private$use_rcpp)) return(NULL)
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			group_id = as.integer(ctx$group_id)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				group_id = group_id,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm_start_params = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start
					fast_gaussian_lmm_cpp(
						X = X_fit,
						y = y,
						group_id = group_id,
						warm_start_params = warm_start_params,
						estimate_only = FALSE,
						maxit = 300L,
						eps_g = 1e-6,
						fixed_idx = j_treat,
						fixed_values = delta,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					as.numeric(get_gaussian_lmm_score_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				observed_information = function(fit){
					as.matrix(get_gaussian_lmm_fisher_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				fisher_information = function(fit){
					as.matrix(get_gaussian_lmm_fisher_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				information = function(fit){
					as.matrix(get_gaussian_lmm_fisher_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		},
		# Fast weighted LMM fit for bootstrap iterations — avoids glmmTMB entirely.
		weighted_rcpp_estimate = function(row_weights, estimate_only = TRUE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			X_fit = private$create_design_matrix()
			if ("treatment" %in% colnames(X_fit))
				colnames(X_fit)[colnames(X_fit) == "treatment"] = "w"
			X_fit = as.matrix(X_fit)
			j_T_r = which(colnames(X_fit) == "w")
			if (length(j_T_r) == 0L) j_T_r = 2L

			# GLS fast path for non-studentised bootstrap (opt-in, default FALSE).
			# Asymptotically valid via plug-in principle; not exact in finite samples.
			if (estimate_only && isTRUE(private$use_gls_fast_path_bootstrap) && !is.null(private$cached_vc_params)) {
				gls_b = tryCatch(
					EDI:::fast_gaussian_lmm_gls_cpp(
						X            = X_fit,
						y            = as.numeric(private$y),
						group_id     = as.integer(group_id),
						log_sigma_e  = private$cached_vc_params[1L],
						log_sigma_b  = private$cached_vc_params[2L],
						weights      = as.numeric(row_weights)
					),
					error = function(e) NULL
				)
				if (!is.null(gls_b) && all(is.finite(gls_b))) {
					beta_hat_T_gls = as.numeric(gls_b[j_T_r])
					if (is.finite(beta_hat_T_gls) && abs(beta_hat_T_gls) <= private$max_abs_reasonable_coef) {
						return(beta_hat_T_gls)
					}
				}
			}

			start_len = ncol(X_fit) + 2L
			ws_args = private$get_optimal_warm_start_config(start_len)
			fit = tryCatch(
				EDI:::fast_gaussian_lmm_cpp(
					X             = X_fit,
					y             = as.numeric(private$y),
					group_id      = as.integer(group_id),
					warm_start_params = ws_args$start_params,
					estimate_only = estimate_only,
					maxit         = 300L,
					eps_g         = 1e-6,
					optimization_alg = private$optimization_alg,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					weights       = as.numeric(row_weights)
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
			beta_hat_T = as.numeric(fit$b[j_T_r])
			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) return(NULL)
			if (estimate_only) return(beta_hat_T)
			vcov = fit$vcov
			se = if (!is.null(vcov) && is.finite(vcov[j_T_r, j_T_r]) && vcov[j_T_r, j_T_r] > 0)
				sqrt(vcov[j_T_r, j_T_r]) else NA_real_
			list(beta = beta_hat_T, se = se)
		},
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		simulate_under_lik_null = function(spec, delta, null_fit){
			# fast_gaussian_lmm_cpp bundles all params in $b: c(betas, log_sigma, log_tau)
			p      = ncol(spec$X)
			b_all  = as.numeric(null_fit$b)
			if (length(b_all) < p + 2L) return(NULL)
			b_null    = b_all[seq_len(p)]
			sigma_u   = exp(b_all[p + 1L])
			sigma_e   = exp(b_all[p + 2L])
			X = spec$X
			group_id = spec$group_id
			n = nrow(X)
			K = max(group_id)
			u_g = rnorm(K, 0, sigma_u)
			mu = as.numeric(X %*% b_null) + u_g[group_id]
			y_sim = rnorm(n, mu, sigma_e)
			j = spec$j
			full_res = tryCatch(
				fast_gaussian_lmm_cpp(
					X = X, y = y_sim, group_id = group_id,
					estimate_only = FALSE, maxit = 300L, eps_g = 1e-6,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$b[j])) return(NULL)
			list(
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					tryCatch(
						fast_gaussian_lmm_cpp(
							X = X, y = y_sim, group_id = group_id,
							warm_start_params = start %||% as.numeric(full_res$b),
							estimate_only = FALSE, maxit = 300L, eps_g = 1e-6,
							fixed_idx = j, fixed_values = d,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik %||% fit$neg_ll) }
			)
		}
	))
)
