#' GEE Inference for KK Designs with Ordinal Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{multgee})
#' for ordinal responses under a KK matching-on-the-fly design using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKGEE = R6::R6Class("InferenceOrdinalKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = utils::modifyList(InferenceMixinKKGEEShared$public, list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				if (!check_package_installed("multgee")){
					stop("Package 'multgee' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_gee_shared(des_obj, use_rcpp = FALSE)
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_gee_dispatch(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared_gee_dispatch(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared_gee_dispatch(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Recomputes the KK ordinal GEE treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (length(row_weights) > 0L && all(is.finite(row_weights)) &&
			    (max(row_weights) - min(row_weights)) <= sqrt(.Machine$double.eps)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$beta_hat_T = beta_hat_T
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$df = Inf
					private$cached_values$summary_table = NULL
					return(private$cached_values$beta_hat_T)
				}
			}
			pred_df = private$gee_predictors_df()
			ok = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
			if (!any(ok)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			X_fit = as.matrix(pred_df[ok, , drop = FALSE])
			y_fit = as.numeric(private$y[ok])
			n_params = ncol(X_fit) + length(sort(unique(y_fit))) - 1L
			mod = tryCatch(
				fast_ordinal_regression_weighted_cpp(
					X = X_fit,
					y = y_fit,
					weights = as.numeric(row_weights[ok]),
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (!is.null(mod) && !is.null(mod$params)) {
				private$set_fit_warm_start(as.numeric(mod$params), "params", fisher = mod$fisher_information)
			}
			beta_hat_T = if (is.null(mod) || length(mod$b) < 1L) NA_real_ else as.numeric(mod$b[1L])
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = Inf
			private$cached_values$summary_table = NULL
			private$cached_values$beta_hat_T
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKGEEShared$private), list(
		gee_response_type = function() "ordinal",
		gee_family        = function() stats::binomial(link = "logit"),
		# Ordinal response requires ordLORgee, not geeglm.
		shared_gee_dispatch = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			pred_df = private$gee_predictors_df()
			dat = data.frame(y = factor(private$y, ordered = TRUE), pred_df, group_id = group_id)
			dat = dat[order(dat$group_id), ]
			id_sorted = dat$group_id
			
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			formula_gee = stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))
			
			bstart = private$get_fit_warm_start_for_length("beta", length(fixed_terms) + nlevels(dat$y) - 1L)
			
			mod = tryCatch({
				utils::capture.output(m <- suppressMessages(suppressWarnings(
					multgee::ordLORgee(
						formula_gee,
						data   = dat,
						id     = id_sorted,
						LORstr = "uniform",
						link   = "logit",
						bstart = bstart
					)
				)))
				m
			}, error = function(e) NULL)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("ordinal_kk_gee_fit_unavailable")
				return(invisible(NULL))
			}
			beta = stats::coef(mod)
			private$set_fit_warm_start(beta, "beta")
			
			j_treat = private$gee_treatment_index(beta)
			private$cached_values$beta_hat_T = as.numeric(beta[j_treat])
			if (estimate_only) return(invisible(NULL))
			vcov_robust = tryCatch(stats::vcov(mod), error = function(e) NULL)
			if (is.null(vcov_robust)) {
				private$cached_values$s_beta_hat_T = NA_real_
			} else {
				private$cached_values$s_beta_hat_T = sqrt(as.numeric(vcov_robust[j_treat, j_treat]))
			}
			private$cached_values$df = Inf
			private$cached_values$summary_table = summary(mod)$coefficients
		}
	))
)
#' GLMM Inference for KK Designs with Ordinal Response
#'
#' Fits a cumulative-logit mixed model (proportional odds) for ordinal responses
#' under a KK matching-on-the-fly design. The random intercept per matched pair is
#' integrated out via Gauss-Hermite quadrature.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an internal
#' Rcpp/L-BFGS routine that requires no external packages. Set \code{use_rcpp = FALSE}
#' to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKGLMM = R6::R6Class("InferenceOrdinalKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = utils::modifyList(as.list(InferenceMixinKKGLMMShared$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use internal Rcpp.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_glmm_shared(des_obj)
			private$use_rcpp = use_rcpp
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
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
		},
		#' @description Recomputes the KK ordinal GLMM treatment estimate under
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
					private$cached_values$df = Inf
					private$cached_values$summary_table = NULL
					return(private$cached_values$beta_hat_T)
				}
			}
			pred_df = private$glmm_predictors_df()
			ok = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
			if (!any(ok)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			X_fit = as.matrix(pred_df[ok, , drop = FALSE])
			y_fit = as.numeric(private$y[ok])
			n_params = ncol(X_fit) + length(sort(unique(y_fit))) - 1L
			mod = tryCatch(
				fast_ordinal_regression_weighted_cpp(
					X = X_fit,
					y = y_fit,
					weights = as.numeric(row_weights[ok]),
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (!is.null(mod) && !is.null(mod$params)) {
				private$set_fit_warm_start(as.numeric(mod$params), "params", fisher = mod$fisher_information)
			}
			beta_hat_T = if (is.null(mod) || length(mod$b) < 1L) NA_real_ else as.numeric(mod$b[1L])
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = Inf
			private$cached_values$summary_table = NULL
			private$cached_values$beta_hat_T
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKGLMMShared$private), list(
		use_rcpp = TRUE,
		glmm_response_type  = function() "ordinal",
		glmm_family         = function() glmmTMB::cumulative(link = "logit"),
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		shared = function(estimate_only = FALSE){
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				private$shared_glmm_tmb(estimate_only)
			}
		},
		shared_rcpp = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()
			private$cached_mod = NULL
			private$cached_values$likelihood_test_context = NULL
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			# X WITHOUT intercept (cutpoints serve as intercepts)
			if (ncol(as.matrix(private$X)) > 0){
				X_fit = as.matrix(private$glmm_predictors_df())  # [w, cov1, ...]
			} else {
				X_fit = matrix(private$w, ncol = 1L, dimnames = list(NULL, "w"))
			}
			# Convert y to 1-indexed integers
			y_levels = sort(unique(private$y))
			K = length(y_levels)
			y = as.integer(match(private$y, y_levels))
			n_alpha = K - 1L
			# Treatment is always the first column of X_fit (j_T = 0, 0-based)
			j_T = 0L
			
			start_len = n_alpha + ncol(X_fit) + 1L
			warm_start = private$get_fit_warm_start_for_length("params", start_len)
			
			# Warm start from fixed-effects ordinal MLE to avoid divergence if no cache
			start = if (!is.null(warm_start)) warm_start else tryCatch({
				nore = fast_ordinal_regression_cpp(X_fit, as.numeric(y) - 1L)
				alpha_direct = as.numeric(nore$alpha)  # K-1 direct cutpoints
				beta_nore    = as.numeric(nore$b)      # p betas
				# Convert direct alphas to log-diff parameterization
				alpha_par = numeric(n_alpha)
				if (n_alpha >= 1L) alpha_par[1L] = alpha_direct[1L]
				if (n_alpha >= 2L) {
					for (k in 2L:n_alpha) {
						diff_k = alpha_direct[k] - alpha_direct[k - 1L]
						alpha_par[k] = if (diff_k > 0) log(diff_k) else 0.0
					}
				}
				c(alpha_par, beta_nore, -3.0)  # log_sigma = -3 (small random effect)
			}, error = function(e) NULL)
			
			fit = tryCatch(
				fast_ordinal_glmm_cpp(
					X          = X_fit,
					y          = y,
					group_id   = as.integer(group_id),
					K          = K,
					j_T        = j_T,
					smart_cold_start = private$smart_cold_start_default,
					estimate_only = estimate_only,
					warm_start_params = start,
					eps_g      = 1e-3,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(start_len),
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				private$cache_nonestimable_estimate("kk_glmm_rcpp_failed")
				return(invisible(NULL))
			}
			# b is the beta vector (no cutpoints); treatment is at index j_T+1 (1-based R)
			beta_hat_T = as.numeric(fit$b[j_T + 1L])
			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				private$cache_nonestimable_estimate("kk_glmm_rcpp_nonestimable")
				return(invisible(NULL))
			}
			private$cached_mod = fit
			full_params = as.numeric(c(fit$alpha, fit$b, fit$log_sigma))
			private$set_fit_warm_start(full_params, "params", fisher = fit$fisher_information)

			private$cached_values$likelihood_test_context = list(
				X = X_fit,
				y = y,
				group_id = as.integer(group_id),
				K = K,
				j_treat = length(fit$alpha) + 1L,
				n_gh = 20L,
				start = full_params
			)
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_T
			if (!is.null(ssq) && is.finite(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				j_in_full = length(fit$alpha) + 1L
				hess = tryCatch(
					get_ordinal_glmm_hessian_cpp(X_fit, y, as.integer(group_id), full_params, K, n_gh = 20L),
					error = function(e) NULL
				)
				se_fallback = NA_real_
				if (!is.null(hess) && is.matrix(hess) && nrow(hess) >= j_in_full) {
					vcov_hess = tryCatch(solve(hess), error = function(e) NULL)
					if (!is.null(vcov_hess) && is.finite(vcov_hess[j_in_full, j_in_full]) && vcov_hess[j_in_full, j_in_full] > 0) {
						se_fallback = sqrt(vcov_hess[j_in_full, j_in_full])
					}
				}
				private$cached_values$s_beta_hat_T = se_fallback
			}
		},
		get_likelihood_test_spec = function(){
			if (!isTRUE(private$use_rcpp)) return(NULL)
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.integer(ctx$y)
			group_id = as.integer(ctx$group_id)
			K = as.integer(ctx$K)
			j_treat = as.integer(ctx$j_treat)
			n_gh = as.integer(ctx$n_gh %||% 20L)
			list(
				X = X_fit,
				y = y,
				group_id = group_id,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					run_fit = function(s){
						n_params = length(ctx$start)
						tryCatch(
							fast_ordinal_glmm_cpp(
								X = X_fit,
								y = y,
								group_id = group_id,
								K = K,
								j_T = 0L,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = FALSE,
								n_gh = n_gh,
								max_abs_log_sigma = 8.0,
								maxit = 300L,
								eps_g = 1e-3,
								warm_start_params = s,
								warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
								optimization_alg = private$optimization_alg %||% "lbfgs",
								fixed_idx = j_treat,
								fixed_values = delta
							),
							error = function(e) NULL
						)
					}
					warm_start = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start
					fit = run_fit(warm_start)
					# If warm start caused failure, retry with the canonical default start
					if (is.null(fit) || !isTRUE(fit$converged)) {
						if (!identical(warm_start, ctx$start)) {
							fit2 = run_fit(ctx$start)
							if (!is.null(fit2) && isTRUE(fit2$converged)) fit = fit2
						}
					}
					if (!is.null(fit)) {
						fit$params = tryCatch(as.numeric(c(fit$alpha, fit$b, fit$log_sigma)), error = function(e) NULL)
					}
					fit
				},
				extract_start = function(fit){
					as.numeric(c(fit$alpha, fit$b, fit$log_sigma))
				},
				score = function(fit){
					as.numeric(get_ordinal_glmm_score_cpp(X_fit, y, group_id, as.numeric(c(fit$alpha, fit$b, fit$log_sigma)), K, n_gh = n_gh))
				},
				observed_information = function(fit){
					as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% get_ordinal_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(c(fit$alpha, fit$b, fit$log_sigma)), K, n_gh = n_gh))
				},
				fisher_information = function(fit){
					as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% get_ordinal_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(c(fit$alpha, fit$b, fit$log_sigma)), K, n_gh = n_gh))
				},
				information = function(fit){
					as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% get_ordinal_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(c(fit$alpha, fit$b, fit$log_sigma)), K, n_gh = n_gh))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		}
	))
)
