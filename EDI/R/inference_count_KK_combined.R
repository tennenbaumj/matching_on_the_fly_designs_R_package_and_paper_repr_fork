#' GLMM Inference for KK Designs with Count Response
#'
#' Fits a Poisson GLMM for count responses under a KK matching-on-the-fly design.
#' The random intercept per matched pair is integrated out via Gauss-Hermite quadrature.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an internal
#' Rcpp routine. Set \code{use_rcpp = FALSE} to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountKKGLMM = R6::R6Class("InferenceCountKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = utils::modifyList(as.list(InferenceMixinKKGLMMShared$public), list(
		#' @description Initialize a KK Poisson GLMM inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the internal Rcpp Poisson GLMM.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			self$set_optimization_alg(optimization_alg, allow_irls = TRUE)
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
		#' @description Recomputes the KK Poisson-GLMM treatment estimate under
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
					private$cached_values$df = Inf
					private$cached_values$summary_table = NULL
					return(private$cached_values$beta_hat_T)
				}
			}
			beta_hat_T = private$compute_weighted_glmm_bootstrap_estimate(row_weights)
			private$cached_values$beta_hat_T = as.numeric(beta_hat_T)[1L]
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = Inf
			private$cached_values$summary_table = NULL
			private$cached_values$beta_hat_T
		},
		#' @description Computes a Wald confidence interval.
		#' @param alpha Significance level.
		compute_wald_confidence_interval = function(alpha = 0.05){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes a Wald two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_wald_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes a likelihood-ratio confidence interval, conservatively
		#'   widened against the design-aware Wald interval.
		#' @param alpha Significance level.
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			ci_model = super$compute_lik_ratio_confidence_interval(alpha = alpha)
			ci_design = self$compute_wald_confidence_interval(alpha = alpha)
			.conservative_kk_onelik_ci(ci_model, ci_design, alpha = alpha)
		},
		#' @description Computes a likelihood-ratio two-sided p-value, conservatively
		#'   calibrated against the design-aware Wald p-value.
		#' @param delta Null treatment effect value.
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			p_model = super$compute_lik_ratio_two_sided_pval(delta = delta)
			p_design = self$compute_wald_two_sided_pval(delta = delta)
			.conservative_kk_onelik_pval(p_model, p_design)
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			switch(
				self$get_testing_type(),
				wald = self$compute_wald_confidence_interval(alpha = alpha),
				lik_ratio = self$compute_lik_ratio_confidence_interval(alpha = alpha),
				self$compute_wald_confidence_interval(alpha = alpha)
			)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			switch(
				self$get_testing_type(),
				wald = self$compute_wald_two_sided_pval(delta = delta),
				lik_ratio = self$compute_lik_ratio_two_sided_pval(delta = delta),
				self$compute_wald_two_sided_pval(delta = delta)
			)
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKGLMMShared$private), list(
		use_rcpp = TRUE,
		glmm_response_type = function() "count",
		compute_weighted_glmm_bootstrap_estimate = function(row_weights){
			if (!isTRUE(private$use_rcpp)) {
				return(callSuper())
			}
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			# drop rows with zero or non-finite weight
			ok = is.finite(row_weights) & row_weights > 0 & is.finite(as.numeric(private$y))
			if (!any(ok)) return(NA_real_)
			if (ncol(as.matrix(private$X)) > 0) {
				X_fit = private$create_design_matrix()
			} else {
				X_fit = cbind(`(Intercept)` = 1, w = private$w)
			}
			X_fit = as.matrix(X_fit)[ok, , drop = FALSE]
			y_ok = as.numeric(private$y)[ok]
			gid_ok = as.integer(group_id)[ok]
			rw_ok = as.numeric(row_weights)[ok]
			j_T = 1L
			n_params = ncol(X_fit) + 1L
			fit = tryCatch(
				fast_poisson_glmm_cpp(
					X        = X_fit,
					y        = y_ok,
					group_id = gid_ok,
					j_T      = j_T,
					row_weights = rw_ok,
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					smart_cold_start  = private$smart_cold_start_default,
					estimate_only     = TRUE,
					optimization_alg  = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (!is.null(fit) && isTRUE(fit$converged)) {
				beta_hat_T = as.numeric(fit$b[j_T + 1L])
				if (is.finite(beta_hat_T) && abs(beta_hat_T) <= private$max_abs_reasonable_coef)
					return(beta_hat_T)
			}
			# fall back to glmmTMB weighted path
			for (predictors_df in private$glmm_predictors_df_candidates()) {
				mod = private$fit_weighted_glmm_on_data(predictors_df, row_weights = row_weights, se = FALSE)
				if (!private$.is_usable_glmm_fit(mod, se = FALSE)) next
				beta = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (!is.null(beta) && "w" %in% names(beta) && is.finite(beta["w"]))
					return(as.numeric(beta["w"]))
			}
			NA_real_
		},
		glmm_family        = function() stats::poisson(link = "log"),
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		get_supported_testing_types_impl = function(){
			if (isTRUE(private$use_rcpp)) c("wald", "lik_ratio") else "wald"
		},
		compute_score_two_sided_pval_impl = function(delta){
			stop(class(self)[1], " does not support score p-values.", call. = FALSE)
		},
		compute_score_confidence_interval_impl = function(alpha){
			stop(class(self)[1], " does not support score confidence intervals.", call. = FALSE)
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			stop(class(self)[1], " does not support gradient p-values.", call. = FALSE)
		},
		compute_gradient_confidence_interval_impl = function(alpha){
			stop(class(self)[1], " does not support gradient confidence intervals.", call. = FALSE)
		},
		shared = function(estimate_only = FALSE){
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				private$shared_glmm_tmb(estimate_only)
			}
			if (!estimate_only) .inflate_kk_onelik_standard_error_with_jackknife(private, self)
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
			if (ncol(as.matrix(private$X)) > 0){
				X_fit = private$create_design_matrix()
			} else {
				X_fit = cbind(`(Intercept)` = 1, w = private$w)
			}
			if ("treatment" %in% colnames(X_fit))
				colnames(X_fit)[colnames(X_fit) == "treatment"] = "w"
			X_fit = as.matrix(X_fit)
			j_T_r = which(colnames(X_fit) == "w")
			if (length(j_T_r) == 0L) j_T_r = 2L
			j_T = as.integer(j_T_r - 1L)
			
			n_params = ncol(X_fit) + 1L
			fit = tryCatch(
				fast_poisson_glmm_cpp(
					X        = X_fit,
					y        = as.numeric(private$y),
					group_id = as.integer(group_id),
					j_T      = j_T,
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					smart_cold_start = private$smart_cold_start_default,
					estimate_only    = estimate_only,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)

			if (is.null(fit) || !isTRUE(fit$converged)) {
				# Rcpp failed; fall back to glmmTMB
				return(private$shared_glmm_tmb(estimate_only = estimate_only))
			}
			beta_hat_T = as.numeric(fit$b[j_T_r])
			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				return(private$shared_glmm_tmb(estimate_only = estimate_only))
			}
			private$cached_mod = fit
			private$set_fit_warm_start(as.numeric(c(fit$b, fit$log_sigma)), "params", fisher = fit$fisher_information)
			private$cached_values$likelihood_test_context = list(
				X = X_fit,
				y = as.numeric(private$y),
				group_id = as.integer(group_id),
				j_treat = as.integer(j_T_r),
				start = as.numeric(c(fit$b, fit$log_sigma))
			)
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		},
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
					fast_poisson_glmm_cpp(
						X = X_fit,
						y = y,
						group_id = group_id,
						warm_start_params = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start,
						warm_start_fisher_info = private$get_fit_warm_start_fisher(length(ctx$start)),
						smart_cold_start = private$smart_cold_start_default,
						j_T = j_treat - 1L,
						estimate_only = FALSE,
						n_gh = 20L,
						maxit = 300L,
						eps_g = 1e-6,
						fixed_idx = j_treat,
						fixed_values = delta,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){
					as.numeric(c(fit$b, fit$log_sigma))
				},
				score = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.numeric(get_poisson_glmm_score_cpp(X_fit, y, group_id, params))
				},
				observed_information = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.matrix(get_poisson_glmm_hessian_cpp(X_fit, y, group_id, params))
				},
				fisher_information = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.matrix(get_poisson_glmm_hessian_cpp(X_fit, y, group_id, params))
				},
				information = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.matrix(get_poisson_glmm_hessian_cpp(X_fit, y, group_id, params))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		},
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$b)
			sigma_u = exp(as.numeric(null_fit$log_sigma))
			X = spec$X
			group_id = spec$group_id
			n = nrow(X)
			K = max(group_id)
			u = rnorm(K, 0, sigma_u)
			eta = as.numeric(X %*% b_null) + u[group_id]
			y_sim = as.integer(rpois(n, exp(pmin(eta, 20))))
			j = spec$j
			full_res = tryCatch(
				fast_poisson_glmm_cpp(
					X = X, y = as.numeric(y_sim), group_id = group_id,
					j_T = j - 1L,
					smart_cold_start = private$smart_cold_start_default,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$b[j])) return(NULL)
			list(
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					tryCatch(
						fast_poisson_glmm_cpp(
							X = X, y = as.numeric(y_sim), group_id = group_id,
							j_T = j - 1L,
							warm_start_params = start %||% as.numeric(c(full_res$b, full_res$log_sigma)),
							estimate_only = FALSE,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = private$smart_cold_start_default,
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
