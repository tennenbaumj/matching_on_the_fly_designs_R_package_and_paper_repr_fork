#' Stereotype Logit Regression Inference for Ordinal Responses
#'
#' Fits a stereotype logit regression for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @name InferenceOrdinalStereotypeLogitRegr
#' @export
InferenceOrdinalStereotypeLogitRegr = R6::R6Class("InferenceOrdinalStereotypeLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a stereotype logit inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		#' @param smart_cold_start_default Whether to use smart cold starts.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Recomputes the stereotype-logit treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			X_fit = private$build_design_matrix()
			if (!is.null(private$best_Xmm_colnames)) {
				keep = c("treatment", intersect(private$best_Xmm_colnames, colnames(X_fit)))
				X_fit = X_fit[, keep, drop = FALSE]
			}
			fit = weighted_ordinal_bootstrap_surrogate_fit(X_fit, private$y, row_weights, method = "logistic")
			if (is.null(fit) || !is.finite(fit$beta_hat)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$summary_table = NULL
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		supports_likelihood_tests = function(){ TRUE },
		best_Xmm_colnames = NULL,
		get_complexity_tier = function() "heavy",
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			if (length(X_cols) == 0L){
				X = as.matrix(private$w)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}

			n_params = (length(sort(unique(private$y))) - 1L) + ncol(X)
			if (length(sort(unique(private$y))) >= 3) n_params = n_params + (length(sort(unique(private$y))) - 2L)
			
			ws_args = private$get_backend_warm_start_args(n_params)
			ws_fisher = ws_args$warm_start_fisher_info
			res = tryCatch(
				fast_stereotype_logit_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_params = ws_args$start_params,
					warm_start_fisher_info = ws_fisher,
					estimate_only = TRUE
				),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[1])){
				return(NA_real_)
			}
			private$set_fit_warm_start(as.numeric(res$params), "params", fisher = ws_fisher)
			as.numeric(res$b[1])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		get_bootstrap_worker_spec = function(){
			private$shared(estimate_only = FALSE)
			list(
				X_full = private$build_design_matrix(),
				best_X_cols = private$best_Xmm_colnames,
				fit_fun = function(X_fit, keep){
					K = length(sort(unique(private$y)))
					n_params = (K - 1L) + ncol(X_fit)
					if (K >= 3) n_params = n_params + (K - 2L)
					ws_args = private$get_backend_warm_start_args(n_params)
					res = fast_stereotype_logit_cpp(
						X_fit, private$y,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info
					)
					list(b = res$b, ssq_b_j = NA_real_, params = res$params, fisher_information = res$fisher_information)
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					K = length(sort(unique(private$y)))
					n_params = (K - 1L) + ncol(X_fit)
					if (K >= 3) n_params = n_params + (K - 2L)
					ws_args = private$get_backend_warm_start_args(n_params)
					if (estimate_only) {
						res = fast_stereotype_logit_cpp(
							X_fit, private$y,
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(b = res$b, ssq_b_j = NA_real_, params = res$params, fisher_information = res$fisher_information)
					} else {
						res = fast_stereotype_logit_with_var_cpp(
							X_fit, private$y,
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(b = res$b, ssq_b_j = res$ssq_b_1, params = res$params, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 1L || !is.finite(mod$b[1])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$params, "params", fisher = attempt$fit$fisher_information)
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), "treatment")
				if (!estimate_only) {
					n_alpha = length(sort(unique(private$y))) - 1L
					private$cached_values$likelihood_test_context = list(
						X = attempt$X,
						j_treat = as.integer(n_alpha + 1L),
						full_params = attempt$fit$params,
						full_neg_loglik = attempt$fit$neg_loglik
					)
				} else {
					private$cached_values$likelihood_test_context = NULL
				}
				list(b = c(0, attempt$fit$b[1]), ssq_b_2 = attempt$fit$ssq_b_j)
			} else {
				private$cached_values$likelihood_test_context = NULL
				NULL
			}
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			full_fit = list(params = ctx$full_params, neg_loglik = ctx$full_neg_loglik)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					n_params = length(ctx$full_params)
					res = tryCatch(
						fast_stereotype_logit_cpp(
							X_fit, y,
							fixed_idx = j_treat, fixed_values = delta,
							warm_start_params = start %||% private$get_fit_warm_start_for_length("params", n_params),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0L) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
				},
				extract_start = function(fit){ as.numeric(fit$params) },
				score = function(fit){
					get_stereotype_logit_score_cpp(X_fit, y, as.numeric(fit$params))
				},
				observed_information = function(fit){
					-get_stereotype_logit_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				fisher_information = function(fit){
					fit$fisher_information %||% (-get_stereotype_logit_hessian_cpp(X_fit, y, as.numeric(fit$params)))
				},
				information = function(fit){
					fit$information %||% fit$fisher_information %||% (-get_stereotype_logit_hessian_cpp(X_fit, y, as.numeric(fit$params)))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			params_null = as.numeric(null_fit$b)
			n_params    = length(params_null)
			K           = as.integer(spec$K)
			n_alpha     = K - 1L
			n_gamma     = K - 2L
			p           = n_params - n_alpha - n_gamma
			X_fit       = spec$X
			j           = spec$j
			n           = nrow(X_fit)
			if (p < 1L) return(NULL)

			alphas = params_null[seq_len(n_alpha)]
			betas  = params_null[(n_alpha + 1L):(n_alpha + p)]
			gammas = params_null[(n_alpha + p + 1L):n_params]

			eg   = exp(gammas)
			seg  = sum(eg)
			phi  = c(0, cumsum(eg) / (1 + seg))

			eta  = as.numeric(X_fit %*% betas)
			y_sim = integer(n)
			for (i in seq_len(n)){
				logits = c(0, alphas + phi[seq(2L, K)] * eta[i])
				cum_p  = plogis(logits)
				cat_p  = pmax(c(cum_p[1L], diff(cum_p), 1 - cum_p[n_alpha + 1L]), 0)
				s = sum(cat_p)
				y_sim[i] = if (is.finite(s) && s > 0) sample.int(K, 1L, prob = cat_p / s) else 1L
			}
			if (length(unique(y_sim)) < K) return(NULL)

			ws   = private$get_fit_warm_start_for_length("params", n_params) %||% params_null
			full = tryCatch(
				fast_stereotype_logit_cpp(
					X = X_fit, y = y_sim, K = K,
					j_T = 0L, estimate_only = FALSE,
					warm_start_params = ws,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
				),
				error = function(e) NULL
			)
			if (is.null(full) || !isTRUE(full$converged)) return(NULL)
			list(
				full_fit = full,
				fit_null = function(d, start = NULL){
					ws2 = start %||% private$get_fit_warm_start_for_length("params", n_params) %||% params_null
					f2  = tryCatch(
						fast_stereotype_logit_cpp(
							X = X_fit, y = y_sim, K = K,
							j_T = 0L, estimate_only = FALSE,
							warm_start_params = ws2,
							warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
							fixed_idx = j, fixed_values = d
						),
						error = function(e) NULL
					)
					if (is.null(f2) || !isTRUE(f2$converged)) return(NULL)
					f2
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik %||% fit$neg_ll)
			)
		},
				build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X = cbind(treatment = private$w, X_cov)
			}
			X
		}
	)
)

#' Continuation Ratio Regression Inference for Ordinal Responses
#'
#' Fits a continuation ratio regression for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @description Fits a continuation-ratio ordinal regression and reports the
#'   treatment effect estimate on the model's coefficient scale.
#' @export
InferenceOrdinalContRatioRegr = R6::R6Class("InferenceOrdinalContRatioRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a continuation ratio inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		#' @param smart_cold_start_default Whether to use smart cold starts.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Recomputes the continuation-ratio treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			X_fit = private$build_design_matrix()
			if (!is.null(private$best_Xmm_colnames)) {
				keep = c("treatment", intersect(private$best_Xmm_colnames, colnames(X_fit)))
				X_fit = X_fit[, keep, drop = FALSE]
			}
			fit = weighted_ordinal_bootstrap_surrogate_fit(X_fit, private$y, row_weights, method = "logistic")
			if (is.null(fit) || !is.finite(fit$beta_hat)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$summary_table = NULL
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		supports_likelihood_tests = function(){ TRUE },
		best_Xmm_colnames = NULL,
		get_complexity_tier = function() "heavy",
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			if (length(X_cols) == 0L){
				X = as.matrix(private$w)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}

			n_params = (length(sort(unique(private$y))) - 1L) + ncol(X)
			ws_args = private$get_backend_warm_start_args(n_params)
			res = tryCatch(
				fast_continuation_ratio_regression_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_fisher_info = ws_args$warm_start_fisher_info
				),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[length(res$b)])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$fisher_information)
			as.numeric(res$b[length(res$b)])
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					n_params = ncol(X_fit) + length(sort(unique(private$y))) - 1L
					ws_args = private$get_backend_warm_start_args(n_params)
					if (estimate_only) {
						res = fast_continuation_ratio_regression_cpp(
							X_fit, private$y,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(b = res$b, ssq_b_j = NA_real_, fisher_information = res$fisher_information)
					} else {
						res = fast_continuation_ratio_regression_with_var_cpp(
							X_fit, private$y,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(
							b = res$b,
							ssq_b_j = res$ssq_b_j %||% res$ssq_b_1,
							ssq_b_2 = res$ssq_b_2 %||% res$ssq_b_j %||% res$ssq_b_1,
							params = res$params,
							neg_loglik = res$neg_loglik,
							fisher_information = res$fisher_information
						)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 1L || !is.finite(mod$b[1])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$b, "beta", fisher = attempt$fit$fisher_information)
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), "treatment")
				if (!estimate_only) {
					n_alpha = length(sort(unique(private$y))) - 1L
					private$cached_values$likelihood_test_context = list(
						X = attempt$X,
						j_treat = as.integer(n_alpha + 1L),
						full_params = attempt$fit$params,
						full_neg_loglik = attempt$fit$neg_loglik
					)
				} else {
					private$cached_values$likelihood_test_context = NULL
				}
				list(
					b = c(0, attempt$fit$b[1]),
					ssq_b_2 = attempt$fit$ssq_b_2 %||% attempt$fit$ssq_b_j
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
				NULL
			}
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			full_fit = list(params = ctx$full_params, neg_loglik = ctx$full_neg_loglik)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					n_params = length(ctx$full_params)
					res = tryCatch(
						fast_continuation_ratio_regression_cpp(
							X_fit, y,
							fixed_idx = j_treat, fixed_values = delta,
							warm_start_beta = start %||% private$get_fit_warm_start_for_length("beta", n_params),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0L) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
				},
				extract_start = function(fit){ as.numeric(fit$params) },
				score = function(fit){
					get_continuation_ratio_regression_score_cpp(X_fit, y, as.numeric(fit$params))
				},
				observed_information = function(fit){
					-get_continuation_ratio_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				fisher_information = function(fit){
					fit$fisher_information %||% (-get_continuation_ratio_regression_hessian_cpp(X_fit, y, as.numeric(fit$params)))
				},
				information = function(fit){
					fit$information %||% fit$fisher_information %||% (-get_continuation_ratio_regression_hessian_cpp(X_fit, y, as.numeric(fit$params)))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			params_null = as.numeric(null_fit$b)
			n_params    = length(params_null)
			K           = as.integer(spec$K)
			n_alpha     = K - 1L
			X_fit       = spec$X
			j           = spec$j
			n           = nrow(X_fit)
			p           = n_params - n_alpha

			alphas = params_null[seq_len(n_alpha)]
			betas  = params_null[(n_alpha + 1L):n_params]
			eta    = as.numeric(X_fit %*% betas)

			y_sim = integer(n)
			for (i in seq_len(n)){
				k = 1L
				assigned = FALSE
				while (k <= n_alpha){
					if (runif(1) < plogis(alphas[k] + eta[i])){
						y_sim[i] = k
						assigned  = TRUE
						break
					}
					k = k + 1L
				}
				if (!assigned) y_sim[i] = K
			}
			if (length(unique(y_sim)) < K) return(NULL)

			ws   = private$get_fit_warm_start_for_length("beta", n_params) %||% params_null
			full = tryCatch(
				fast_continuation_ratio_regression_cpp(
					X = X_fit, y = y_sim, K = K,
					j_T = 0L, estimate_only = FALSE,
					warm_start_beta = ws,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
				),
				error = function(e) NULL
			)
			if (is.null(full) || !isTRUE(full$converged)) return(NULL)
			list(
				full_fit = full,
				fit_null = function(d, start = NULL){
					ws2 = start %||% private$get_fit_warm_start_for_length("beta", n_params) %||% params_null
					f2  = tryCatch(
						fast_continuation_ratio_regression_cpp(
							X = X_fit, y = y_sim, K = K,
							j_T = 0L, estimate_only = FALSE,
							warm_start_beta = ws2,
							warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
							fixed_idx = j, fixed_values = d
						),
						error = function(e) NULL
					)
					if (is.null(f2) || !isTRUE(f2$converged)) return(NULL)
					f2
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik %||% fit$neg_ll)
			)
		},
				build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X = cbind(treatment = private$w, X_cov)
			}
			X
		}
	)
)
