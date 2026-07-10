#' Poisson Regression Inference for Count Responses
#'
#' Fits a Poisson log-link regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountPoisson$new(seq_des)
#' inf$compute_estimate()
#' }
#' \donttest{
#' inf$set_seed(1)
#' inf$compute_lik_ratio_bootstrap_two_sided_pval(delta = 0, B = 9, show_progress = FALSE)
#' }
#' @export
InferenceCountPoisson = R6::R6Class("InferenceCountPoisson",
	lock_objects = FALSE,
	inherit = InferenceCountLikelihood,
	public = list(
		#' @description Initialize a Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose               Whether to print progress messages.
		#' @param harden                Whether to apply robustness measures.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL, harden = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes an asymptotic confidence interval using the configured test.
		#' @param alpha Significance level. Default 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			switch(
				private$testing_type,
				wald = self$compute_wald_confidence_interval(alpha = alpha),
				score = self$compute_score_confidence_interval(alpha = alpha),
				gradient = self$compute_gradient_confidence_interval(alpha = alpha),
				lik_ratio = self$compute_lik_ratio_confidence_interval(alpha = alpha)
			)
		},
		#' @description Computes an asymptotic two-sided p-value using the configured test.
		#' @param delta Null treatment effect. Default 0.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			switch(
				private$testing_type,
				wald = self$compute_wald_two_sided_pval(delta = delta),
				score = self$compute_score_two_sided_pval(delta = delta),
				gradient = self$compute_gradient_two_sided_pval(delta = delta),
				lik_ratio = self$compute_lik_ratio_two_sided_pval(delta = delta)
			)
		},
		#' @description Computes a design-conservative Wald confidence interval.
		#' @param alpha Significance level. Default 0.05.
		compute_wald_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			ci_model = private$compute_wald_confidence_interval_impl(alpha)
			private$design_conservative_ci(ci_model, alpha = alpha)
		},
		#' @description Computes a design-conservative Wald two-sided p-value.
		#' @param delta Null treatment effect. Default 0.
		compute_wald_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			p_model = private$compute_wald_two_sided_pval_impl(delta)
			private$design_conservative_pval(p_model, delta = delta)
		},
		#' @description Computes a design-conservative score confidence interval.
		#' @param alpha Significance level. Default 0.05.
		compute_score_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			ci_model = private$compute_score_confidence_interval_impl(alpha)
			private$design_conservative_ci(ci_model, alpha = alpha)
		},
		#' @description Computes a design-conservative score two-sided p-value.
		#' @param delta Null treatment effect. Default 0.
		compute_score_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			p_model = private$compute_score_two_sided_pval_impl(delta)
			private$design_conservative_pval(p_model, delta = delta)
		},
		#' @description Computes a design-conservative likelihood-ratio confidence interval.
		#' @param alpha Significance level. Default 0.05.
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			ci_model = private$compute_lik_ratio_confidence_interval_impl(alpha)
			private$design_conservative_ci(ci_model, alpha = alpha)
		},
		#' @description Computes a design-conservative likelihood-ratio two-sided p-value.
		#' @param delta Null treatment effect. Default 0.
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			p_model = private$compute_lik_ratio_two_sided_pval_impl(delta)
			private$design_conservative_pval(p_model, delta = delta)
		},
		#' @description Computes a design-conservative gradient confidence interval.
		#' @param alpha Significance level. Default 0.05.
		compute_gradient_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			ci_model = private$compute_gradient_confidence_interval_impl(alpha)
			private$design_conservative_ci(ci_model, alpha = alpha)
		},
		#' @description Computes a design-conservative gradient two-sided p-value.
		#' @param delta Null treatment effect. Default 0.
		compute_gradient_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			p_model = private$compute_gradient_two_sided_pval_impl(delta)
			private$design_conservative_pval(p_model, delta = delta)
		},
		#' @description Computes a design-conservative parametric LR-bootstrap p-value.
		#' @param delta Null treatment effect. Default 0.
		#' @param B Number of bootstrap replicates.
		#' @param show_progress Whether to show progress.
		#' @param min_number_usable_samples Minimum usable bootstrap samples.
		#' @param max_attempts_per_replicate Maximum attempts per replicate.
		compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			p_model = super$compute_lik_ratio_bootstrap_two_sided_pval(
				delta = delta,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate
			)
			private$design_conservative_pval(p_model, delta = delta)
		},
		#' @description Computes a design-conservative parametric LR-bootstrap CI.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of bootstrap replicates.
		#' @param show_progress Whether to show progress.
		#' @param min_number_usable_samples Minimum usable bootstrap samples.
		#' @param max_attempts_per_replicate Maximum attempts per replicate.
		#' @param root_tolerance Root tolerance.
		#' @param max_root_iterations Maximum root iterations.
		compute_lik_ratio_bootstrap_confidence_interval = function(alpha = 0.05, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L, root_tolerance = NULL, max_root_iterations = 8L){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			ci_model = super$compute_lik_ratio_bootstrap_confidence_interval(
				alpha = alpha,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate,
				root_tolerance = root_tolerance,
				max_root_iterations = max_root_iterations
			)
			private$design_conservative_ci(ci_model, alpha = alpha)
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_data = private$get_X()
			if (is.null(X_data) || ncol(X_data) == 0) {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					res = fast_poisson_regression_weighted_cpp(
						X = X_fit,
						y = private$y,
						weights = row_weights,
						warm_start_beta = ws_args$warm_start_beta,
						warm_start_weights = ws_args$warm_start_weights,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default
					)
					XtWX = res$XtWX
					ssq_b_j = NA_real_
					if (!estimate_only && !is.null(XtWX) && is.matrix(XtWX) && j_treat <= nrow(XtWX)) {
						inv_XtWX = tryCatch(solve(XtWX), error = function(e) NULL)
						if (!is.null(inv_XtWX) && is.finite(inv_XtWX[j_treat, j_treat]) && inv_XtWX[j_treat, j_treat] > 0) ssq_b_j = inv_XtWX[j_treat, j_treat]
					}
					list(b = res$b, XtWX = XtWX, ssq_b_j = ssq_b_j, j_treat = j_treat)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && length(mod$b) >= mod$j_treat && is.finite(mod$b[mod$j_treat])
				}
			)
			private$cached_mod = attempt$fit
			if (is.null(attempt$fit) || is.null(attempt$fit$b)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			j_treat = attempt$fit$j_treat %||% 2L
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[j_treat])
			ssq = attempt$fit$ssq_b_j
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
			private$set_fit_warm_start(
				as.numeric(attempt$fit$b),
				"beta",
				fisher = attempt$fit$XtWX
			)
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		best_X_colnames = NULL,
		poisson_X_full_cache = NULL,
		poisson_w_cache = NULL,
		get_complexity_tier = function() "medium",
		design_jackknife_pval = function(delta = 0){
			tryCatch({
				p = self$compute_jackknife_wald_two_sided_pval(delta = delta)
				p = as.numeric(p)[1L]
				if (is.finite(p) && p >= 0 && p <= 1) p else NA_real_
			}, error = function(e) NA_real_)
		},
		design_jackknife_ci = function(alpha = 0.05){
			tryCatch({
				ci = as.numeric(self$compute_jackknife_wald_confidence_interval(alpha = alpha))
				if (length(ci) >= 2L && all(is.finite(ci[1:2])) && ci[1L] <= ci[2L]) ci[1:2] else c(NA_real_, NA_real_)
			}, error = function(e) c(NA_real_, NA_real_))
		},
		design_conservative_pval = function(model_p, delta = 0){
			model_p = as.numeric(model_p)[1L]
			design_p = private$design_jackknife_pval(delta = delta)
			if (is.finite(model_p) && is.finite(design_p)) return(max(model_p, design_p))
			if (is.finite(model_p)) return(model_p)
			if (is.finite(design_p)) return(design_p)
			NA_real_
		},
		design_conservative_ci = function(model_ci, alpha = 0.05){
			model_ci = as.numeric(model_ci)
			if (length(model_ci) < 2L) model_ci = c(NA_real_, NA_real_)
			model_ci = model_ci[1:2]
			design_ci = private$design_jackknife_ci(alpha = alpha)
			ci = if (all(is.finite(model_ci)) && all(is.finite(design_ci))) {
				c(min(model_ci[1L], design_ci[1L]), max(model_ci[2L], design_ci[2L]))
			} else if (all(is.finite(model_ci))) {
				model_ci
			} else if (all(is.finite(design_ci))) {
				design_ci
			} else {
				c(NA_real_, NA_real_)
			}
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_X_colnames
			X_data = private$get_X()

			if (length(X_cols) == 0L){
				X = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(1, treatment = private$w, X_cov)
			}
			ws_args = private$get_backend_warm_start_args(ncol(X))
			res = tryCatch(
				fast_poisson_regression_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_weights = ws_args$warm_start_weights,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default,
					estimate_only = TRUE
				),
				error = function(e) NULL
			)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$XtWX)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_lik_ratio_param_bootstrap = function(){
			TRUE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		supports_fisher_information = function(){
			TRUE
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null     = as.numeric(null_fit$b)
			mu         = pmax(exp(as.numeric(spec$X %*% b_null)), 0)
			y_sim      = private$simulate_param_boot_poisson_y(mu)
			if (is.null(y_sim)) return(NULL)
			X_fit      = spec$X
			j          = spec$j

			# Parametric bootstrap: use observed fit as anchor
			ws_args = private$get_backend_warm_start_args(ncol(X_fit))
			full_fit_b = tryCatch(
				fast_poisson_regression_cpp(
					X = X_fit, y = y_sim,
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_weights = ws_args$warm_start_weights,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(full_fit_b) || length(full_fit_b$b) < j || !is.finite(full_fit_b$b[j])) return(NULL)
			list(
				worker_data = list(y = y_sim),
				full_fit = full_fit_b,
				fit_null = function(d, start = NULL){
					ws_args_null = private$get_backend_warm_start_args(ncol(X_fit))
					tryCatch(
						fast_poisson_regression_with_var_cpp(
							X = X_fit, y = y_sim, j = j,
							warm_start_beta = start %||% full_fit_b$b,
							warm_start_weights = ws_args_null$warm_start_weights,
							warm_start_fisher_info = ws_args_null$warm_start_fisher_info,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit){
					eta_f = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y_sim * eta_f - exp(eta_f) - lgamma(y_sim + 1))
				}
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					fast_poisson_regression_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						warm_start_beta = start %||% ws_args$warm_start_beta,
						warm_start_weights = ws_args$warm_start_weights,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_cold_start = private$smart_cold_start_default
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_poisson_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y * eta - exp(eta) - lgamma(y + 1))
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			if (is.null(private$poisson_X_full_cache) || !identical(private$w, private$poisson_w_cache)) {
				X_data = private$get_X()
				private$poisson_X_full_cache = if (is.null(X_data) || ncol(X_data) == 0) {
					cbind(`(Intercept)` = 1, treatment = private$w)
				} else {
					cbind(`(Intercept)` = 1, treatment = private$w, X_data)
				}
				private$poisson_w_cache = private$w
			}
			X_full = private$poisson_X_full_cache
			
			if (!private$harden) {
				ws_args = private$get_backend_warm_start_args(ncol(X_full))
				if (estimate_only) {
					res = fast_poisson_regression_cpp(
						X = X_full, y = private$y,
						warm_start_beta = ws_args$warm_start_beta,
						warm_start_weights = ws_args$warm_start_weights,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = TRUE
					)
					res$ssq_b_j = NA_real_
					res$j_treat = 2L
				} else {
					res = fast_poisson_regression_with_var_cpp(
						X = X_full, y = private$y, j = 2L,
						warm_start_beta = ws_args$warm_start_beta,
						warm_start_weights = ws_args$warm_start_weights,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default
					)
					res$j_treat = 2L
				}
				private$best_X_colnames = setdiff(colnames(X_full), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(X = X_full, j_treat = 2L, full_neg_loglik = res$neg_ll)
				return(res)
			}
			
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					if (estimate_only) {
						res = fast_poisson_regression_cpp(
							X = X_fit, y = private$y,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_weights = ws_args$warm_start_weights,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = TRUE
						)
						list(b = res$b, XtWX = res$XtWX, w = res$w, ssq_b_j = NA_real_, j_treat = j_treat, neg_ll = res$neg_ll)
					} else {
						res = fast_poisson_regression_with_var_cpp(
							X = X_fit, y = private$y, j = j_treat,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_weights = ws_args$warm_start_weights,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default
						)
						res$j_treat = j_treat
						res
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod)) return(FALSE)
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = which(attempt$keep == 2L),
					full_neg_loglik = attempt$fit$neg_ll
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			if (is.null(w_mat)) return(NULL)
			X_covars = private$X
			log_transform = transform_responses == "log"
			compute_poisson_distr_parallel_cpp(X_covars, as.numeric(y), w_mat, as.numeric(delta), log_transform, private$n_cpp_threads(ncol(w_mat)))
		}
	)
)
