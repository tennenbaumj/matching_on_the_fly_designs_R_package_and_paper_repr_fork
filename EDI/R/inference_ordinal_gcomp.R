#' G-Computation Mean-Difference Inference for Ordinal Responses
#'
#' Fits a proportional-odds working model for an ordinal outcome using treatment
#' and, optionally, all recorded covariates, then estimates the marginal mean
#' difference by standardizing predicted mean ranks under all-treated and
#' all-control assignments over the empirical covariate distribution.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalGCompMeanDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalGCompMeanDiff = R6::R6Class("InferenceOrdinalGCompMeanDiff",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description Initialize the g-computation (G-Comp) inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an ordinal response.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL}
		#' (default), the formula from the design object is used and its pre-computed design
		#' matrix is reused. If a formula is provided, a new design matrix is constructed from
		#' the design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		},
		#' @description Computes the g-computation (G-Comp) treatment-effect estimate (mean difference).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$md
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			# Keep weighted replicate evaluation side-effect free so it cannot poison
			# the ordinary cached fit or its warm-start / model-selection state.
			saved_cached_values = private$cached_values
			saved_best_X_colnames = private$best_X_colnames
			saved_fit_warm_start = private$fit_warm_start
			saved_fit_warm_start_type = private$fit_warm_start_type
			saved_fit_warm_start_fisher = private$fit_warm_start_fisher
			on.exit({
				private$cached_values = saved_cached_values
				private$best_X_colnames = saved_best_X_colnames
				private$fit_warm_start = saved_fit_warm_start
				private$fit_warm_start_type = saved_fit_warm_start_type
				private$fit_warm_start_fisher = saved_fit_warm_start_fisher
			}, add = TRUE)
			as.numeric(private$weighted_gcomp_md_from_row_weights(row_weights))[1L]
		},
		#' @description Computes a 1 - \code{alpha} confidence interval for the G-Comp mean difference.
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (!private$has_finite_md_se()){
				warning(
					"Ordinal G-computation: falling back to bootstrap because ",
					"delta-method standard error is unavailable."
				)
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			z_val = stats::qnorm(1 - alpha / 2)
			ci = private$cached_values$md + c(-1, 1) * z_val * private$cached_values$se_md
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		#' @description Computes a two-sided Wald p-value for the G-Comp mean difference.
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (!private$has_finite_md_se()){
				warning(
					"Ordinal G-computation: falling back to bootstrap because ",
					"delta-method standard error is unavailable."
				)
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			z_val = (private$cached_values$md - delta) / private$cached_values$se_md
			2 * stats::pnorm(-abs(z_val))
		},
		#' @description Computes a Wald confidence interval for the G-Comp mean difference.
		#' @param alpha The significance level (default 0.05).
		compute_wald_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (!private$has_finite_md_se()){
				warning(
					"Ordinal G-computation: falling back to bootstrap because ",
					"delta-method standard error is unavailable."
				)
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			z_val = stats::qnorm(1 - alpha / 2)
			ci = private$cached_values$md + c(-1, 1) * z_val * private$cached_values$se_md
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		#' @description Computes a Wald two-sided p-value for the G-Comp mean difference.
		#' @param delta The null treatment effect (default 0).
		compute_wald_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (!private$has_finite_md_se()){
				warning(
					"Ordinal G-computation: falling back to bootstrap because ",
					"delta-method standard error is unavailable."
				)
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			z_val = (private$cached_values$md - delta) / private$cached_values$se_md
			2 * stats::pnorm(-abs(z_val))
		}
	),
	private = list(
		best_X_colnames = NULL,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_X_colnames
			X_data = private$get_X()
			X_fit = if (length(X_cols) == 0L){
				matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				cbind(treatment = private$w, X_cov)
			}
			n_params = (length(sort(unique(private$y))) - 1L) + ncol(X_fit)
			fit = tryCatch(
				fast_ordinal_regression_cpp(
					X = X_fit, y = as.numeric(private$y),
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(fit) || length(fit$b) == 0 || is.null(fit$alpha)){
				return(NA_real_)
			}
			private$set_fit_warm_start(fit$params, "params", fisher = fit$fisher_information)
			res = gcomp_ordinal_proportional_odds_post_fit_cpp(
				X_fit = X_fit,
				coef_hat = as.numeric(fit$b),
				alpha_hat = as.numeric(fit$alpha),
				j_treat = 1L
			)
			as.numeric(res$md)
		},
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		},
		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (estimate_only && !is.null(private$cached_values$md)) return(invisible(NULL))
			X_full = private$build_design_matrix()
			X_fit = X_full[, -1, drop = FALSE]
			j_treat = 1
			fit = tryCatch(
				if (estimate_only) {
					n_params = ncol(X_fit) + length(sort(unique(private$y))) - 1L
					fast_ordinal_regression_cpp(
						X = X_fit,
						y = as.numeric(private$y),
						warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
						smart_cold_start = private$smart_cold_start_default
					)
				} else {
					fast_ordinal_regression_with_var_cpp(X = X_fit, y = as.numeric(private$y))
				},
				error = function(e) NULL
			)
			if (is.null(fit) || length(fit$b) == 0 || is.null(fit$alpha)){
				private$cached_values$mean1 = NA_real_
				private$cached_values$mean0 = NA_real_
				private$cached_values$md = NA_real_
				private$cached_values$se_md = NA_real_
				return(invisible(NULL))
			}
			private$best_X_colnames = setdiff(colnames(X_fit), "treatment")
			coef_hat = as.numeric(fit$b)
			alpha_hat = as.numeric(fit$alpha)
			res = gcomp_ordinal_proportional_odds_post_fit_cpp(
				X_fit = X_fit,
				coef_hat = coef_hat,
				alpha_hat = alpha_hat,
				j_treat = j_treat
			)
			private$cached_values$mean1 = res$mean1
			private$cached_values$mean0 = res$mean0
			private$cached_values$md = res$md
			if (estimate_only) {
				private$cached_values$beta_hat_T = private$cached_values$md
				private$cached_values$s_beta_hat_T = NA_real_
				private$set_fit_warm_start(fit$params, "params", fisher = fit$fisher_information)
				return(invisible(NULL))
			}
			theta = c(alpha_hat, coef_hat)
			vcov_mat = if (!is.null(fit$vcov)) as.matrix(fit$vcov) else NULL
			if (!is.null(vcov_mat) && length(theta) > 0L &&
				 nrow(vcov_mat) == length(theta) && ncol(vcov_mat) == length(theta)){
				grad = private$compute_md_gradient(X_fit, theta, length(alpha_hat), j_treat)
				var_md = as.numeric(crossprod(grad, vcov_mat %*% grad))
				private$cached_values$se_md = if (is.finite(var_md) && var_md > 0) sqrt(var_md) else NA_real_
			} else {
				private$cached_values$se_md = NA_real_
			}
			private$cached_values$beta_hat_T = private$cached_values$md
			private$cached_values$s_beta_hat_T = private$cached_values$se_md
		},
		compute_md_gradient = function(X_fit, theta, n_alpha, j_treat, base_step = 1e-6){
			n_params = length(theta)
			grad = numeric(n_params)
			for (j in seq_len(n_params)){
				step = max(base_step, base_step * (1 + abs(theta[j])))
				theta_plus = theta
				theta_minus = theta
				theta_plus[j] = theta[j] + step
				theta_minus[j] = theta[j] - step
				md_plus = private$compute_md_from_theta(X_fit, theta_plus, n_alpha, j_treat)
				md_minus = private$compute_md_from_theta(X_fit, theta_minus, n_alpha, j_treat)
				grad[j] = (md_plus - md_minus) / (2 * step)
			}
			grad
		},
		compute_md_from_theta = function(X_fit, theta, n_alpha, j_treat){
			alpha_vec = theta[seq_len(n_alpha)]
			coef_vec = theta[(n_alpha + 1):length(theta)]
			as.numeric(
				gcomp_ordinal_proportional_odds_post_fit_cpp(
					X_fit = X_fit,
					coef_hat = coef_vec,
					alpha_hat = alpha_vec,
					j_treat = j_treat
				)$md
			)
		},
		weighted_gcomp_md_from_row_weights = function(row_weights){
			X_full = private$build_design_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			if (is.null(reduced$X) || !is.finite(reduced$j_treat)) return(NA_real_)
			X_fit = reduced$X[, -1, drop = FALSE]
			j_treat = reduced$j_treat - 1L
			ok = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
			if (sum(ok) <= ncol(X_fit)) return(NA_real_)
			X_fit_ok = X_fit[ok, , drop = FALSE]
			y_ok = as.numeric(private$y[ok])
			n_params = ncol(X_fit_ok) + length(sort(unique(y_ok))) - 1L
			try_weighted_fit = function(start_params = NULL, start_fisher = NULL){
				tryCatch(
					fast_ordinal_regression_weighted_cpp(
						X = X_fit_ok,
						y = y_ok,
						weights = as.numeric(row_weights[ok]),
						warm_start_params = start_params,
						warm_start_fisher_info = start_fisher,
						smart_cold_start = private$smart_cold_start_default
					),
					error = function(e) NULL
				)
			}
			mod = NULL
			start_candidates = list(
				list(
					params = private$get_fit_warm_start_for_length("params", n_params),
					fisher = private$get_fit_warm_start_fisher(n_params)
				)
			)
			if (is.null(start_candidates[[1L]]$params)) {
				unweighted_fit = tryCatch(
					fast_ordinal_regression_cpp(
						X = X_fit_ok,
						y = y_ok,
						warm_start_params = NULL,
						warm_start_fisher_info = NULL,
						smart_cold_start = private$smart_cold_start_default
					),
					error = function(e) NULL
				)
				if (!is.null(unweighted_fit) && !is.null(unweighted_fit$params)) {
					start_candidates[[length(start_candidates) + 1L]] = list(
						params = as.numeric(unweighted_fit$params),
						fisher = unweighted_fit$fisher_information %||% NULL
					)
				}
			}
			start_candidates[[length(start_candidates) + 1L]] = list(params = NULL, fisher = NULL)
			for (cand in start_candidates) {
				mod = try_weighted_fit(cand$params, cand$fisher)
				if (!is.null(mod) && !is.null(mod$b) && length(mod$b) > 0L && !is.null(mod$alpha)) break
			}
			if (is.null(mod)) return(NA_real_)
			coef_hat = as.numeric(mod$b)
			alpha_hat = as.numeric(mod$alpha)
			if (!length(coef_hat) || !length(alpha_hat) || j_treat < 1L || j_treat > length(coef_hat)) return(NA_real_)
			private$set_fit_warm_start(as.numeric(mod$params), "params", fisher = mod$fisher_information)
			private$best_X_colnames = setdiff(colnames(X_fit_ok), "treatment")
			res = tryCatch(
				gcomp_ordinal_proportional_odds_post_fit_cpp(
					X_fit = X_fit_ok,
					coef_hat = coef_hat,
					alpha_hat = alpha_hat,
					j_treat = j_treat
				),
				error = function(e) NULL
			)
			if (is.null(res)) return(NA_real_)
			as.numeric(res$md)
		},
		has_finite_md_se = function(){
			se = private$cached_values$se_md
			is.finite(se) && se > 0
		},
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T %||% private$cached_values$se_md %||% NA_real_
		},
		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df %||% NA_real_
		}
	)
)
