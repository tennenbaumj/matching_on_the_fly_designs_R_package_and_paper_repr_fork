#' Risk Difference Inference for Incidence Responses
#'
#' Fits a linear probability model via OLS for binary (incidence) responses using
#' the treatment indicator and, optionally, all recorded covariates as
#' predictors. The treatment effect is reported as a risk difference.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidRiskDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidRiskDiff = R6::R6Class("InferenceIncidRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCacheNoParamBootstrap,
	public = list(
				
		#' @description Initialize a risk-difference inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
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
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			if (estimate_only) {
				if (!is.null(private$cached_values$beta_hat_T)) return(private$cached_values$beta_hat_T)
				if (isFALSE(private$harden)) {
					X = private$build_design_matrix()
					fit = stats::lm.fit(x = X, y = as.numeric(private$y))
					b = stats::coef(fit)
					j = 2L
					private$cached_values$beta_hat_T = if (length(b) >= j && is.finite(b[j])) {
						private$best_X_colnames = setdiff(colnames(X), c("(Intercept)", "treatment"))
						as.numeric(b[j])
					} else NA_real_
					return(private$cached_values$beta_hat_T)
				}
			}
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
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				fit_fun = function(X_fit, keep){
					w_fit = row_weights
					ok = is.finite(w_fit) & w_fit > 0 & is.finite(private$y)
					if (sum(ok) <= ncol(X_fit)) return(NULL)
					stats::lm.wfit(
						x = X_fit[ok, , drop = FALSE],
						y = as.numeric(private$y[ok]),
						w = as.numeric(w_fit[ok])
					)
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = which(keep == 2L)
					!is.null(mod) &&
						is.finite(j_treat) &&
						length(mod$coefficients) >= j_treat &&
						is.finite(mod$coefficients[j_treat])
				}
			)
			if (is.null(attempt$fit)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			j_treat = which(attempt$keep == 2L)
			private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$coefficients[j_treat])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		best_X_colnames = NULL,
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
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
			res = tryCatch(fast_ols_cpp(X = X, y = as.numeric(private$y)), error = function(e) NULL)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		get_supported_testing_types_impl = function(){
			"wald"
		},
		generate_mod = function(estimate_only = FALSE){
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = private$build_design_matrix(),
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					if (estimate_only) {
						res = stats::lm.fit(x = X_fit, y = as.numeric(private$y))
						b = as.numeric(stats::coef(res))
						list(b = b, beta_hat_T = as.numeric(b[j_treat]), ssq_b_j = NA_real_, j_treat = j_treat)
					} else {
						res = fast_ols_with_var_cpp(X = X_fit, y = private$y, j = j_treat)
						res$j_treat = j_treat
						res$beta_hat_T = as.numeric(res$b[j_treat])
						res
					}
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
			}
			attempt$fit
		}
	)
	)
