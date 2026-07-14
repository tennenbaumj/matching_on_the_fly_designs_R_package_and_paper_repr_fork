#' Fractional Logit Inference for Proportion Responses
#'
#' Fits a fractional logistic regression (quasi-binomial) for proportion responses
#' (constrained to [0, 1]) using the treatment indicator and, optionally, all
#' recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropFractionalLogit$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropFractionalLogit = R6::R6Class("InferencePropFractionalLogit",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCacheNoParamBootstrap,
	public = list(
		#' @description Initialize a fractional-logit inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param harden  		Whether to apply robustness measures.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, harden = TRUE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
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
					fit = glm.fit(x = X, y = as.numeric(private$y), family = quasibinomial())
					b = fit$coefficients
					private$cached_values$beta_hat_T = if (length(b) >= 2L && is.finite(b[2L])) {
						private$best_X_colnames = setdiff(colnames(X), c("(Intercept)", "treatment"))
						as.numeric(b[2L])
					} else NA_real_
					return(private$cached_values$beta_hat_T)
				}
			}
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					res = tryCatch(
						fast_logistic_regression_weighted_cpp(
							X = X_fit,
							y = as.numeric(private$y),
							weights = as.numeric(row_weights),
							warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					list(b = res$b, fisher_information = res$fisher_information, ssq_b_2 = NA_real_)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && length(mod$b) >= 2L && is.finite(mod$b[2L])
				}
			)
			private$cached_mod = attempt$fit
			if (is.null(attempt$fit) || is.null(attempt$fit$b) || length(attempt$fit$b) < 2L || !is.finite(attempt$fit$b[2L])) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[2L])
			private$cached_values$s_beta_hat_T = NA_real_
			private$set_fit_warm_start(
				as.numeric(attempt$fit$b),
				"beta",
				fisher = attempt$fit$fisher_information
			)
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		best_X_colnames = NULL,
		supports_likelihood_tests = function(){
			FALSE
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
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			res = fast_logistic_regression_cpp(X = X, y = as.numeric(private$y), estimate_only = TRUE)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			
			attempt = if (private$harden) {
				private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = 2L,
					fit_fun = function(X_fit){
						ws_args = private$get_backend_warm_start_args(ncol(X_fit))
						if (estimate_only) {
							res = fast_logistic_regression_cpp(
								X = X_fit, 
								y = private$y, 
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_weights = ws_args$warm_start_weights,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = TRUE
							)
							list(b = res$b, ssq_b_2 = NA_real_)
						} else {
							fast_logistic_regression_with_var_cpp(
								X = X_fit, 
								y = private$y,
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_weights = ws_args$warm_start_weights,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							)
						}
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
						if (max(abs(mod$b), na.rm = TRUE) > 100) return(FALSE)
						if (estimate_only) return(TRUE)
						is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
					}
				)
			} else {
				list(
					X = X_full,
					keep = seq_len(ncol(X_full)),
					fit = {
						ws_args = private$get_backend_warm_start_args(ncol(X_full))
						if (estimate_only) {
							res = fast_logistic_regression_cpp(
								X = X_full, 
								y = private$y, 
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_weights = ws_args$warm_start_weights,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = TRUE
							)
							list(b = res$b, ssq_b_2 = NA_real_)
						} else {
							fast_logistic_regression_with_var_cpp(
								X = X_full, 
								y = private$y,
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_weights = ws_args$warm_start_weights,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							)
						}
					}
				)
			}
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
			}
			attempt$fit
		},
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
