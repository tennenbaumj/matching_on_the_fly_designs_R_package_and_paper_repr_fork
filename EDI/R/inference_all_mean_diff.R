#' Mean Difference Inference for Continuous Responses
#'
#' Fits a simple mean difference for continuous responses using the treatment
#' indicator. Note that warm starts are disabled for this class as the simple
#' mean difference is a closed-form estimator and does not benefit from initialization.
#'
#' @examples
#' \dontrun{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
#'
#' seq_des_inf = InferenceAllSimpleMeanDiff$new(seq_des)
#' seq_des_inf$compute_estimate()
#' seq_des_inf$compute_asymp_confidence_interval()
#' seq_des_inf$compute_asymp_two_sided_pval()
#' }
#' @export
InferenceAllSimpleMeanDiff = R6::R6Class("InferenceAllSimpleMeanDiff",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = list(
		#' @description Initialize a simple mean-difference inference object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages. Default \code{FALSE}.
		#' @param max_resample_attempts Maximum number of times a single bootstrap replicate
		#'   may be redrawn when the drawn sample fails validity screening. If all attempts
		#'   fail the replicate is recorded as \code{NA}, silently reducing the effective \code{B}.
		#'   Must be a positive integer. Default \code{50L}.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, max_resample_attempts = 50L, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertCount(max_resample_attempts, positive = TRUE)
			}
			super$initialize(des_obj = des_obj, verbose = verbose, harden = TRUE, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$fit_warm_start_enabled = FALSE
			private$max_resample_attempts = max_resample_attempts
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
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes the appropriate estimate for mean difference
		#'
		#' @return    The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			if (is.null(private$cached_values$beta_hat_T)){
				private$cached_values$yTs = private$y[private$w == 1]
				private$cached_values$yCs = private$y[private$w == 0]
				if (length(private$cached_values$yTs) == 0 || length(private$cached_values$yCs) == 0) {
					return(NA_real_)
				}
				private$cached_values$beta_hat_T = mean(private$cached_values$yTs) - mean(private$cached_values$yCs)
			}
			if (!estimate_only && is.null(private$cached_values$s_beta_hat_T)) {
				private$shared(estimate_only = FALSE)
			}
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			keep = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
			if (!any(keep)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			y_w  = private$y[keep]
			w_w  = private$w[keep]
			rw_w = row_weights[keep]
			rw_T = rw_w[w_w == 1]; y_T = y_w[w_w == 1]
			rw_C = rw_w[w_w == 0]; y_C = y_w[w_w == 0]
			mean_t = sum(y_T * rw_T) / sum(rw_T)
			mean_c = sum(y_C * rw_C) / sum(rw_C)
			private$cached_values$beta_hat_T = mean_t - mean_c
			if (!estimate_only) {
				if (length(y_T) >= 2L && length(y_C) >= 2L) {
					n_eff_T = sum(rw_T)^2 / sum(rw_T^2)
					n_eff_C = sum(rw_C)^2 / sum(rw_C^2)
					s_T_sq  = sum(rw_T * (y_T - mean_t)^2) / sum(rw_T) / (n_eff_T - 1)
					s_C_sq  = sum(rw_C * (y_C - mean_c)^2) / sum(rw_C) / (n_eff_C - 1)
					se = sqrt(s_T_sq + s_C_sq)
					df = (s_T_sq + s_C_sq)^2 / (s_T_sq^2 / (n_eff_T - 1) + s_C_sq^2 / (n_eff_C - 1))
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
					private$cached_values$df = if (is.finite(df) && df > 0) df else NA_real_
				} else {
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$df = NA_real_
				}
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
			}
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$cached_values$beta_hat_T)){
				yTs = private$y[private$w == 1]
				yCs = private$y[private$w == 0]
				if (length(yTs) == 0 || length(yCs) == 0) return(NA_real_)
				private$cached_values$beta_hat_T = mean(yTs) - mean(yCs)
			}
			private$cached_values$beta_hat_T
		},
		max_resample_attempts = 50L,
		get_standard_error = function(){
			if (is.null(private$cached_values$s_beta_hat_T)) private$shared()
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			if (is.null(private$cached_values$df)) private$shared()
			private$cached_values$df
		},
		compute_fast_bootstrap_distr = function(B, ...) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			if (private$is_KK) return(NULL)
			args = list(...)
			n = args[[1]]
			y = args[[2]]
			dead = args[[3]]
			w = args[[4]]
			y_mat = matrix(NA_real_, nrow = n, ncol = B)
			w_mat = matrix(NA_integer_, nrow = n, ncol = B)
			for (b in 1:B) {
				attempt = 1
				repeat {
					i_b = sample(n, n, replace = TRUE)
					w_b = w[i_b]
					if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) {
						if (!private$any_censoring) break
						dead_b_temp = dead[i_b]
						if (any(dead_b_temp[w_b == 1] == 1) && any(dead_b_temp[w_b == 0] == 1) && min(y[i_b]) > 0) break
					}
					attempt = attempt + 1
					if (attempt > private$max_resample_attempts) break
				}
				if (attempt <= private$max_resample_attempts) {
					y_mat[, b] = y[i_b]
					w_mat[, b] = w_b
				}
			}
			res = numeric(B)
			for (b in 1:B) {
				wb = w_mat[, b]
				yb = y_mat[, b]
				res[b] = mean(yb[wb == 1], na.rm = TRUE) - mean(yb[wb == 0], na.rm = TRUE)
			}
			return(res)
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			res = compute_simple_mean_diff_parallel_cpp(as.numeric(y), w_mat, as.numeric(delta), private$n_cpp_threads(ncol(w_mat)))
			return(res)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			
			if (is.null(private$cached_values$beta_hat_T)){
				self$compute_estimate()
			}
			nT = length(private$cached_values$yTs)
			nC = length(private$cached_values$yCs)
			if (nT <= 1 || nC <= 1) { 
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return() 
			}
			s_1_sq = var(private$cached_values$yTs) / nT
			s_2_sq = var(private$cached_values$yCs) / nC
			private$cached_values$s_beta_hat_T = sqrt(s_1_sq + s_2_sq)
			private$cached_values$df = (s_1_sq + s_2_sq)^2 / (s_1_sq^2 / (nT - 1) + s_2_sq^2 / (nC - 1))
			
			private$cached_values$likelihood_test_context = list(
				X = cbind(1, private$w),
				j = 2L,
				full_fit = list(b = c(mean(private$cached_values$yCs), private$cached_values$beta_hat_T), 
								vt = var(private$cached_values$yTs), vc = var(private$cached_values$yCs))
			)
		},
		supports_lik_ratio_param_bootstrap = function() FALSE,
		supports_likelihood_tests = function() FALSE,
		get_supported_testing_types_impl = function(){
			"wald"
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$b)
			vt = spec$full_fit$vt; vc = spec$full_fit$vc
			w = spec$X[, 2]; n = length(w)
			y_sim = numeric(n)
			y_sim[w == 1] = b_null[1] + b_null[2] + rnorm(sum(w == 1), 0, sqrt(vt))
			y_sim[w == 0] = b_null[1] + rnorm(sum(w == 0), 0, sqrt(vc))
			
			list(
				full_fit = list(b = c(mean(y_sim[w==0]), mean(y_sim[w==1]) - mean(y_sim[w==0])), vt = var(y_sim[w==1]), vc = var(y_sim[w==0])),
				fit_null = function(d, start = NULL){
					m_joint = mean(y_sim - w * d)
					list(b = c(m_joint, d), vt = var(y_sim[w==1]), vc = var(y_sim[w==0]))
				},
				neg_loglik = function(fit){
					mu = fit$b[1] + fit$b[2]*w
					sum((y_sim - mu)^2)
				}
			)
		},
		get_likelihood_test_spec = function(){
			NULL
		}
	)
)
