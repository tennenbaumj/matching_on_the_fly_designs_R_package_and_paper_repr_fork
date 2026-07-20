#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalRestrictedMeanDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalRestrictedMeanDiff = R6::R6Class("InferenceSurvivalRestrictedMeanDiff",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description Initialize the Inference object.
		#'
		#' @param des_obj The design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose If TRUE, print additional information.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL) {
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
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
		#' @description Computes the appropriate estimate for mean difference
		#'
		#' @return  The setting-appropriate (see description) numeric estimate of the treatment effect
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			if (is.null(private$cached_values$beta_hat_T)){
				private$cached_values$beta_hat_T = get_survival_stat_diff(
					private$y,
					private$dead,
					private$w,
					"restricted_mean"
				)
			}
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			private$cached_values$beta_hat_T = private$weighted_survival_stat_diff(
				row_weights,
				requested_stat = "restricted_mean"
			)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes a 1-alpha level frequentist confidence interval
		#' differently for all response types, estimate types, and
		#' test types.
		#'
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal.
		#' Hence these confidence intervals are asymptotically valid
		#' and thus approximate for any sample size.
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#'
		#' @return  A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (is.null(private$cached_values$beta_hat_T)){
				self$compute_estimate()
			}
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$compute_s_beta_hat_T()
			}
			if (is.na(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0) {
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes a 2-sided p-value via the log rank test
		#'
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		#'
		#' @return  The approximate frequentist p-value
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (delta == 0){
				if (is.null(private$cached_values$s_beta_hat_T)){
					private$compute_s_beta_hat_T()
				}
				if (is.na(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0) {
					return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
				}
				z_beta_hat_T = private$cached_values$beta_hat_T / private$cached_values$s_beta_hat_T
				2 * min(stats::pnorm(z_beta_hat_T), 1 - stats::pnorm(z_beta_hat_T))
			} else {
				if (should_run_asserts()) {
					stop("TO-DO")
				}
				NA_real_
			}
		},
		#' @description Computes a 1-alpha level frequentist confidence interval for the randomization test
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#' @param  r  	The number of randomization vectors. The default is 501.
		#' @param  pval_epsilon  		The bisection algorithm tolerance. The default is 0.005.
		#' @param  show_progress  	Show a text progress indicator.
		#' @param ci_search_control Unused.
		#' @return  A 1 - alpha sized frequentist confidence interval
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			stop("Randomization confidence intervals are not supported for InferenceSurvivalRestrictedMeanDiff due to inconsistent estimator units on the transformed scale (estimates time difference, but randomization test searches for log-time ratio).")
		}
	),
	private = list(
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), TRUE, mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
		weighted_survival_stat_for_group = function(y, dead, row_weights, requested_stat = c("median", "restricted_mean")){
			requested_stat = match.arg(requested_stat)
			keep = is.finite(y) & is.finite(dead) & is.finite(row_weights) & row_weights > 0
			if (!any(keep)) return(NA_real_)
			y = y[keep]
			dead = dead[keep]
			row_weights = as.numeric(row_weights[keep])
			fit = tryCatch(
				survival::survfit(
					survival::Surv(y, dead) ~ 1,
					weights = row_weights
				),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NA_real_)
			if (requested_stat == "median") {
				q = tryCatch(stats::quantile(fit, probs = 0.5), error = function(e) NULL)
				med = if (!is.null(q)) as.numeric(q$quantile) else NA_real_
				return(if (length(med)) med[1L] else NA_real_)
			}
			tau = max(y)
			times = c(0, fit$time)
			surv_vals = c(1, fit$surv)
			if (!length(times) || !length(surv_vals)) return(NA_real_)
			area = 0
			for (i in seq_len(length(times) - 1L)) {
				area = area + surv_vals[i] * (times[i + 1L] - times[i])
			}
			if (length(times) >= 1L) {
				area = area + surv_vals[length(surv_vals)] * (tau - times[length(times)])
			}
			as.numeric(area)
		},
		weighted_survival_stat_diff = function(row_weights, requested_stat = c("median", "restricted_mean")){
			requested_stat = match.arg(requested_stat)
			idx_t = private$w == 1
			idx_c = private$w == 0
			if (!any(idx_t) || !any(idx_c)) return(NA_real_)
			stat_t = private$weighted_survival_stat_for_group(
				private$y[idx_t], private$dead[idx_t], row_weights[idx_t], requested_stat = requested_stat
			)
			stat_c = private$weighted_survival_stat_for_group(
				private$y[idx_c], private$dead[idx_c], row_weights[idx_c], requested_stat = requested_stat
			)
			if (!is.finite(stat_t) || !is.finite(stat_c)) return(NA_real_)
			as.numeric(stat_t - stat_c)
		},
		compute_s_beta_hat_T = function(){
			se_val = get_restricted_mean_se_diff(
				private$y,
				private$dead,
				private$w
			)
			if (is.na(se_val) || se_val <= 0) {
				warning("Restricted mean SE is non-positive or NA; MLE p-value/CI unavailable.")
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$s_beta_hat_T = se_val
		}
	)
)
