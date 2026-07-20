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
#' inf = InferenceSurvivalKMDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKMDiff = R6::R6Class("InferenceSurvivalKMDiff",
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
		#'
		#' @examples
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = InferenceSurvivalKMDiff$new(seq_des)
		#' seq_des_inf$compute_estimate()
		#'
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			private$cached_values$beta_hat_T = private$weighted_survival_stat_diff(row_weights, requested_stat = "median")
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes a (1 - alpha)-level confidence interval for the difference in Kaplan-Meier
		#' median survival times (treatment minus control).
		#'
		#' The Brookmeyer-Crowley confidence interval is obtained for each group's median
		#' separately via \code{survival::survfit} (using a log-log transformation of the
		#' survival function by default). The per-group SE is back-calculated from the CI
		#' half-width as \eqn{\hat\sigma_i = (\text{upper}_i - \text{lower}_i) / (2 z_{\alpha/2})}.
		#' The two groups are independent by design, so the SE of the difference is
		#' \eqn{\sqrt{\hat\sigma_T^2 + \hat\sigma_C^2}}, and the CI is
		#' \eqn{(\hat{m}_T - \hat{m}_C) \pm z_{\alpha/2} \cdot \sqrt{\hat\sigma_T^2 +
		#' \hat\sigma_C^2}}.
		#'
		#' Falls back to \code{compute_bootstrap_confidence_interval} when either group's
		#' median is not estimable (i.e., the Kaplan-Meier curve does not reach 0.5) or
		#' when the Brookmeyer-Crowley CI bounds are \code{NA}.
		#'
		#' @param alpha           The significance level; the confidence level is 1 - \code{alpha}.
		#'   Default is 0.05.
		#'
		#' @return  A numeric vector of length 2 giving the (lower, upper) confidence bounds
		#' 			for the difference in median survival times, on the original time scale.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (is.na(private$cached_values$s_beta_hat_T)) {
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes a Wald-style 2-sided p-value based on the median difference
		#' and its back-calculated standard error.
		#'
		#' @param delta The null difference to test against. Default is 0.
		#'
		#' @return  The approximate frequentist p-value
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (is.na(private$cached_values$s_beta_hat_T)) {
				return(self$compute_bootstrap_two_sided_pval(delta = delta))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes a 2-sided p-value via the log rank test
		#'
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		#'
		#' @return  The approximate frequentist p-value
		compute_asymp_log_rank_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
				if (delta != 0){
					stop("Log-rank p-value for non-zero delta is not yet implemented.")
				}
			}
			survival_obj = survival::Surv(private$y, private$dead)
			surv_diff = survival::survdiff(survival_obj ~ private$w)
			surv_diff$pvalue
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
			stop("Randomization confidence intervals are not supported for InferenceSurvivalKMDiff due to inconsistent estimator units on the transformed scale.")
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
				as.numeric(delta), FALSE, mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			
			# Point estimate via fast Rcpp
			private$cached_values$beta_hat_T = get_survival_stat_diff(
				private$y,
				private$dead,
				private$w,
				"median"
			)
			if (estimate_only) return(invisible(NULL))
			
			# Variance components via survfit (Brookmeyer-Crowley back-calculation)
			y    = private$y
			dead = private$dead
			w    = private$w
			alpha = 0.05 # standard anchor for SE back-calculation
			
			# Use a single grouped survfit call to match canonical R performance path
			fit = tryCatch(survival::survfit(survival::Surv(y, dead) ~ w, conf.int = 1 - alpha), error = function(e) NULL)
			
			if (is.null(fit)){
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}
			
			q = stats::quantile(fit, 0.5)
			strata_names = rownames(q$lower)
			
			# Extract bounds for treatment (w=1) and control (w=0)
			idx_T = which(strata_names == "w=1")
			idx_C = which(strata_names == "w=0")
			
			if (length(idx_T) != 1L || length(idx_C) != 1L) {
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}
			
			lo_T = as.numeric(q$lower[idx_T, 1])
			hi_T = as.numeric(q$upper[idx_T, 1])
			lo_C = as.numeric(q$lower[idx_C, 1])
			hi_C = as.numeric(q$upper[idx_C, 1])
			
			if (!is.finite(lo_T) || !is.finite(hi_T) || !is.finite(lo_C) || !is.finite(hi_C)){
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}
			
			z = stats::qnorm(1 - alpha / 2)
			private$cached_values$s_beta_hat_T = sqrt(((hi_T - lo_T) / (2 * z))^2 + ((hi_C - lo_C) / (2 * z))^2)
		},
		weighted_survival_stat_for_group = function(y, dead, row_weights, requested_stat = c("median", "restricted_mean")){
			requested_stat = match.arg(requested_stat)
			keep = is.finite(y) & is.finite(dead) & is.finite(row_weights) & row_weights > 0
			if (!any(keep)) return(NA_real_)
			y = y[keep]
			dead = dead[keep]
			row_weights = as.numeric(row_weights[keep])
			if (requested_stat == "median") {
				return(private$weighted_km_median(y, dead, row_weights))
			}
			fit = tryCatch(
				survival::survfit(
					survival::Surv(y, dead) ~ 1,
					weights = row_weights
				),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NA_real_)
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
		weighted_km_median = function(y, dead, row_weights){
			ord = order(y)
			y = as.numeric(y[ord])
			dead = as.integer(dead[ord])
			row_weights = as.numeric(row_weights[ord])
			survival_prob = 1.0
			unique_times = 0
			survival_probs = 1
			i = 1L
			n = length(y)
			while (i <= n) {
				current_time = y[i]
				j = i
				event_weight = 0
				at_risk_weight = sum(row_weights[i:n])
				while (j <= n && y[j] == current_time) {
					if (dead[j] == 1L) {
						event_weight = event_weight + row_weights[j]
					}
					j = j + 1L
				}
				if (event_weight > 0 && at_risk_weight > 0) {
					survival_prob = survival_prob * (1 - event_weight / at_risk_weight)
					unique_times = c(unique_times, current_time)
					survival_probs = c(survival_probs, survival_prob)
				}
				i = j
			}
			for (k in seq_along(survival_probs)) {
				if (survival_probs[k] < 0.5) {
					if (k > 1L) {
						p1 = survival_probs[k - 1L]
						p2 = survival_probs[k]
						t1 = unique_times[k - 1L]
						t2 = unique_times[k]
						return(as.numeric(t1 + (t2 - t1) * (0.5 - p1) / (p2 - p1)))
					}
					return(as.numeric(unique_times[k]))
				}
			}
			Inf
		}
	)
)
