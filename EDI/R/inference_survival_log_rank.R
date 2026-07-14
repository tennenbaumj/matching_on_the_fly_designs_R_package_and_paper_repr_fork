#' Log-Rank Inference for Survival Data with Censoring
#'
#' Non-parametric all-subject inference for survival outcomes supporting right
#' censoring, based on the standard two-sample log-rank test. The treatment effect
#' estimate is the difference in mean martingale residuals between the treatment and
#' control groups under the pooled null hazard. The p-value uses the classic
#' log-rank score statistic with its hypergeometric tie-adjusted variance.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- InferenceSurvivalLogRank$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalLogRank = R6::R6Class("InferenceSurvivalLogRank",
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
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE) {
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		},
		#' @description Computes the treatment-effect estimate on the martingale-residual mean-difference scale.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$compute_shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			private$cached_values$beta_hat_T = private$weighted_logrank_mean_difference(row_weights)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$logrank_score = NA_real_
			private$cached_values$logrank_var = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes a (1 - alpha)-level confidence interval based on the asymptotic
		#' normality of the martingale-residual mean-difference estimate.
		#' Falls back to bootstrap if the estimated standard error is unavailable.
		#' @param alpha Significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes a Wald-style 2-sided p-value by inverting the confidence interval.
		#'
		#' @param delta The null difference to test against. Default is 0.
		#'
		#' @return  The approximate frequentist p-value
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$invert_ci_to_find_two_sided_pval_for_treatment_effect(delta = delta)
		},
		#' @description Computes the standard two-sided log-rank p-value for a zero treatment effect.
		#' @param delta Null treatment effect to test against. Only \code{0} is supported.
		compute_asymp_log_rank_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$compute_shared()
			if (should_run_asserts()) {
				if (delta != 0){
					stop("Testing non-zero delta is not yet implemented for InferenceSurvivalLogRank.")
				}
			}
			if (!is.finite(private$cached_values$logrank_var) || private$cached_values$logrank_var <= 0){
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			chisq_stat = private$cached_values$logrank_score ^ 2 / private$cached_values$logrank_var
			stats::pchisq(chisq_stat, df = 1, lower.tail = FALSE)
		},
		#' @description Randomization confidence intervals are not supported for this class because
		#' the martingale-residual score scale is not commensurate with the transformed
		#' time-ratio null used by the randomization CI algorithm.
		#' @param alpha Unused.
		#' @param r Unused.
		#' @param pval_epsilon Unused.
		#' @param show_progress Unused.
		#' @param ci_search_control Unused.
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			stop("Randomization confidence intervals are not supported for InferenceSurvivalLogRank due to inconsistent estimator units on the log-rank score scale.")
		}
	),
	private = list(
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# survival sharp-null shift is multiplicative (delta on the log scale)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_logrank_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
		weighted_logrank_mean_difference = function(row_weights){
			keep = is.finite(private$y) & is.finite(private$dead) & is.finite(row_weights) & row_weights > 0
			if (!any(keep)) return(NA_real_)
			y = private$y[keep]
			dead = private$dead[keep]
			w = private$w[keep]
			row_weights = as.numeric(row_weights[keep])
			surv_obj = survival::Surv(y, dead)
			cox_null = tryCatch(
				survival::coxph(surv_obj ~ 1, weights = row_weights),
				error = function(e) NULL
			)
			if (is.null(cox_null)) return(NA_real_)
			M = tryCatch(
				as.numeric(stats::residuals(cox_null, type = "martingale")),
				error = function(e) NULL
			)
			if (is.null(M) || length(M) != length(y) || !all(is.finite(M))) return(NA_real_)
			idx_t = w == 1
			idx_c = w == 0
			if (!any(idx_t) || !any(idx_c)) return(NA_real_)
			wt_t = row_weights[idx_t]
			wt_c = row_weights[idx_c]
			if (sum(wt_t) <= 0 || sum(wt_c) <= 0) return(NA_real_)
			mean_t = sum(wt_t * M[idx_t]) / sum(wt_t)
			mean_c = sum(wt_c * M[idx_c]) / sum(wt_c)
			as.numeric(mean_t - mean_c)
		},
		compute_shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			logrank_stats = tryCatch(
				fast_logrank_stats_cpp(
					w = as.integer(private$w),
					y = as.numeric(private$y),
					dead = as.integer(private$dead)
				),
				error = function(e) NULL
			)
			if (is.null(logrank_stats)){
				private$cached_values$beta_hat_T = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$logrank_score = NA_real_
				private$cached_values$logrank_var = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T = as.numeric(logrank_stats$beta_hat)
			private$cached_values$s_beta_hat_T = as.numeric(logrank_stats$se_beta_hat)
			private$cached_values$logrank_score = as.numeric(logrank_stats$score)
			private$cached_values$logrank_var = as.numeric(logrank_stats$var_score)
		}
	)
)
