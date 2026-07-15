#' Newcombe Risk-Difference Inference for Binary Responses
#'
#' Fits the Newcombe hybrid score method (Method 10) for the risk difference in a
#' two-arm binary trial. This method constructs a confidence interval for the
#' difference between two independent proportions by combining Wilson score intervals
#' for each group.
#'
#' @details
#' This class is unadjusted and assumes independent samples (e.g. from a Bernoulli). It
#' ignores any matched-pair structure if present. For matched data, use
#' \code{InferenceIncidKKNewcombeRiskDiff}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidNewcombeRiskDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidNewcombeRiskDiff = R6::R6Class("InferenceIncidNewcombeRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description Initialize a Newcombe risk-difference inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes the observed risk-difference estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			i_t = private$w == 1
			i_c = private$w == 0
			w_t = sum(row_weights[i_t])
			w_c = sum(row_weights[i_c])
			if (!is.finite(w_t) || !is.finite(w_c) || w_t <= 0 || w_c <= 0) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			p_t = sum(row_weights[i_t] * as.numeric(private$y[i_t])) / w_t
			p_c = sum(row_weights[i_c] * as.numeric(private$y[i_c])) / w_c
			private$cached_values$counts = list(
				n_t = w_t, n_c = w_c,
				x_t = sum(row_weights[i_t] * as.numeric(private$y[i_t])),
				x_c = sum(row_weights[i_c] * as.numeric(private$y[i_c])),
				p_t = p_t, p_c = p_c
			)
			private$cached_values$beta_hat_T = p_t - p_c
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes a 1 - \code{alpha} Newcombe confidence interval.
		#' @param alpha The significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			counts = private$cached_values$counts
			if (is.null(counts)) return(c(NA_real_, NA_real_))
			
			ci = newcombe_independent_ci_cpp(counts$x_t, counts$n_t, counts$x_c, counts$n_c, alpha)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		#' @description Computes a two-sided p-value by inverting the Newcombe interval.
		#' @param delta The null risk difference.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
			}
			private$shared()
			counts = private$cached_values$counts
			if (is.null(counts) || counts$n_t == 0 || counts$n_c == 0) return(NA_real_)
			
			# Invert the Newcombe CI to find the largest alpha such that delta is on the boundary
			p_fn = function(a) {
				ci = newcombe_independent_ci_cpp(counts$x_t, counts$n_t, counts$x_c, counts$n_c, a)
				if (delta < private$cached_values$beta_hat_T) ci[1] - delta else ci[2] - delta
			}
			
			res = tryCatch(stats::uniroot(p_fn, interval = c(1e-10, 1 - 1e-10))$root, error = function(e) NA_real_)
			if (!is.finite(res)) return(1.0) # If root not found, delta likely far inside
			res
		}
	),
	private = list(
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		get_counts = function(){
			i_t = private$w == 1
			i_c = private$w == 0
			n_t = sum(i_t)
			n_c = sum(i_c)
			x_t = sum(private$y[i_t])
			x_c = sum(private$y[i_c])
			list(
				n_t = n_t, n_c = n_c,
				x_t = x_t, x_c = x_c,
				p_t = if (n_t > 0) x_t / n_t else NA_real_,
				p_c = if (n_c > 0) x_c / n_c else NA_real_
			)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			counts = private$get_counts()
			private$cached_values$counts = counts
			private$cached_values$beta_hat_T = counts$p_t - counts$p_c
			if (estimate_only) return(invisible(NULL))
		},
		get_supported_testing_types_impl = function(){
			character(0)
		}
	)
)
