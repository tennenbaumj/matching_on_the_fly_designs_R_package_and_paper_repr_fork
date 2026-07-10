#' GEE Inference for KK Designs with Count Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using an internal Rcpp
#' solver or \pkg{geepack}) for Poisson (count) responses under a KK 
#' matching-on-the-fly design using the treatment indicator and, optionally,
#' all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountPoissonKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountPoissonKKGEE = R6::R6Class("InferenceCountPoissonKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = utils::modifyList(as.list(InferenceMixinKKGEEShared$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param use_rcpp Whether to use the internal Rcpp solver.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, smart_cold_start_default = NULL){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_gee_shared(des_obj, use_rcpp = use_rcpp, model_formula = model_formula)
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_gee_dispatch(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$compute_kk_gee_jackknife_wald_confidence_interval(alpha = alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$compute_kk_gee_jackknife_wald_two_sided_pval(delta = delta)
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (length(row_weights) > 0L && all(is.finite(row_weights)) &&
			    (max(row_weights) - min(row_weights)) <= sqrt(.Machine$double.eps)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$beta_hat_T = beta_hat_T
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$df = Inf
					private$cached_values$summary_table = NULL
					private$cached_values$nonestimable = FALSE
					private$cached_values$nonestimable_reason = NULL
					private$cached_values$nonestimable_stage = NULL
					return(private$cached_values$beta_hat_T)
				}
			}
			beta_hat_T = private$fit_weighted_gee_with_fallback(row_weights)
			private$cached_values$beta_hat_T = as.numeric(beta_hat_T)[1L]
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = Inf
			private$cached_values$summary_table = NULL
			private$cached_values$nonestimable = !is.finite(private$cached_values$beta_hat_T)
			private$cached_values$nonestimable_reason = if (is.finite(private$cached_values$beta_hat_T)) NULL else "weighted_gee_estimate_unavailable"
			private$cached_values$nonestimable_stage = if (is.finite(private$cached_values$beta_hat_T)) NULL else "estimate"
			private$cached_values$beta_hat_T
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKGEEShared$private), list(
		gee_response_type = function() "count",
		gee_family        = function() stats::poisson(link = "log"),
		shared_gee_dispatch = function(estimate_only = FALSE) private$shared_gee_default(estimate_only)
	))
)
