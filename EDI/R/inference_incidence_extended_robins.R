#' Extended Robins Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a block-stratified standard error.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneRandomBlockSize$new(n = 20, response_type = 'incidence', strata_cols = 'x1')
#' for (i in 1:20) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = factor(rep(1:2, 10)[i], levels=1:2)))
#' }
#' seq_des$add_all_subject_responses(rbinom(20, 1, 0.5))
#' inf = InferenceIncidExtendedRobins$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceIncidExtendedRobins = R6::R6Class("InferenceIncidExtendedRobins",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description Computes an approximate confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$get_standard_error()
			super$compute_asymp_confidence_interval(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			private$get_standard_error()
			super$compute_asymp_two_sided_pval(delta)
		},
		#' @description Initialize Extended Robins blocked-design incidence inference.
		#' @param des_obj A completed design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose Logical. Whether to print progress messages.
		#' @return A new \code{InferenceIncidExtendedRobins} object.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (!des_obj$is_blocking_design()) {
				stop("InferenceIncidExtendedRobins requires a blocking design with equal block sizes and even allocation.")
			}
			if (des_obj$get_prob_T() != 0.5) {
				stop("InferenceIncidExtendedRobins requires a blocking design with even allocation.")
			}
			block_ids = des_obj$get_block_ids()
			block_sizes = as.integer(table(block_ids))
			if (length(block_sizes) > 1L && any(block_sizes != block_sizes[1L])) {
				stop("InferenceIncidExtendedRobins requires a blocking design with equal block sizes.")
			}

			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),
	private = list(
		supports_lik_ratio_param_bootstrap = function() FALSE,
		supports_likelihood_tests = function() FALSE,
		get_supported_testing_types_impl = function(){
			"wald"
		},
		get_standard_error = function(){
			if (!is.null(private$cached_values$robins_s_beta_hat_T)) {
				se = private$cached_values$robins_s_beta_hat_T
				if (is.finite(se) && se > 0) return(se)
				private$cache_nonestimable_se("extended_robins_standard_error_unavailable")
				return(NA_real_)
			}
			private$cached_values$robins_s_beta_hat_T = compute_extended_robins_block_se_cpp(
				private$des_obj_priv_int$y,
				private$des_obj$get_w(),
				private$des_obj$get_block_ids(),
				private$des_obj_priv_int$n
			)
			if (!is.finite(private$cached_values$robins_s_beta_hat_T) || private$cached_values$robins_s_beta_hat_T <= 0) {
				private$cached_values$robins_s_beta_hat_T = NA_real_
				private$cache_nonestimable_se("extended_robins_standard_error_unavailable")
				return(NA_real_)
			}
			private$cached_values$s_beta_hat_T = private$cached_values$robins_s_beta_hat_T
			private$cached_values$df = NA_real_
			private$cached_values$robins_s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			NA_real_
		}
	)
)
