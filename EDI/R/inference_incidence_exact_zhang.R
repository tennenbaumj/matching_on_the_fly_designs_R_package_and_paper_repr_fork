#' Exact Zhang Incidence Inference
#'
#' @examples
#' \dontrun{
#' # Example for InferenceIncidenceExactZhang
#' }
#' @export
InferenceIncidenceExactZhang = R6::R6Class("InferenceIncidenceExactZhang",
	lock_objects = FALSE,
	inherit = InferenceExact,
	public = list(
		#' @description Initialize exact Zhang incidence inference.
		#' @param des_obj A completed design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		#' @return A new \code{InferenceIncidenceExactZhang} object.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Compute the Zhang incidence treatment estimate.
		#' @param estimate_only Ignored for this estimator.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			stats = zhang_get_exact_stats(self)
			zhang_incid_treatment_estimate(stats)
		},
		#' @description Compute an exact Zhang confidence interval.
		#' @param alpha Significance level.
		#' @param pval_epsilon Bisection tolerance for the inversion routine.
		#' @param type Exact inference type.
		#' @param args_for_type Optional arguments keyed by exact type.
		#' @return A confidence interval.
		compute_exact_confidence_interval = function(alpha = 0.05, pval_epsilon = 0.005, type = NULL, args_for_type = NULL){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			}
			exact_type = private$resolve_exact_type(type)
			exact_args = private$normalize_exact_inference_args(
				exact_type,
				args_for_type = args_for_type,
				pval_epsilon = pval_epsilon
			)
			private$compute_exact_confidence_interval_by_type(exact_type, alpha, exact_args)
		}
	),
	private = list(
		default_exact_type = "Zhang",
		supports_bayesian_bootstrap = function() FALSE,
		resolve_exact_type = function(type){
			if (is.null(type)) type = private$default_exact_type
			if (should_run_asserts()) {
				assertChoice(type, c("Zhang"))
			}
			type
		},
		normalize_exact_inference_args = function(type, args_for_type = NULL, pval_epsilon = NULL){
			zhang_normalize_exact_inference_args(type, args_for_type = args_for_type, pval_epsilon = pval_epsilon)
		},
		assert_exact_inference_params = function(type, args_for_type){
			zhang_assert_exact_inference_params(self, type, args_for_type)
		},
		compute_exact_confidence_interval_by_type = function(type, alpha, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Zhang = zhang_ci_exact_combined(
					self,
					alpha = alpha,
					pval_epsilon = args_for_type[[type]]$pval_epsilon,
					combination_method = args_for_type[[type]]$combination_method
				)
			)
		},
		compute_exact_two_sided_pval_for_treatment_effect_by_type = function(type, delta, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Zhang = zhang_pval_exact_combined(
					self,
					delta_0 = delta,
					combination_method = args_for_type[[type]]$combination_method
				)
			)
		}
	)
)

# Backwards-compatible alias for the shortened legacy name.
InferenceIncidExactZhang = InferenceIncidenceExactZhang
