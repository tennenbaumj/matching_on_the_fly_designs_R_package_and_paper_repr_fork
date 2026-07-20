#' Conditional Logistic Plus GLMM IVWC Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKCondLogitPlusGLMMOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferenceIncidKKCondLogitPlusGLMMIVWC = R6::R6Class("InferenceIncidKKCondLogitPlusGLMMIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKCondLogitPlusGLMM,
	public = list(
	),
	private = list(
		combine_reservoir_into_glmm = function() FALSE
	)
)

#' Conditional Logistic Plus GLMM Combined-Likelihood Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKCondLogitPlusGLMMOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKCondLogitPlusGLMMOneLik = R6::R6Class("InferenceIncidKKCondLogitPlusGLMMOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKCondLogitPlusGLMM,
	public = list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param max_abs_reasonable_coef Cap for reasonable coefficient estimates.
		#' @param max_abs_reasonable_se Cap for reasonable treatment standard errors.
		#' @param max_abs_log_sigma Cap for reasonable log random effect variance.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart optimizer start values.
		#' @param optimization_alg Character. Optimization algorithm (default "lbfgs").
		initialize = function(des_obj, model_formula = NULL, max_abs_reasonable_coef = 50, max_abs_reasonable_se = 1.25, max_abs_log_sigma = 8, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			super$initialize(des_obj, model_formula = model_formula, max_abs_reasonable_coef = max_abs_reasonable_coef, max_abs_reasonable_se = max_abs_reasonable_se, max_abs_log_sigma = max_abs_log_sigma, verbose = verbose, smart_cold_start_default = smart_cold_start_default, optimization_alg = optimization_alg)
		}
	),
	private = list(
		combine_reservoir_into_glmm = function() TRUE
	)
)
