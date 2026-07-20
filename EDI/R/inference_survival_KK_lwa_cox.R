#' LWA-style Marginal Cox IVWC Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using a marginal Cox model with Lee-Wei-Amato style cluster-robust variance for
#' matched pairs and standard Cox regression for reservoir subjects.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKLWACoxPHIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferenceSurvivalKKLWACoxPHIVWC = R6::R6Class("InferenceSurvivalKKLWACoxPHIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKLWACoxIVWC,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)
#' LWA-style Marginal Cox Combined-Likelihood Inference for KK Designs
#'
#' Fits a combined-likelihood Cox model for KK matching-on-the-fly designs with
#' survival responses using a marginal approach over all subjects.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKLWACoxPHOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKLWACoxPHOneLik = R6::R6Class("InferenceSurvivalKKLWACoxPHOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKLWACoxOneLik,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)
