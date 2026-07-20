#' Rank Regression Inference for Survival Responses under KK Designs
#'
#' Fits a multivariate Gehan-Wilcoxon rank regression for survival outcomes under
#' a KK matching-on-the-fly design. The model adjusts for the treatment indicator
#' and, optionally, all recorded covariates.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKRankRegrIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferenceSurvivalKKRankRegrIVWC = R6::R6Class("InferenceSurvivalKKRankRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKSurvivalRankRegrIVWC,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj = des_obj, model_formula = model_formula, verbose = verbose)
		}
	),
	private = list(
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		}
	)
)
