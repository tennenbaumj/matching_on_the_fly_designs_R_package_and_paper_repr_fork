#' Inference based on Maximum Likelihood for KK designs
#'
#' Inference for mean difference. Note that warm starts are disabled for this class
#' as the Bai adjusted t-test is a closed-form estimator and does not benefit from initialization.
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceBaiAdjustedTKK14$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceBaiAdjustedTKK14 = R6::R6Class("InferenceBaiAdjustedTKK14",
	lock_objects = FALSE,
	inherit = InferenceBaiAdjustedT,
	public = list(

	),

	private = list(
	distance = function(avg1, avg2){
		sum((avg1 - avg2)^2)
	}
	)
)
