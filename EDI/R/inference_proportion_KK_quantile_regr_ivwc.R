#' Quantile Regression Compound Estimator for KK Matching-on-the-Fly Designs (Proportion
#' Outcomes)
#'
#' A variance-weighted compound quantile regression estimator for KK matching-on-the-fly
#' designs with proportion responses. Inference is performed on the \strong{logit (log-odds)
#' scale}: responses \eqn{y \in (0,1)} are transformed via \eqn{\text{logit}(y) = \log(y/(1-y))}
#' before quantile regression.
#'
#' The estimator combines:
#' \enumerate{
#'   \item Quantile regression on logit-scale within-pair differences
#'     \eqn{\text{logit}(y_T) - \text{logit}(y_C)} (matched pairs)
#'   \item Quantile regression of \eqn{\text{logit}(y)} on treatment and covariates (reservoir)
#' }
#' using the same variance-weighted combination logic as the OLS compound estimator.
#'
#' The estimated treatment effect is a \strong{log-odds-ratio shift} at quantile \code{tau}.
#' At \code{beta_T = 1} (one log-odds-ratio unit of treatment effect), the population
#' treatment effect on the logit scale is exactly 1, so no \code{skip_ci} is needed.
#'
#' \strong{Default quantile: \code{tau = 0.5} (median regression).}
#' To target a different quantile --- for example the 25th or 75th percentile --- pass
#' \code{tau = 0.25} or \code{tau = 0.75} to the constructor:
#' \preformatted{
#'   inf = InferencePropKKQuantileRegrIVWC$
#'   new(seq_des, tau = 0.75)
#' }
#' Any value strictly between 0 and 1 is accepted.
#'
#' Standard errors use Powell's "nid" sandwich estimator (non-iid), falling back to "iid"
#' on failure. Asymptotic z-based inference is used throughout.
#'
#' This class requires the \pkg{quantreg} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{quantreg} before using this class.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropKKQuantileRegrIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferencePropKKQuantileRegrIVWC = R6::R6Class("InferencePropKKQuantileRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKQuantileRegrIVWC,
	public = list(
		#' @description Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param tau                             The quantile level for regression on the logit
		#'   scale, strictly between 0 and 1.
		#' 							The default \code{tau = 0.5} estimates the median log-odds-ratio treatment effect.
		#' 							Pass a different value (e.g. \code{tau = 0.25} or \code{tau = 0.75}) to target a
		#' 							different percentile of the treatment effect distribution.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{FALSE}.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, tau = 0.5,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
			}
			super$initialize(des_obj = des_obj, model_formula = model_formula, tau = tau, transform_y_fn = qlogis, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$y = .sanitize_proportion_response(private$y, interior = TRUE)
			if (should_run_asserts()) {
				assertNumeric(private$y, any.missing = FALSE, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			}
			# Rebuild KK summary data after sanitizing the response, otherwise the
			# superclass cache can retain raw 0/1 values and qlogis() will produce
			# non-finite values during the quantile fit and its randomization refits.
			private$cached_values$KKstats = NULL
			private$compute_basic_match_data()
		}
	)
)
