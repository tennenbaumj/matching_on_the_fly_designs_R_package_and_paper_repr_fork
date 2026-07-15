#' Wei's Urn Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for Wei's Urn sequential
#' experimental design.
#' This design uses an adaptive urn model where the probability of assignment to a group
#' decreases as the number of subjects in that group increases.
#'
#' @examples
#' seq_des = DesignSeqOneByOneUrn$new(n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' @export
DesignSeqOneByOneUrn = R6::R6Class("DesignSeqOneByOneUrn",
	inherit = DesignSeqOneByOne,
	public = list(
		#'
		#' @description Initialize an Urn sequential experimental design
		#'
		#' @param alpha The initial number of balls of each type (Treatment/Control) in the urn.
		#' @param beta The number of balls of the opposite type to add to the urn after an assignment.
		#' @param  response_type 	The data type of response values.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param  n  		The sample size.
		#' @param verbose A flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#' @return  A new `DesignSeqOneByOneUrn` object
		#'
		initialize = function(
						alpha = 1,
						beta = 1,
						response_type,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,

						verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
			if (should_run_asserts()) {
				assertNumber(alpha, lower = 0)
				assertNumber(beta, lower = 0)
			}

			super$initialize(response_type, 0.5, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			
			private$alpha = alpha
			private$beta = beta
		},
		#' @description Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			# Probability of Treatment based on current counts
			# In Wei's Urn, P(T) = (alpha + beta * nC) / (2 * alpha + beta * (nT + nC))
			# where nT and nC are current counts of assigned subjects.
			nT = sum(private$w == 1, na.rm = TRUE)
			nC = sum(private$w == 0, na.rm = TRUE)
			
			prob_T = (private$alpha + private$beta * nC) / (2 * private$alpha + private$beta * (nT + nC))
			
			rbinom(1, 1, prob_T)
		}
	),
	private = list(
		alpha = NULL,
		beta = NULL
	)
)
