#' A completely randomized / Bernoulli Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#'
#' @examples
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' @export
DesignSeqOneByOneBernoulli = R6::R6Class("DesignSeqOneByOneBernoulli",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description Initialize a Bernoulli sequential experimental design
		#'
		#' @param  response_type 	The data type of response values which must be one of the following:
		#' 								"continuous",
		#' 								"incidence",
		#' 								"proportion",
		#' 								"count",
		#' 								"survival",
		#' 								"ordinal".
		#' @param  prob_T  The probability of the treatment assignment. This defaults to \code{0.5}.
		#' @param include_is_missing_as_a_new_feature     If missing data is present in a variable,
		#'   should we include another dummy variable for its missingness? The default is \code{TRUE}.
		#' @param  n  		The sample size (if fixed). Default is \code{NULL}.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#' @return  A new `DesignSeqOneByOneBernoulli` object
		#'
		initialize = function(
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,

						verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
		},
		#' @description Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			rbinom(1, 1, private$prob_T)
		}
	),
	private = list(
		draw_ws_raw = function(r = 100){
			generate_permutations_bernoulli_cpp(as.integer(private$t), as.integer(r), as.numeric(private$prob_T))$w_mat
		}
	)
)
