#' An incomplete / balanced completely randomized Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments.
#'
#' @examples
#' seq_des = DesignSeqOneByOneiBCRD$new(n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' @export
DesignSeqOneByOneiBCRD = R6::R6Class("DesignSeqOneByOneiBCRD",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description Initialize a balanced sequential experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param  prob_T  Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param  n  		The sample size.
		#' @param verbose A flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return  A new `DesignSeqOneByOneiBCRD` object
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
		#' @description Adds a subject and assigns treatment.
		#' Sets a single-block match vector once all subjects have arrived.
		#' @param x_new A data frame with one row representing the new subject's covariates.
		#' @return The treatment assignment (0 or 1).
		add_one_subject_to_experiment_and_assign = function(x_new){
			w_t = super$add_one_subject_to_experiment_and_assign(x_new)
			private$m = rep(1L, private$t)
			w_t
		},
		#' @description Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			nT = sum(private$w == 1, na.rm = TRUE)
			nC = sum(private$w == 0, na.rm = TRUE)
			
			if (is.null(private$n)){
				#if n is not fixed, we cannot really ensure balance at the end, 
				#so we just use Bernoulli
				private$assign_wt_Bernoulli()
			} else {
				#if n is fixed, we use the remaining slots
				nT_rem = round(private$n * private$prob_T) - nT
				nC_rem = (private$n - round(private$n * private$prob_T)) - nC
				
				if (nT_rem <= 0) return(0)
				if (nC_rem <= 0) return(1)
				
				rbinom(1, 1, nT_rem / (nT_rem + nC_rem))
			}
		}
	),
	private = list(
		draw_ws_raw = function(r = 100){
			generate_permutations_ibcrd_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		}
	)
)
