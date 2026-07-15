#' A Factorial Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed factorial experimental design.
#' This design handles multiple treatment factors and balances assignments across
#' all factor combinations.
#'
#' @examples
#' des = DesignFixedFactorial$new(n = 12, response_type = 'continuous', factors = list(f1 = 2, f2 = 2))
#' des$add_all_subjects_to_experiment(data.frame(x=1:12))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedFactorial = R6::R6Class("DesignFixedFactorial",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a factorial fixed experimental design
		#'
		#' @param factors         A list where names are factor names and values are number of
		#'   levels (e.g. list(A=2, B=2)).
		#' @param response_type 	The data type of response values.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedFactorial` object
		initialize = function(
				factors,
				response_type,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,

				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
			if (should_run_asserts()) {
				assertList(factors, types = "numeric", min.len = 1)
			}
			# We don't use prob_T in the standard way here, as we have multiple factors
			# But base Design needs it. We'll set it to 0.5.
			super$initialize(response_type, 0.5, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$factors = factors
			
			# Precompute all combinations
			private$combinations = expand.grid(lapply(factors, function(l) 1:l))
			private$num_combinations = nrow(private$combinations)
		},
		#' @description Assign treatment (factor combination) to all subjects.
		#' For factorial designs the w vector stores factor combination indices.
		#' @param w_precomputed Optional pre-computed factor combination index vector.
		assign_w_to_all_subjects = function(w_precomputed = NULL){
			if (!is.null(w_precomputed)) {
				private$w[1:self$get_n()] = as.integer(w_precomputed)
			} else {
				private$w[1:self$get_n()] = private$draw_ws_raw(1)[, 1]
			}
		},
		#' @description Draw multiple treatment assignment vectors according to balanced factorial randomization.
		#' NOTE: Factorial designs return factor combination indices (1..num_combinations), not {-1,+1}.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r with factor combination indices.
		draw_ws_according_to_design = function(r = 100){
			private$draw_ws_raw(r)
		},
		#' @description Get the factor combination index assigned to each subject.
		#' NOTE: For factorial designs this returns raw factor indices, not {-1,+1}.
		#'
		#' @return 		An integer vector of factor combination indices.
		get_w = function(){
			private$w
		},
		#' @description Get the data frame of factor assignments for each subject.
		#'
		#' @return A data frame with n rows and columns corresponding to factors.
		get_w_factorial = function(){
			w_idx = private$w
			if (length(w_idx) == 0 || any(is.na(w_idx))) return(NULL)
			private$combinations[w_idx, , drop = FALSE]
		}
	),
	private = list(
		factors = NULL,
		combinations = NULL,
		num_combinations = NULL,
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			# Standard balanced factorial: each combination appears n / num_combinations times
			base_alloc = rep(1:private$num_combinations, length.out = n)
			w_mat = matrix(NA_integer_, nrow = n, ncol = r)
			for (j in 1:r){
				w_mat[, j] = sample(base_alloc)
			}
			w_mat
		}
	)
)
