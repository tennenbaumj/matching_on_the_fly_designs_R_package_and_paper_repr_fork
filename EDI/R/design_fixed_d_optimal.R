#' A D-optimal Search Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed D-optimal experimental design.
#' This design searches for an allocation that maximizes the determinant of the information matrix
#' (equivalent to minimizing the variance of the parameter estimates).
#'
#' @examples
#' des = DesignFixedDOptimal$new(n = 10, response_type = 'continuous')
#' des$add_all_subjects_to_experiment(data.frame(x1 = rnorm(10)))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedDOptimal = R6::R6Class("DesignFixedDOptimal",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a D-optimal search fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T  The probability of the treatment assignment. Must be 0.5 for exchange search.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedDOptimal` object
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
			if (should_run_asserts()) {
				if (prob_T != 0.5){
					stop("D-optimal exchange search currently only supports even treatment allocation (prob_T = 0.5)")
				}
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$uses_covariates = TRUE
		}
	),
	private = list(
		P = NULL, # Projection matrix
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}
			# Precompute projection matrix P = Z0 (Z0' Z0)^-1 Z0'
			if (is.null(private$P)){
				X = private$X[1:n, , drop = FALSE]
				Z0 = cbind(1, X)
				Z0_qr = qr(Z0)
				Q = qr.Q(Z0_qr)
				private$P = Q %*% t(Q)
			}
			# Use C++ speedup
			res = d_optimal_search_cpp(private$P, as.integer(r), as.integer(round(n * private$prob_T)))
			w_mat = res
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	)
)
