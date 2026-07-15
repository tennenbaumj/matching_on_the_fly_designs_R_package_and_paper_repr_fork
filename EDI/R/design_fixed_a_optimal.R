#' An A-optimal Search Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed A-optimal experimental design.
#' This design searches for an allocation that minimizes the trace of the inverse
#' information matrix
#' (equivalent to minimizing the average variance of the parameter estimates).
#'
#' @examples
#' des = DesignFixedAOptimal$new(n = 10, response_type = 'continuous')
#' des$add_all_subjects_to_experiment(data.frame(x1 = rnorm(10)))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedAOptimal = R6::R6Class("DesignFixedAOptimal",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize an A-optimal search fixed experimental design
		#'
		#' @param response_type "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param prob_T Probability of treatment assignment (default 0.5).
		#' @param include_is_missing_as_a_new_feature Flag for missingness indicators.
		#' @param n Sample size (if fixed).
		#' @param verbose Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return A new `DesignFixedAOptimal` object
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
					stop("A-optimal exchange search currently only supports even treatment allocation (prob_T = 0.5)")
				}
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$uses_covariates = TRUE
		}
	),
	private = list(
		P = NULL, # Projection matrix
		H = NULL, # Trace-inverse kernel matrix
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}
			if (is.null(private$P)){
				X = private$X[1:n, , drop = FALSE]
				Z0 = cbind(1, X)

				# P = Z0 (Z0'Z0)^-1 Z0'
				Z0_qr = qr(Z0)
				Q = qr.Q(Z0_qr)
				private$P = Q %*% t(Q)

				# H = Z0 (Z0'Z0)^-2 Z0'
				R = qr.R(Z0_qr)
				M = tryCatch(solve(t(R) %*% R), error = function(e) MASS::ginv(t(R) %*% R))
				H_kernel = M %*% M
				private$H = Z0 %*% H_kernel %*% t(Z0)
			}
			# Use C++ speedup
			res = a_optimal_search_cpp(private$P, private$H, as.integer(r), as.integer(round(n * private$prob_T)))
			w_mat = res
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	)
)
