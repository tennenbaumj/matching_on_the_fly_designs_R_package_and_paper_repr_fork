#' A Greedy Search Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed greedy experimental design.
#' Uses a native C++ (RcppEigen) greedy swap search for balanced allocations.
#'
#' @examples
#' \dontrun{
#' des = DesignFixedGreedy$new(n = 10, response_type = 'continuous')
#' }
#' @export
DesignFixedGreedy = R6::R6Class("DesignFixedGreedy",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a greedy search fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T  The probability of the treatment assignment. Must be 0.5.
		#' @param objective 	The objective function to use. Default is "mahal_dist".
		#' @param n_iter Number of swap iterations. \code{Inf} (default) uses exhaustive
		#'   best-improvement search guaranteed to reach a strict local optimum. A positive
		#'   integer runs that many stochastic random-pair iterations with patience-based
		#'   early stopping.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedGreedy` object
		#'
		initialize = function(
				response_type,
				prob_T = 0.5,
				objective = "mahal_dist",
				n_iter = Inf,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
			if (should_run_asserts()) {
				if (prob_T != 0.5){
					stop("Greedy designs currently only support even treatment allocation (prob_T = 0.5)")
				}
			}
			# GED availability is checked lazily in draw_ws_according_to_design so
			# that workers using pre-computed w vectors never trigger a JVM load.
			if (should_run_asserts()) {
				if (!is.infinite(n_iter) && (!is.numeric(n_iter) || length(n_iter) != 1L || n_iter <= 0 || n_iter != floor(n_iter)))
					stop("n_iter must be Inf or a positive integer")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$objective = objective
			private$n_iter    = n_iter
			private$uses_covariates = TRUE
		},
		#' @description Whether this design supports batch pregeneration of treatment vectors.
		#'
		#' @return \code{TRUE}.
		supports_batch_w_pregeneration = function() TRUE
	),
	private = list(
		objective = NULL,
		n_iter    = NULL,
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
			}
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			if (should_run_asserts()) {
				if (n %% 2 != 0){
					stop("DesignFixedGreedy requires an even number of subjects.")
				}
			}
			if (is.null(private$X) || ncol(private$X) == 0){
				w_mat = replicate(r, sample(c(rep(1, n / 2), rep(0, n / 2))))
				return(private$validate_allocation_matrix(w_mat, n = n, r = r))
			}
			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X          = private$X[1:n, , drop = FALSE]
			cpp_n_iter = if (is.infinite(private$n_iter)) -1L else as.integer(private$n_iter)
			w_mat      = greedy_design_search_cpp(
				X_raw     = X,
				r         = as.integer(r),
				objective = private$objective,
				n_iter    = cpp_n_iter
			)
			# greedy_design_search_cpp already returns n x r
			storage.mode(w_mat) = "numeric"
			private$validate_allocation_matrix(w_mat, n = n, r = r)
		},
		validate_allocation_matrix = function(w_mat, n, r){
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (should_run_asserts()) {
				if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < r) {
					stop("DesignFixedGreedy returned an unexpected allocation matrix shape.")
				}
			}
			w_mat = w_mat[, seq_len(r), drop = FALSE]
			storage.mode(w_mat) = "numeric"
			if (should_run_asserts()) {
				if (any(!is.finite(w_mat)) || any(is.na(w_mat))) {
					stop("DesignFixedGreedy returned non-finite treatment assignments.")
				}
				if (any(!(w_mat %in% c(0, 1)))) {
					stop("DesignFixedGreedy returned an invalid treatment assignment matrix.")
				}
			}
			treated_counts = colSums(w_mat)
			if (should_run_asserts()) {
				if (any(treated_counts != n / 2)) {
					stop("DesignFixedGreedy returned an unbalanced allocation.")
				}
			}
			w_mat
		}
	)
)
