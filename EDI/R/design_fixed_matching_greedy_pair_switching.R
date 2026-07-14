#' A Fixed Design Combining Binary Matching and Greedy Pair Switching
#'
#' An R6 class encapsulating a fixed experimental design that first computes
#' binary matches and then improves the matched design with greedy pair switching
#' using a native C++ implementation.
#'
#' @examples
#' \dontrun{
#' des = DesignFixedMatchingGreedyPairSwitching$new(n = 10, response_type = 'continuous')
#' }
#' @export
DesignFixedMatchingGreedyPairSwitching = R6::R6Class("DesignFixedMatchingGreedyPairSwitching",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a fixed design that performs binary matching followed by greedy pair switching.
		#'
		#' @param response_type The data type of response values.
		#' @param prob_T The probability of treatment assignment. Must be \code{0.5}.
		#' @param include_is_missing_as_a_new_feature Flag for missingness indicators.
		#' @param n The sample size.
		#' @param verbose A flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param model_formula A formula object.
		#' @param objective The imbalance objective. Either \code{"mahal_dist"} (default) or \code{"abs_sum_diff"}.
		#' @param n_iter Number of swap iterations. \code{Inf} (default) uses exhaustive
		#'   best-improvement search guaranteed to reach a strict local optimum. A positive
		#'   integer runs that many stochastic random-pair iterations with patience-based
		#'   early stopping.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return A new \code{DesignFixedMatchingGreedyPairSwitching} object.
		initialize = function(
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n,
				verbose = FALSE,
				objective = "mahal_dist",
				n_iter = Inf,
				missingness_method = "impute",
				model_formula = ~ .,
				seed = NULL
			) {
			if (should_run_asserts()) {
				if (prob_T != 0.5) {
					stop("DesignFixedMatchingGreedyPairSwitching only supports balanced designs (prob_T = 0.5).")
				}
			}
			if (should_run_asserts()) {
				if (!is.infinite(n_iter) && (!is.numeric(n_iter) || length(n_iter) != 1L || n_iter <= 0 || n_iter != floor(n_iter)))
					stop("n_iter must be Inf or a positive integer")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, model_formula, seed = seed)
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
		bms = NULL,
		ensure_pair_structure_computed = function(){
			if (is.null(private$bms)) {
				n = self$get_n()
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
				if (is.null(private$X) || ncol(private$X) == 0) return(invisible(NULL))
				X = private$X[1:n, , drop = FALSE]
				private$bms = compute_binary_match_structure(X, mahal_match = (private$objective == "mahal_dist"))
			}
			invisible(NULL)
		},
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			# The greedy search only flips assignments within binary-match pairs
			# (pair_cur_t in greedy_design_search_cpp), so w always has exactly one
			# treated subject per pair: the exchangeable resampling unit is the pair.
			private$ensure_pair_structure_computed()
			if (is.null(private$bms)) {
				return(list(i_b = sample_int_replace_cpp(private$t, private$t), m_vec_b = NULL))
			}
			pair_rows = private$bms$indicies_pairs
			storage.mode(pair_rows) = "integer"
			draw_matching_bootstrap_sample_cpp(
				i_reservoir = integer(0),
				pair_rows = pair_rows,
				n_reservoir = 0L
			)
		},
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			if (should_run_asserts()) {
				if (n %% 4 != 0) {
					stop("DesignFixedMatchingGreedyPairSwitching requires n divisible by 4.")
				}
			}
			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X = private$X[1:n, , drop = FALSE]
			if (is.null(private$bms)) {
				private$bms = compute_binary_match_structure(X, mahal_match = (private$objective == "mahal_dist"))
			}
			pairs_mat = private$bms$indicies_pairs
			storage.mode(pairs_mat) = "integer"
			cpp_n_iter = if (is.infinite(private$n_iter)) -1L else as.integer(private$n_iter)
			w_mat = greedy_design_search_cpp(
				X_raw          = X,
				r              = as.integer(r),
				objective      = private$objective,
				n_iter         = cpp_n_iter,
				indicies_pairs = pairs_mat
			)
			storage.mode(w_mat) = "numeric"
			private$validate_allocation_matrix(w_mat, n = n, r = r)
		},
		validate_allocation_matrix = function(w_mat, n, r){
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (should_run_asserts()) {
				if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < r) {
					stop("DesignFixedMatchingGreedyPairSwitching returned an unexpected allocation matrix shape.")
				}
			}
			w_mat = w_mat[, seq_len(r), drop = FALSE]
			storage.mode(w_mat) = "numeric"
			if (should_run_asserts()) {
				if (any(!is.finite(w_mat)) || any(is.na(w_mat))) {
					stop("DesignFixedMatchingGreedyPairSwitching returned non-finite treatment assignments.")
				}
				if (any(!(w_mat %in% c(0, 1)))) {
					stop("DesignFixedMatchingGreedyPairSwitching returned an invalid treatment assignment matrix.")
				}
				treated_counts = colSums(w_mat)
				if (any(treated_counts != n / 2)) {
					stop("DesignFixedMatchingGreedyPairSwitching returned an unbalanced allocation.")
				}
			}
			w_mat
		}
	)
)
