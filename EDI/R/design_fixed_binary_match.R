#' A Binary Match Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed binary match
#' experimental design.
#' This design pairs subjects based on covariate distances and randomizes within pairs.
#' Uses non-bipartite matching via \pkg{nbpMatching} for distance computation.
#'
#' @examples
#' des = DesignFixedBinaryMatch$new(n = 10, response_type = 'continuous')
#' des$add_all_subjects_to_experiment(data.frame(x1 = rnorm(10)))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedBinaryMatch = R6::R6Class("DesignFixedBinaryMatch",
	inherit = DesignFixed,
	public = list(
		#' @description Returns TRUE so the framework pre-generates all w vectors
		#'   once per cell (paying the nbpMatching cost once, reusing across reps).
		supports_batch_w_pregeneration = function() TRUE,
		#' @description Initialize a binary match fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T  The probability of the treatment assignment. Must be 0.5.
		#' @param mahal_match 	Match using Mahalanobis distance (else use Euclidean). Default is \code{TRUE}.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param m Optional integer vector of explicit matched-pair identifiers, one per
		#'   subject. If supplied, `n` must also be supplied, `length(m)` must equal `n`,
		#'   all values must be positive, and each pair ID must occur exactly twice.
		#'   This bypasses the package-computed matching step.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedBinaryMatch` object
		#'
		initialize = function(
				response_type,
				prob_T = 0.5,
				mahal_match = TRUE,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				m = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
				) {
				if (should_run_asserts()) {
					if (prob_T != 0.5){
						stop("Binary match designs only support even treatment allocation (prob_T = 0.5)")
					}
					if (!is.null(m)) {
						if (is.null(n)) {
							stop("When supplying m to DesignFixedBinaryMatch$new(), n must also be supplied.")
						}
						if (length(m) != as.integer(n)) {
							stop("When supplying m to DesignFixedBinaryMatch$new(), length(m) must equal n.")
						}
					}
					# nbpMatching availability is checked lazily in ensure_matching_structure_computed
					# so that workers using pre-computed w vectors never load it unnecessarily.
				}
				super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
				private$blocking_capable = TRUE
				private$matching_capable = TRUE
				private$mahal_match = mahal_match
				private$uses_covariates = TRUE
				if (!is.null(m)) {
					self$set_m(m)
					private$set_binary_match_structure_from_m(m)
				}
			},
		#' @description Assign treatment to all subjects. Ensures matching structure is computed
		#'   even when a pre-computed w vector is injected (bypassing draw_ws_according_to_design).
		#' @param w_precomputed Optional numeric vector of length n.
		assign_w_to_all_subjects = function(w_precomputed = NULL){
			private$ensure_matching_structure_computed()
			super$assign_w_to_all_subjects(w_precomputed)
		}
	),
	private = list(
		mahal_match = NULL,
		bms = NULL,
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
				self$assert_all_subjects_arrived()
			}
			private$ensure_matching_structure_computed()
			n = self$get_n()
			if (is.null(private$m)){
				# No covariates: fall back to BCRD
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}
			# Use the in-house Rcpp search
			w_mat = draw_binary_match_assignments_cpp(
				private$bms$indicies_pairs,
				as.integer(n),
				as.integer(r),
				as.integer(self$num_cores)
			)
			private$validate_allocation_matrix(w_mat, n = n, r = r)
		},
		validate_allocation_matrix = function(w_mat, n, r){
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (should_run_asserts()) {
				if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < 1L) {
					stop("resultsBinaryMatchSearch returned an unexpected allocation matrix shape.")
				}
			}
			storage.mode(w_mat) = "numeric"
			if (should_run_asserts()) {
				if (any(!is.finite(w_mat)) || any(is.na(w_mat))) {
					stop("resultsBinaryMatchSearch returned non-finite treatment assignments.")
				}
				if (any(!(w_mat %in% c(0, 1)))) {
					stop("resultsBinaryMatchSearch returned an invalid treatment assignment matrix.")
				}
				treated_counts = colSums(w_mat)
				if (any(treated_counts != n / 2)) {
					stop("resultsBinaryMatchSearch returned an unbalanced allocation.")
				}
			}
			if (ncol(w_mat) < r){
				w_mat = w_mat[, rep(seq_len(ncol(w_mat)), length.out = r), drop = FALSE]
			}
			w_mat[, seq_len(r), drop = FALSE]
		},
			ensure_matching_structure_computed = function(){
				n = self$get_n()
				if (is.null(private$bms)) {
					if (is.null(private$X) || ncol(private$X) == 0) {
						stop("no covariates provided to run the binary matching algorithm")
					}
					X = private$X[1:n, , drop = FALSE]
					private$bms = compute_binary_match_structure(X, mahal_match = private$mahal_match)
					# Build pair-ID vector m where m[i] = pair index for subject i
					m_vec = integer(n)
					pairs = private$bms$indicies_pairs
					for (i in seq_len(nrow(pairs))){
						m_vec[pairs[i, 1]] = i
						m_vec[pairs[i, 2]] = i
					}
					private$m = m_vec
					private$reset_matching_caches()
				}
				invisible(NULL)
			},
			set_binary_match_structure_from_m = function(m){
				m_vec = as.integer(m)
				if (should_run_asserts()) {
					assertIntegerish(m_vec, lower = 1, any.missing = FALSE, len = self$get_n())
					pair_ids = sort(unique(m_vec))
					pair_sizes = tabulate(match(m_vec, pair_ids), nbins = length(pair_ids))
					if (any(pair_sizes != 2L)) {
						stop("Explicit m for DesignFixedBinaryMatch must define matched pairs only: each pair ID must occur exactly twice.")
					}
				}
				pair_ids = sort(unique(m_vec))
				pairs = t(vapply(pair_ids, function(pid) {
					which(m_vec == pid)
				}, integer(2L)))
				private$bms = list(indicies_pairs = pairs)
				private$reset_matching_caches()
				invisible(NULL)
			}
		)
	)
