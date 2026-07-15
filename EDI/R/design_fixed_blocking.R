#' A stratified blocking Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed stratified
#' blocking experimental design.
#'
#' @examples
#' des = DesignFixedBlocking$new(n = 20, response_type = 'continuous', strata_cols = 'x2', equal_block_sizes = FALSE)
#' X = data.frame(x1 = rnorm(20), x2 = factor(rep(1:2, 10)))
#' des$add_all_subjects_to_experiment(X)
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedBlocking = R6::R6Class("DesignFixedBlocking",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a fixed stratified blocking experimental design
		#'
		#' @param strata_cols A character vector of column names to use for stratification.
		#'   If `NULL` (the default), all available covariate columns are used.
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param  prob_T  Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param  n  		The sample size.
		#' @param preferred_num_bins_for_continuous_covariate The number of quantile bins to use for continuous strata. Default is 2.
		#' @param B_target The target number of blocks. Columns from `strata_cols`
		#'   are added greedily in order, each column being included only if it does not push
		#'   the total number of unique blocks beyond this target. For categorical covariates
		#'   their natural levels are used; for continuous covariates
		#'   `preferred_num_bins_for_continuous_covariate` quantile bins are used. Earlier columns
		#'   are always preferred over later ones. The default is `floor(sqrt(n))` when `n`
		#'   is known at construction time, or is resolved to `floor(sqrt(n))` when subjects
		#'   are added. Set `B_target = NULL` to use all columns unconditionally.
		#'   Set `exact_num_blocks = TRUE` to hard fail if the final key construction
		#'   does not produce exactly `B_target` blocks.
			#' @param exact_num_blocks Whether to require the greedy key construction to produce
			#'   exactly `B_target` blocks. Default `FALSE`.
			#' @param equal_block_sizes Whether to require all blocks to have the same number of
			#'   subjects. Default `TRUE`. When `TRUE` and both `n` and `B_target` are known at
			#'   construction time, an error is raised immediately if `n` is not divisible by
			#'   `B_target`. A second check fires when subjects are added: if the covariate-based
			#'   strata produce unequal block counts the design errors at that point. Set to
			#'   `FALSE` to allow unequal blocks (note that `InferenceIncidCMH` and
			#'   `InferenceIncidExtendedRobins` still require equal block sizes regardless).
			#' @param m Optional integer vector of explicit block identifiers, one per subject.
			#'   If supplied, `n` must also be supplied and `length(m)` must equal `n`.
			#'   The constructor then records this blocking structure immediately via
			#'   `set_m()`, bypassing covariate-derived strata construction.
			#' @param verbose A flag for verbosity.
			#' @param missingness_method How to handle missing values in covariates.
			#' @param design_formula A formula object.
			#' @param seed Integer seed for reproducibility.
		#'
		#' @return  A new `DesignFixedBlocking` object
			initialize = function(
						strata_cols = NULL,
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						preferred_num_bins_for_continuous_covariate = 2,
							B_target = if (!is.null(n)) max(1L, floor(sqrt(n))) else NA_integer_,
							exact_num_blocks = FALSE,
							equal_block_sizes = TRUE,
							m = NULL,
							verbose = FALSE,
					missingness_method = "impute",
					design_formula = ~ .,
					seed = NULL) {
				if (should_run_asserts()) {
					if (!is.null(strata_cols)) assertCharacter(strata_cols, min.len = 1)
				assertCount(preferred_num_bins_for_continuous_covariate, positive = TRUE)
				if (!is.null(B_target) && !is.na(B_target)) assertCount(B_target, positive = TRUE)
					assertLogical(exact_num_blocks, len = 1)
					assertLogical(equal_block_sizes, len = 1)
					if (!is.null(m)) {
						if (is.null(n)) {
							stop("When supplying m to DesignFixedBlocking$new(), n must also be supplied.")
						}
						if (length(m) != as.integer(n)) {
							stop("When supplying m to DesignFixedBlocking$new(), length(m) must equal n.")
						}
					}
					if (isTRUE(equal_block_sizes) && !is.null(n) && !is.null(B_target) && !is.na(B_target)) {
						if (n %% B_target != 0L) {
							stop("equal_block_sizes = TRUE requires n to be divisible by B_target, but n = ",
							n, " is not divisible by B_target = ", B_target, ".")
					}
					private$assert_min_block_size(n, B_target)
				}
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$blocking_capable = TRUE
			private$strata_cols = strata_cols
			private$preferred_num_bins_for_continuous_covariate = preferred_num_bins_for_continuous_covariate
				private$B_target = B_target
				private$exact_num_blocks = exact_num_blocks
				private$equal_block_sizes = equal_block_sizes
				private$uses_covariates = TRUE
				if (!is.null(m)) {
					self$set_m(m)
				}
			}
	),
	private = list(
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}

			strata_keys = private$get_strata_keys()

			# Use randomizr::block_ra for canonical stratified blocking if available,
			# or fallback to our C++ implementation.
			if (check_package_installed("randomizr")) {
				w_mat = replicate(r, as.numeric(as.character(randomizr::block_ra(blocks = strata_keys, prob = private$prob_T))))
				storage.mode(w_mat) = "numeric"
				return(w_mat)
			}

			unique_keys = unique(strata_keys)
			strata_indices = lapply(unique_keys, function(key) which(strata_keys == key))

			res = generate_permutations_blocking_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T),
				strata_indices
			)
			w_mat = res$w_mat
			storage.mode(w_mat) = "numeric"
			w_mat
		},
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			strata_keys = private$get_strata_keys()
			if (is.null(bootstrap_type) || bootstrap_type == "within_blocks") {
				strata_ids = match(strata_keys, unique(strata_keys))
				list(i_b = stratified_bootstrap_indices_cpp(as.integer(strata_ids)), m_vec_b = NULL)
			} else {
				group_id = match(strata_keys, unique(strata_keys))
				i_b = resample_group_rows_cpp(as.integer(group_id), length(unique(group_id)))
				list(i_b = as.integer(i_b), m_vec_b = NULL)
			}
		}
	)
)
