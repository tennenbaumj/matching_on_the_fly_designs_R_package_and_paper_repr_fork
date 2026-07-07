#' A stratified permuted block Sequential Design (SPBR)
#'
#' An R6 Class encapsulating the data and functionality for a stratified permuted
#' block sequential experimental design.
#' This design ensures balance within specified strata using blocks of a fixed size.
#'
#' @examples
#' seq_des = DesignSeqOneByOneSPBR$new(strata_cols = 'x1', n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = factor(1, levels=1:2)))
#' @export
DesignSeqOneByOneSPBR = R6::R6Class("DesignSeqOneByOneSPBR",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description Initialize a stratified permuted block sequential experimental design
		#'
		#' @param strata_cols A character vector of column names to use for stratification.
		#' @param block_size The size of the permuted blocks.
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param  prob_T  Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param  n  		The sample size.
		#' @param verbose A flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param model_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return  A new `DesignSeqOneByOneSPBR` object
		initialize = function(
						strata_cols,
						block_size = 4,
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,

						verbose = FALSE,
				missingness_method = "impute",
				model_formula = ~ .,
				seed = NULL
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, model_formula, seed = seed)
			private$blocking_capable = TRUE
			private$strata_cols = strata_cols
			private$block_size = as.integer(block_size)
			private$uses_covariates = TRUE
			private$strata_states = new.env(parent = emptyenv())
			
			if (should_run_asserts()) {
				if (abs(block_size * prob_T - round(block_size * prob_T)) > 1e-10) {
					stop("block_size must result in an integer number of treatment assignments (block_size * prob_T).")
				}
			}
		},
		#' @description Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			x_new = private$Xraw[private$t, ]
			key = private$get_strata_key(x_row = x_new)
			
			if (is.null(private$strata_states[[key]]) || length(private$strata_states[[key]]) == 0) {
				n_T = round(private$block_size * private$prob_T)
				n_C = private$block_size - n_T
				new_block = sample(c(rep(1, n_T), rep(0, n_C)))
				private$strata_states[[key]] = new_block
			}
			
			block = private$strata_states[[key]]
			w_t = block[1]
			private$strata_states[[key]] = block[-1]
			w_t
		}
	),
	private = list(
		strata_cols = NULL,
		block_size = NULL,
		strata_states = NULL,
		draw_ws_raw = function(r = 100){
			strata_keys = vapply(1:private$t, function(i) {
				private$get_strata_key(private$Xraw[i, ])
			}, character(1))

			generate_permutations_spbr_cpp(
				as.character(unname(strata_keys)),
				as.integer(private$block_size),
				as.numeric(private$prob_T),
				as.integer(r)
			)$w_mat
		},
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			strata_keys = vapply(1:private$t, function(i) {
				private$get_strata_key(private$Xraw[i, ])
			}, character(1))
			if (is.null(bootstrap_type) || bootstrap_type == "within_blocks") {
				strata_ids = match(strata_keys, unique(strata_keys))
				list(i_b = stratified_bootstrap_indices_cpp(as.integer(strata_ids)), m_vec_b = NULL)
			} else {
				group_id = match(strata_keys, unique(strata_keys))
				i_b = resample_group_rows_cpp(as.integer(group_id), length(unique(group_id)))
				list(i_b = as.integer(i_b), m_vec_b = NULL)
			}
		},
		get_strata_key = function(x_row) {
			vals = vapply(private$strata_cols, function(col) {
				val = x_row[[col]]
				if (is.na(val)) "NA" else as.character(val)
			}, character(1))
			paste(vals, collapse = "|")
		}
	)
)
