#' Pocock & Simon's Minimization Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for a Pocock & Simon
#' sequential experimental design.
#' This design minimizes the imbalance across treatments for multiple covariates.
#'
#' @examples
#' seq_des = DesignSeqOneByOnePocockSimon$new(strata_cols = 'x1', n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = factor(1, levels=1:2)))
#' @export
DesignSeqOneByOnePocockSimon = R6::R6Class("DesignSeqOneByOnePocockSimon",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description Initialize a Pocock & Simon sequential experimental design
		#'
		#' @param strata_cols     The names of the covariates to be used for minimization. These
		#'   must be factor or categorical variables.
		#' @param weights 		A numeric vector of weights for each covariate. Defaults to 1 for all.
		#' @param p_best          The probability of assigning the treatment that minimizes the
		#'   imbalance. Defaults to 0.8.
		#' @param response_type 	The data type of response values.
		#' @param prob_T  The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignSeqOneByOnePocockSimon` object
		#'
		initialize = function(
				strata_cols,
				weights = NULL,
				p_best = 0.8,
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
				assertCharacter(strata_cols, min.len = 1)
				assertNumeric(p_best, lower = 0.5, upper = 1)
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			
			private$strata_cols = strata_cols
			private$p_best = p_best
			private$uses_covariates = TRUE
			
			if (is.null(weights)){
				private$weights = rep(1, length(strata_cols))
			} else {
				if (should_run_asserts()) {
					assertNumeric(weights, len = length(strata_cols), lower = 0)
				}
				private$weights = weights
			}
		},
		#' @description Assign the next subject to a treatment group using minimization.
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			private$ensure_factor_metadata()
			subject_levels_idx = private$get_subject_levels_idx(private$Xraw[private$t, ])
			
			if (is.null(private$counts)){
				private$counts = matrix(0, nrow = private$num_levels_total, ncol = 2)
			}
			
			# Call Rcpp function that assigns and updates counts in-place
			pocock_simon_assign_and_update_cpp(
				private$counts, 
				as.integer(subject_levels_idx), 
				private$weights, 
				private$p_best, 
				private$prob_T
			)
		}
	),
	private = list(
		draw_ws_raw = function(r = 100){
			private$ensure_factor_metadata()
			n = self$get_n()
			x_levels_matrix = matrix(NA_integer_, nrow = n, ncol = length(private$strata_cols))
			for (i in 1 : n){
				x_levels_matrix[i, ] = private$get_subject_levels_idx(private$Xraw[i, ])
			}
			generate_permutations_pocock_simon_cpp(
				x_levels_matrix,
				as.integer(private$num_levels_total),
				private$weights,
				private$p_best,
				private$prob_T,
				as.integer(r)
			)$w_mat
		},
		p_best = NULL,
		weights = NULL,
		counts = NULL,
		num_levels_total = NULL,
		strata_level_rows = NULL,
		draw_bootstrap_indices = function(bootstrap_type = NULL) {
			list(i_b = sample_int_replace_cpp(private$t, private$t), m_vec_b = NULL)
		},
		ensure_factor_metadata = function(){
			if (is.null(private$strata_level_rows)) private$strata_level_rows = vector("list", length(private$strata_cols))
			if (length(private$strata_level_rows) != length(private$strata_cols)) {
				private$strata_level_rows = vector("list", length(private$strata_cols))
			}
			names(private$strata_level_rows) = private$strata_cols
			next_row = 1L
			for (col in private$strata_cols) {
				row_map = private$strata_level_rows[[col]]
				if (is.null(row_map)) row_map = integer(0)
				col_vals = private$Xraw[[col]]
				col_keys = ifelse(is.na(col_vals), "NA", as.character(col_vals))
				new_levels = setdiff(unique(col_keys), names(row_map))
				if (length(new_levels) > 0L) {
					new_rows = seq.int(next_row, length.out = length(new_levels))
					names(new_rows) = new_levels
					row_map = c(row_map, new_rows)
				}
				private$strata_level_rows[[col]] = row_map
				if (length(row_map) > 0L) next_row = max(unname(row_map)) + 1L
			}
			private$num_levels_total = max(0L, next_row - 1L)
			if (!is.null(private$counts) && nrow(private$counts) < private$num_levels_total) {
				expanded = matrix(0, nrow = private$num_levels_total, ncol = 2L)
				expanded[seq_len(nrow(private$counts)), ] = private$counts
				private$counts = expanded
			}
		},
		get_subject_levels_idx = function(x_row){
			private$ensure_factor_metadata()
			vapply(private$strata_cols, function(col) {
				key = if (is.na(x_row[[col]])) "NA" else as.character(x_row[[col]])
				row_idx = private$strata_level_rows[[col]][[key]]
				if (should_run_asserts()) {
					if (is.null(row_idx) || !is.finite(row_idx)) {
						stop("Unknown strata level encountered for Pocock-Simon column ", col, ": ", key)
					}
				}
				as.integer(row_idx)
			}, integer(1))
		}
	)
)
