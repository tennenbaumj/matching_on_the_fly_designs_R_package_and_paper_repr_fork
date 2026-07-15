#' A Blocked and Cluster Randomized Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed blocked and
#' cluster randomized experimental design.
#' This design randomizes clusters within specified blocks using the \pkg{randomizr} package.
#'
#' @examples
#' des = DesignFixedBlockedCluster$new(n = 20, response_type = 'continuous', strata_cols = 'x2', cluster_col = 'cl')
#' X = data.frame(x1 = rnorm(20), x2 = factor(rep(1:2, each = 10)), cl = factor(rep(1:10, each = 2)))
#' des$add_all_subjects_to_experiment(X)
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedBlockedCluster = R6::R6Class("DesignFixedBlockedCluster",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a blocked and cluster randomized fixed experimental design
		#'
		#' @param strata_cols 	A character vector of column names to use for stratification (blocks).
		#' @param cluster_col 	The column name in the data that identifies the cluster for each subject.
		#' @param response_type 	The data type of response values.
		#' @param prob_T  The probability of the treatment assignment for each cluster.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param preferred_num_bins_for_continuous_covariate The number of quantile bins to use for continuous strata. Default is 2.
		#' @param num_bins_for_continuous_covariate Deprecated alias for
		#'   `preferred_num_bins_for_continuous_covariate`.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedBlockedCluster` object
		#'
		initialize = function(
				strata_cols,
				cluster_col,
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				preferred_num_bins_for_continuous_covariate = 2,
				num_bins_for_continuous_covariate = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
			if (!is.null(num_bins_for_continuous_covariate)) {
				preferred_num_bins_for_continuous_covariate = num_bins_for_continuous_covariate
			}
			if (should_run_asserts()) {
				assertCharacter(strata_cols, min.len = 1)
				assertCharacter(cluster_col, len = 1)
				assertCount(preferred_num_bins_for_continuous_covariate, positive = TRUE)
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$blocking_capable = TRUE
			private$strata_cols = strata_cols
			private$cluster_col = cluster_col
			private$preferred_num_bins_for_continuous_covariate = preferred_num_bins_for_continuous_covariate
			private$uses_covariates = TRUE
		}
	),
	private = list(
		cluster_col = NULL,
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}

			strata_keys = private$get_strata_keys()

			cluster_ids = as.character(private$Xraw[[private$cluster_col]])

			# Use randomizr::block_and_cluster_ra for canonical blocked and clustered randomization
			w_mat = replicate(r, as.numeric(as.character(randomizr::block_and_cluster_ra(
				blocks = strata_keys,
				clusters = cluster_ids,
				prob = private$prob_T
			))))
			storage.mode(w_mat) = "numeric"
			w_mat
		},
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			n = private$t
			strata_keys = private$get_strata_keys()
			cluster_ids = as.character(private$Xraw[1:n, ][[private$cluster_col]])
			if (is.null(bootstrap_type) || bootstrap_type == "within_blocks") {
				# Resample clusters within each stratum
				unique_strata = unique(strata_keys)
				i_b = unlist(lapply(unique_strata, function(stratum) {
					stratum_idx = which(strata_keys == stratum)
					stratum_group_id = match(cluster_ids[stratum_idx], unique(cluster_ids[stratum_idx]))
					stratum_idx[resample_group_rows_cpp(as.integer(stratum_group_id), length(unique(stratum_group_id)))]
				}), use.names = FALSE)
			} else {
				# Resample blocks (strata) themselves
				strata_group_id = match(strata_keys, unique(strata_keys))
				i_b = resample_group_rows_cpp(as.integer(strata_group_id), length(unique(strata_group_id)))
			}
			list(i_b = as.integer(i_b), m_vec_b = NULL)
		}
	)
)
