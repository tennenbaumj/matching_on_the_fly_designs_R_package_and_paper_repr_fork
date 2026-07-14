#' A Cluster Randomized Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed cluster
#' randomized experimental design.
#' This design randomizes entire groups (clusters) of subjects together using the
#' \pkg{randomizr} package.
#'
#' @examples
#' des = DesignFixedCluster$new(n = 20, response_type = 'continuous', cluster_col = 'cl')
#' X = data.frame(x = rnorm(20), cl = factor(rep(1:5, each = 4)))
#' des$add_all_subjects_to_experiment(X)
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedCluster = R6::R6Class("DesignFixedCluster",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a cluster randomized fixed experimental design
		#'
		#' @param cluster_col 	The column name in the data that identifies the cluster for each subject.
		#' @param response_type 	The data type of response values.
		#' @param prob_T  The probability of the treatment assignment for each cluster.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param model_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedCluster` object
		#'
		initialize = function(
				cluster_col,
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,

				verbose = FALSE,
				missingness_method = "impute",
				model_formula = ~ .,
				seed = NULL
			) {
			if (should_run_asserts()) {
				assertCharacter(cluster_col, len = 1)
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, model_formula, seed = seed)
			private$cluster_col = cluster_col
			private$uses_covariates = TRUE
		}
	),
	private = list(
		cluster_col = NULL,
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			# Assignment is at the cluster level and outcomes are correlated within
			# clusters, so the exchangeable resampling unit is the cluster, not the row.
			n = private$t
			cluster_ids = as.character(private$Xraw[1:n, ][[private$cluster_col]])
			group_id = match(cluster_ids, unique(cluster_ids))
			i_b = resample_group_rows_cpp(as.integer(group_id), length(unique(group_id)))
			list(i_b = as.integer(i_b), m_vec_b = NULL)
		},
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}
			cluster_ids = as.character(private$Xraw[[private$cluster_col]])
			if (should_run_asserts()) {
				if (any(is.na(cluster_ids))){
					stop("Cluster IDs cannot be missing.")
				}
			}

			# Use randomizr::cluster_ra for canonical cluster randomization
			w_mat = replicate(r, as.numeric(as.character(randomizr::cluster_ra(clusters = cluster_ids, prob = private$prob_T))))
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	)
)
