#' A Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#'
#' @examples
#' seq_des = DesignSeqOneByOneKK14$new(n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' @export
DesignSeqOneByOneKK14 = R6::R6Class("DesignSeqOneByOneKK14",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description Initialize a KK14 sequential experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param  prob_T  Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param  n  		The sample size.
		#' @param verbose A flag for verbosity.
		#' @param lambda The penalty parameter for covariate imbalance.
		#' @param t_0_pct The percentage of subjects to allocate before matching begins.
		#' @param morrison If TRUE, use Morrison's method for matching.
		#' @param p The number of covariates to use for matching.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param model_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return  A new `DesignSeqOneByOneKK14` object
		initialize = function(
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,

						verbose = FALSE,
						lambda = NULL,
						t_0_pct = NULL,
						morrison = FALSE,
						p = NULL,
						missingness_method = "impute",
						model_formula = ~ .,
						seed = NULL
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, model_formula, seed = seed)
			private$blocking_capable = TRUE
			private$matching_capable = TRUE
			private$uses_covariates = TRUE
			private$lambda = if (is.null(lambda)) 0.1 else lambda
			private$t_0_pct = if (is.null(t_0_pct)) 0.35 else t_0_pct
			private$morrison = morrison
			private$p = p
		},
		#' @description Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			wt = 	if (private$too_early_to_match()){
						#we're early, so randomize
						private$m[private$t] = 0 #zero means "reservoir", >0 means match number
						private$assign_wt_Bernoulli()
					} else {
						all_subject_data = private$compute_all_subject_data()
						#compute inverse sample covariance of all past subjects (with eps regularization)
						S_xs_inv = solve(var(all_subject_data$X_prev) + diag(.Machine$double.eps, all_subject_data$rank_prev), tol = .Machine$double.xmin)
						#now find the best match in the reservoir
						reservoir_indices = which(private$m[1 : (private$t - 1)] == 0)
						if (length(reservoir_indices) == 0){
							private$m[private$t] = 0
							return(private$assign_wt_Bernoulli())
						}
						sqd_distances_times_two = compute_proportional_mahal_distances_cpp(
							all_subject_data$xt_prev,
							all_subject_data$X_prev,
							as.integer(reservoir_indices),
							S_xs_inv
						)
						#compute F-test threshold
						F_crit = qf(private$compute_lambda(), all_subject_data$rank_prev, private$t - all_subject_data$rank_prev)
						n = self$get_n()
						T_cutoff_sq = all_subject_data$rank_prev * (n - 1) / (n - all_subject_data$rank_prev) * F_crit
						min_sqd_dist_index = which(sqd_distances_times_two == min(sqd_distances_times_two))
						if (length(min_sqd_dist_index) > 1) min_sqd_dist_index = min_sqd_dist_index[1]
						if (sqd_distances_times_two[min_sqd_dist_index] < T_cutoff_sq){
							#we matched!
							match_num = max(private$m) + 1
							#update previous subject's match number
							private$m[reservoir_indices[min_sqd_dist_index]] = match_num
							#update current subject's match number
							private$m[private$t] = match_num
							#assign opposite
							1 - private$w[reservoir_indices[min_sqd_dist_index]]
						} else { #otherwise, randomize and add it to the reservoir
							private$m[private$t] = 0
							private$assign_wt_Bernoulli()
						}
					}
			if (should_run_asserts()) {
				if (is.na(private$m[private$t])){
					stop("no match data recorded")
				}
			}
			wt
		}
	),
	private = list(
		m = NULL,
		draw_ws_raw = function(r = 100){
			generate_permutations_matching_cpp(
				as.integer(private$m),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		},
		lambda = NULL,
		t_0_pct = NULL,
		morrison = NULL,
		p = NULL,
		compute_lambda = function(){
			private$lambda
		},
		too_early_to_match = function(){
			private$t <= private$t_0_pct * private$n ||
				is.null(private$X) ||
				ncol(as.matrix(private$X)) == 0L
		}
	)
)
