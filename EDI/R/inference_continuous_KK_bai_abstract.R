#' Inference based on Maximum Likelihood for KK designs
#'
#' Inference for mean difference. Note that warm starts are disabled for this class
#' as the Bai adjusted t-test is a closed-form estimator and does not benefit from initialization.
#'
#' @details
#' This class requires the \pkg{nbpMatching} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}. Install it manually with
#' \code{install.packages("nbpMatching")} before using this class.
#'
#' @keywords internal
InferenceBaiAdjustedT = R6::R6Class("InferenceBaiAdjustedT",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompoundNoParamBootstrap,
	public = list(
		#' @description Initialize a sequential experimental design estimation and test object after the
		#' sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose                 A flag indicating whether messages should be displayed
		#'   to the user. Default is \code{FALSE}
		#' @param convex_flag       A flag indicating whether the estimator should use a convex
		#'   combination of the Bai et al
		#' matched pairs estimate with the reservoir estimate, or just the Bai et al estimate by its self.
		#'
		initialize = function(des_obj, model_formula = NULL, verbose = TRUE, convex_flag = FALSE){
			if (should_run_asserts()) {
				if (!check_package_installed("nbpMatching")) {
				stop("Package 'nbpMatching' is required for InferenceBaiAdjustedT. Please install it.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			private$fit_warm_start_enabled = FALSE
			private$convex_flag = convex_flag
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#'
		#' @description Computes the appropriate estimate for compound mean difference across pairs and reservoir
		#'
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#'
		#' seq_des_inf = InferenceAllKKMeanDiffIVWC$new(seq_des)
		#' seq_des_inf$compute_estimate()
		#' }
		#'
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			if (is.null(private$cached_values$KKstats$d_bar)) private$compute_reservoir_and_match_statistics()
			KKstats = private$cached_values$KKstats
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			m = KKstats$m
			reservoir_unusable = !is.finite(nRT) || !is.finite(nRC) || nRT <= 1 || nRC <= 1
			no_matches = !is.finite(m) || m <= 1
			has_matched_est = is.finite(KKstats$d_bar)
			has_reservoir_est = is.finite(KKstats$r_bar)
			private$cached_values$beta_hat_T =
				if (reservoir_unusable && has_matched_est){
					KKstats$d_bar
				} else if (no_matches && has_reservoir_est){
					KKstats$r_bar
				} else if (has_matched_est && has_reservoir_est && isTRUE(private$convex_flag)){
					if (is.null(private$cached_values$s_beta_hat_T)){
						private$shared()
					}
					w_star_bai = KKstats$ssqR / (KKstats$ssqR + private$cached_values$bai_var_d_bar)
					w_star_bai * KKstats$d_bar + (1 - w_star_bai) * KKstats$r_bar
				} else if (has_matched_est){
					KKstats$d_bar
				} else if (has_reservoir_est){
					KKstats$r_bar
				} else {
					NA_real_
				}
			private$cached_values$beta_hat_T
		},
		#' @description Computes a 1-alpha level frequentist confidence interval
		#'
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal
		#' (except in the case
		#' of estimat_type "median difference" where a nonparametric bootstrap confidence
		#' interval (see the \code{controlTest::quantileControlTest} method)
		#' is employed. Hence these confidence intervals are asymptotically valid and thus
		#' approximate for any sample size.
		#'
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		#'
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#'
		#' seq_des_inf = InferenceAllKKMeanDiffIVWC$new(seq_des)
		#' seq_des_inf$compute_asymp_confidence_interval()
		#' }
		#'
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (is.null(private$cached_values$beta_hat_T)){
			self$compute_estimate()
			}
			if (is.null(private$cached_values$s_beta_hat_T)){
			private$shared()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes a 2-sided p-value
		#'
		#' @param delta   The null difference to test against. For any treatment effect at all
		#'   this is set to zero (the default).
		#'
		#' @return 	The approximate frequentist p-value
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#'
		#' seq_des_inf = InferenceAllKKMeanDiffIVWC$new(seq_des)
		#' seq_des_inf$compute_asymp_two_sided_pval()
		#' }
		#'
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (is.null(private$cached_values$beta_hat_T)){
			self$compute_estimate()
			}
			if (is.null(private$cached_values$s_beta_hat_T)){
			private$shared()
			}
			2 * stats::pnorm(
			 -abs(private$cached_values$beta_hat_T / private$cached_values$s_beta_hat_T)
			) #approximate by using N(0, 1) distribution
		}
	),
	private = list(
		convex_flag = NULL,
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			
			# Optimization: Ensure matching stats are calculated once
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			KKstats = private$cached_values$KKstats
			m = KKstats$m
			
			# Ensure pairing of pairs is pre-calculated
			halves = private$compute_halves()
			halves_idx = if (nrow(halves) > 0) suppressWarnings(matrix(as.integer(as.matrix(halves[, c(1, 3)])), ncol = 2)) else matrix(0L, 0, 2)
			
			w_mat = permutations$w_mat
			if (!is.null(w_mat) && !is.integer(w_mat)) {
				storage.mode(w_mat) = "integer"
			}
			m_mat = permutations$m_mat
			nsim_local = ncol(w_mat)
			if (is.null(m_mat)) {
				n_subjects = nrow(w_mat)
				m_vec = as.integer(private$des_obj_priv_int$m)
				if (length(m_vec) != n_subjects) m_vec = rep(0L, n_subjects)
				m_mat = matrix(rep(m_vec, nsim_local), nrow = n_subjects, ncol = nsim_local)
			}
			if (!is.null(m_mat) && !is.integer(m_mat)) {
				storage.mode(m_mat) = "integer"
			}
			
			# The Bai statistic is the treatment estimate itself (beta_hat_T)
			res = compute_bai_distr_parallel_cpp(
				w_mat,
				m_mat,
				as.numeric(y),
				as.numeric(delta),
				halves_idx,
				isTRUE(private$convex_flag),
				private$n_cpp_threads(ncol(w_mat))
			)
			return(res)
		},
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			i = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			i
		},
		shared = function(estimate_only = FALSE){
				if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
				if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			if (is.null(private$cached_values$KKstats$d_bar)) private$compute_reservoir_and_match_statistics()
			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			reservoir_unusable = !is.finite(nRT) || !is.finite(nRC) || nRT <= 1 || nRC <= 1
			no_matches = !is.finite(m) || m <= 1
			ssqD = NA_real_
			if (is.finite(m) && m > 0){
				private$cached_values$bai_var_d_bar = private$compute_bai_variance_for_pairs() / m
				ssqD = private$cached_values$bai_var_d_bar
			}
			ssqR = KKstats$ssqR
			private$cached_values$s_beta_hat_T =
				if (reservoir_unusable){
					if (is.finite(ssqD) && ssqD > 0) sqrt(ssqD) else if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
				} else if (no_matches){
					if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else if (is.finite(ssqD) && ssqD > 0) sqrt(ssqD) else NA_real_
				} else if (isTRUE(private$convex_flag) && is.finite(ssqD) && ssqD > 0 && is.finite(ssqR) && ssqR > 0){
					sqrt(ssqR * ssqD / (ssqR + ssqD))
				} else if (is.finite(ssqD) && ssqD > 0){
					sqrt(ssqD)
				} else if (is.finite(ssqR) && ssqR > 0){
					sqrt(ssqR)
				} else {
					NA_real_
				}
		},
		compute_bai_variance_for_pairs = function(){
			pairs_df = data.frame(
			pair_id = 1 : private$cached_values$KKstats$m,
			yT = private$cached_values$KKstats$yTs_matched,
			yC = private$cached_values$KKstats$yCs_matched,
			d_i = private$cached_values$KKstats$y_matched_diffs
			)
			halves = private$compute_halves()
			delta_sq = mean(pairs_df$d_i)^2
			tau_sq = mean(pairs_df$d_i^2)
			# lambda_squ^2 term
			lambda_squ = 0
			if (nrow(halves) > 0){
			halves_idx = suppressWarnings(matrix(as.integer(as.matrix(halves[, c(1, 3)])), ncol = 2))
			halves_idx = halves_idx[complete.cases(halves_idx), , drop = FALSE]
			if (nrow(halves_idx) > 0){
				lambda_squ = compute_lambda_squ_cpp(pairs_df$d_i, halves_idx)
			}
			}
			v_sq = tau_sq - (lambda_squ + delta_sq) / 2
			# The variance cannot be negative.
			max(v_sq, 1e-8)
		},
		compute_halves = function(){
			if (!is.null(private$cached_values$halves)) return(private$cached_values$halves)
			
			m = private$cached_values$KKstats$m
			if (is.null(m) || !is.finite(m) || m < 2) return(data.frame()) # Cannot make pairs of pairs if there's < 2 pairs
			X = private$get_X()
			pair_avg = compute_pair_averages_cpp(X, private$des_obj_priv_int$m, m)
			weights = private$des_obj_priv_int$covariate_weights
			if (is.null(weights) || length(weights) != ncol(pair_avg)){
			weights = numeric()
			}
			dist_mat = compute_pair_distance_matrix_cpp(pair_avg, weights)
			# Use nbpMatching to find the optimal pairing of the pairs to minimize total distance
			dist_obj = suppressWarnings(nbpMatching::distancematrix(dist_mat))
			match_obj = suppressWarnings(nbpMatching::nonbimatch(dist_obj))
			halves = match_obj$halves
			# If there's an odd number of pairs, remove the "ghost" match
			if (m %% 2 == 1){
			ghost_row = which(halves[, 3] == "ghost")
			if (length(ghost_row) > 0) halves = halves[-ghost_row, ]
			}
			private$cached_values$halves = halves
			halves
		}
	)
)
