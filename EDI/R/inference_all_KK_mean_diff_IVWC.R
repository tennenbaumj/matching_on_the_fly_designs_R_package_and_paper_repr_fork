#' Mean-Difference IVWC Inference for KK Designs
#'
#' Fits a compound mean-difference estimator for KK matching-on-the-fly designs.
#' For matched pairs, it uses the average within-pair difference. For reservoir
#' subjects, it uses the treated-minus-control difference in means. The two
#' estimates are combined via inverse-variance weighting.
#'
#'
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
#' @examples
#' \dontrun{
#' seq_des = DesignSeqOneByOneKK14$new(n = 6, response_type = "continuous")
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
#'
#' seq_des_inf = InferenceAllKKMeanDiffIVWC$
#'   new(seq_des)
#' seq_des_inf$compute_estimate()
#' seq_des_inf$compute_asymp_confidence_interval()
#' seq_des_inf$compute_asymp_two_sided_pval()
#' }
InferenceAllKKMeanDiffIVWC = R6::R6Class("InferenceAllKKMeanDiffIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompoundNoParamBootstrap,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		},
		#'
		#' @description Computes the IVWC mean-difference estimate across pairs and reservoir
		#'
		#' @return  The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Whether likelihood-ratio parametric bootstrap is supported.
		#'
		#' @return \code{TRUE}.
		supports_lik_ratio_param_bootstrap = function() TRUE,
		#' @description Whether likelihood-based tests are supported.
		#'
		#' @return \code{TRUE}.
		supports_likelihood_tests = function() TRUE,
		#' @description Simulate responses under a likelihood null model.
		#'
		#' @param spec A likelihood-test specification list.
		#' @param delta The null treatment effect.
		#' @param null_fit The fitted null model object.
		#'
		#' @return A list containing simulated full and null likelihood components.
		simulate_under_lik_null = function(spec, delta, null_fit){
			# Generative model: Matched pairs + Reservoir Gaussian
			KKstats = private$cached_values$KKstats
			n = private$n
			w = spec$w
			m = KKstats$m
			
			# 1. Simulate matched pairs
			ssq_D = KKstats$ssqD_bar * m 
			b_null = as.numeric(null_fit$b)
			y_diffs_sim = rnorm(m, b_null[2], sqrt(ssq_D))
			
			# 2. Simulate reservoir
			ssq_R = KKstats$ssqR
			y_res_sim = numeric(length(KKstats$y_reservoir))
			w_res = KKstats$w_reservoir
			y_res_sim[w_res == 1] = b_null[1] + b_null[2] + rnorm(sum(w_res == 1), 0, sqrt(ssq_R))
			y_res_sim[w_res == 0] = b_null[1] + rnorm(sum(w_res == 0), 0, sqrt(ssq_R))
			
			list(
				full_fit = list(b = b_null, y_diffs = y_diffs_sim, y_res = y_res_sim),
				fit_null = function(d, start = NULL){
					list(b = c(b_null[1], d), y_diffs = y_diffs_sim, y_res = y_res_sim)
				},
				neg_loglik = function(fit){
					b = fit$b
					ll_m = -0.5 * sum((y_diffs_sim - b[2])^2) / ssq_D
					ll_r = -0.5 * sum((y_res_sim - (b[1] + w_res * b[2]))^2) / ssq_R
					-(ll_m + ll_r)
				}
			)
		},
		#' @description Get the likelihood-test specification.
		#'
		#' @return A likelihood-test specification list, or \code{NULL}.
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			if (is.null(private$cached_values$beta_hat_T)) return(NULL)
			KKstats = private$cached_values$KKstats
			list(
				w = private$w,
				full_fit = list(b = c(mean(KKstats$y_reservoir[KKstats$w_reservoir==0], na.rm=T), private$cached_values$beta_hat_T)),
				j = 2L
			)
		},
		#' @description Computes a 1-alpha level frequentist confidence interval
		#'
		#' The compound estimator is treated as asymptotically normal, so this
		#' interval is based on the estimated standard error of the inverse-variance
		#' weighted combination.
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#'
		#' @return  A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		#'
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes a 2-sided p-value
		#'
		#' @param delta   The null difference to test against. For any treatment effect at all this is
		#'   set to zero (the default).
		#'
		#' @return  The approximate frequentist p-value
		#'
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes a 1-alpha level frequentist confidence interval for the randomization test
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#' @param  r  	The number of randomization vectors. The default is 501.
		#' @param  pval_epsilon  		The bisection algorithm tolerance. The default is 0.005.
		#' @param  show_progress  		Show a text progress indicator.
		#' @param ci_search_control Optional randomization-CI search control list passed through 
		#'   to the base method.
		#' @return  A 1 - alpha sized frequentist confidence interval
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			if (should_run_asserts()) {
				if (private$des_obj_priv_int$response_type %in% c("proportion", "count", "survival")) {
					stop("Randomization confidence intervals are not supported for InferenceAllKKMeanDiffIVWC with proportion, count, or survival response types due to inconsistent estimator units on the transformed scale.")
				}
			}
			super$compute_rand_confidence_interval(alpha = alpha, r = r, pval_epsilon = pval_epsilon, show_progress = show_progress, ci_search_control = ci_search_control)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_fast_bootstrap_distr = function(B, i_reservoir, n_reservoir, m, y, w, m_vec) {
			# Only safe for simple additive/linear combinations right now.
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			n = length(y)
			y_mat = matrix(0.0, nrow = n, ncol = B)
			w_mat = matrix(0L, nrow = n, ncol = B)
			m_mat = matrix(0L, nrow = n, ncol = B)
			for (b in 1:B) {
				# Resample reservoir with replacement
				i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
				# For matched pairs, sample which pairs to include (with replacement)
				if (m > 0) {
					pairs_to_include = sample_int_replace_cpp(m, m)
					i_matched_b = integer(0)
					m_vec_b_matched = integer(0)
					for (new_pair_id in 1:m) {
						original_pair_id = pairs_to_include[new_pair_id]
						pair_indices = which(m_vec == original_pair_id)
						i_matched_b = c(i_matched_b, pair_indices)
						m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
					}
				} else {
					i_matched_b = integer(0)
					m_vec_b_matched = integer(0)
				}
				# Combine reservoir and matched indices
				i_b = c(i_reservoir_b, i_matched_b)
				y_mat[, b] = y[i_b]
				w_mat[, b] = w[i_b]
				m_mat[, b] = c(rep(0L, n_reservoir), m_vec_b_matched)
			}
			res = compute_matching_compound_bootstrap_parallel_cpp(
				w_mat,
				m_mat,
				y_mat,
				private$n_cpp_threads(ncol(y_mat))
			)
			return(res)
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			if (delta != 0) return(NULL)
			n = length(y)
			w_mat = as.matrix(permutations$w_mat)
			storage.mode(w_mat) = "integer"
			m_mat = permutations$m_mat
			if (is.null(m_mat)) {
				m_mat = matrix(0L, nrow = n, ncol = ncol(w_mat))
			} else {
				m_mat = as.matrix(m_mat)
				m_mat[is.na(m_mat)] = 0L
				storage.mode(m_mat) = "integer"
			}
			res = compute_matching_compound_distr_parallel_cpp(
				as.numeric(y),
				w_mat,
				m_mat,
				private$n_cpp_threads(ncol(w_mat))
			)
			return(res)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			if (is.null(private$cached_values$KKstats$d_bar)) private$compute_reservoir_and_match_statistics()
			
			KKstats = private$cached_values$KKstats
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			m = KKstats$m
			reservoir_unusable = !is.finite(nRT) || !is.finite(nRC) || nRT <= 1 || nRC <= 1
			no_matches = !is.finite(m) || m <= 1
			
			if (is.null(private$cached_values$beta_hat_T)){
				has_matched_est = is.finite(KKstats$d_bar)
				has_reservoir_est = is.finite(KKstats$r_bar)
				private$cached_values$beta_hat_T =
					if (reservoir_unusable && has_matched_est){
						KKstats$d_bar
					} else if (no_matches && has_reservoir_est){
						KKstats$r_bar
					} else if (has_matched_est && has_reservoir_est){
						KKstats$w_star * KKstats$d_bar + (1 - KKstats$w_star) * KKstats$r_bar #proper weighting
					} else if (has_reservoir_est){
						KKstats$r_bar
					} else if (has_matched_est){
						KKstats$d_bar
					} else {
						NA_real_
					}
			}

			if (estimate_only) return(invisible(NULL))
			
			if (is.null(private$cached_values$KKstats$ssqD_bar)) private$compute_reservoir_and_match_statistics()
			ssqD = private$cached_values$KKstats$ssqD_bar
			ssqR = private$cached_values$KKstats$ssqR
			
			private$cached_values$s_beta_hat_T =
				if (reservoir_unusable){
					# Only matched pairs are usable; fall back to ssqR if ssqD is degenerate
					if (is.finite(ssqD) && ssqD > 0) sqrt(ssqD) else if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
				} else if (no_matches){
					# No matched pairs
					if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
				} else {
					# Combined: require both components to be positive and finite.
					if (!is.finite(ssqD) || ssqD <= 0) {
						if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
					} else if (!is.finite(ssqR) || ssqR <= 0) {
						sqrt(ssqD)
					} else {
						sqrt(ssqR * ssqD / (ssqR + ssqD))
					}
				}
		}
	)))
)

#' @export
InferenceAllKKCompoundMeanDiff <- InferenceAllKKMeanDiffIVWC
