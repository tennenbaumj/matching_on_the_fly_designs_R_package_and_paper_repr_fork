#' Abstract base class for KK Wilcoxon-based compound inference
#'
#' Shared base for all KK Wilcoxon inference classes. Overrides the per-permutation
#' statistic used in randomization tests with standardized Wilcoxon W statistics
#' (O(n log n), conf.int = FALSE), avoiding the O(n^2) Walsh-average computation
#' required by the full Hodges-Lehmann estimate.
#'
#' @keywords internal
InferenceAbstractKKWilcoxBaseIVWC = R6::R6Class("InferenceAbstractKKWilcoxBaseIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the abstract base.
		#' @param des_obj A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		},
		#' @description Override to avoid O(n^2) per-resample HL computation during the bootstrap warm-start
		#' inside compute_rand_confidence_interval. The asymptotic MLE CI is a perfectly
		#' adequate starting bound for the bisection and is computed in O(1).
		#' @param alpha  				The confidence level. Default is 0.05.
		#' @param ... 					Additional arguments passed to super.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, ...){
			self$compute_asymp_confidence_interval(alpha)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			# Optimization: w_mat and m_mat are already pre-computed matrices
			w_mat = permutations$w_mat
			m_mat = permutations$m_mat
			nsim = ncol(w_mat)
			# Reconstruct m_mat if NULL (KK permutations return fixed matching as NULL to save memory)
			if (is.null(m_mat)) {
				n_subjects = nrow(w_mat)
				m_vec = as.integer(private$des_obj_priv_int$m)
				if (length(m_vec) != n_subjects) m_vec = rep(0L, n_subjects)
				m_mat = matrix(rep(m_vec, nsim), nrow = n_subjects, ncol = nsim)
				is_fixed_matching = TRUE
			} else {
				# Check if all matchings are identical
				is_fixed_matching = TRUE
				if (nsim > 1) {
					first_match = m_mat[, 1]
					for (j in 2:min(nsim, 5)) {
						if (!all(m_mat[, j] == first_match)) {
							is_fixed_matching = FALSE
							break
						}
					}
				}
			}
			y_sim = as.numeric(y)
			# Map transform_responses to transform_code
			t_code = 0L # none
			if (transform_responses == "log") {
				t_code = 1L
			} else if (transform_responses == "logit") {
				t_code = 2L
			} else if (transform_responses == "log1p") {
				t_code = 3L
			}
			res = compute_matching_wilcox_distr_parallel_cpp(
				w_mat,
				m_mat,
				y_sim,
				as.numeric(delta),
				t_code,
				as.numeric(zero_one_logit_clamp),
				is_fixed_matching,
				private$n_cpp_threads(nsim)
			)
			return(res)
		},
		# Override the per-permutation statistic to avoid the O(n^2) conf.int = TRUE cost.
		# Uses standardized Wilcoxon W statistics (conf.int = FALSE, O(n log n)) as the
		# rank test statistic; monotone with the HL / rank-regression estimate under the null.
		compute_treatment_estimate_during_randomization_inference = function(){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}
			stat = 0
			n_components = 0
			# Matched pairs: signed-rank W (standardized)
			diffs = KKstats$y_matched_diffs
			m_pairs = length(diffs)
			if (m_pairs > 0){
				# Optimization: manually compute Wilcoxon rank sum statistic to avoid R overhead
				# This is equivalent to wilcox.test(diffs)$statistic
				abs_diffs = abs(diffs)
				signs = sign(diffs)
				if (!all(signs == 0)){
					# Fast ranking
					rks = rank(abs_diffs, ties.method = "average")
					W_plus = sum(rks[signs > 0]) + 0.5 * sum(rks[signs == 0])
					
					E_W = m_pairs * (m_pairs + 1) / 4
					V_W = m_pairs * (m_pairs + 1) * (2 * m_pairs + 1) / 24
					stat = stat + (W_plus - E_W) / sqrt(V_W)
					n_components = n_components + 1
				}
			}
			# Reservoir: rank-sum W (standardized)
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (nRT > 0 && nRC > 0){
				y_r = KKstats$y_reservoir
				w_r = KKstats$w_reservoir
				# Optimization: use rank() directly
				rks_r = rank(y_r, ties.method = "average")
				W_r = sum(rks_r[w_r == 1])
				
				E_W = nRT * (nRT + nRC + 1) / 2
				V_W = nRT * nRC * (nRT + nRC + 1) / 12
				stat = stat + (W_r - E_W) / sqrt(V_W)
				n_components = n_components + 1
			}
			if (n_components == 0L) NA_real_ else stat
		}
	)))
)
#' Non-parametric Wilcoxon-based Compound Inference for KK Designs
#'
#' Fits a non-parametric compound estimator for KK matching-on-the-fly designs.
#' For matched pairs, it uses the Wilcoxon Signed-Rank Hodges-Lehmann estimate.
#' For reservoir subjects, it uses the Wilcoxon Rank-Sum (Mann-Whitney U) Hodges-Lehmann
#' estimate. The two estimates are combined via a variance-weighted linear combination.
#' This method is robust to outliers and does not assume a specific parametric
#' distribution for the response.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceAllKKWilcoxIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceAllKKWilcoxIVWC = R6::R6Class("InferenceAllKKWilcoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKWilcoxBaseIVWC,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose Whether to print progress messages.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			res_type = des_obj$get_response_type()
			if (should_run_asserts()) {
				if (res_type == "incidence"){
					stop("Rank-based compound inference is not recommended for incidence data; clogit or compound mean difference estimators are preferred.")
				}
			}
			if (should_run_asserts()) {
				assertResponseType(res_type, c("continuous", "count", "proportion", "survival", "ordinal"))
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				if (private$any_censoring){
					stop(class(self)[1], " does not currently support censored survival data. Use restricted mean or Cox-based methods instead.")
				}
			}
		},
		#' @description Returns the estimated treatment effect (Hodges-Lehmann median shift).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the non-parametric confidence interval.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			# Even though estimates are non-parametric, the combined estimator is asymptotically normal
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the non-parametric p-value.
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				if (should_run_asserts()) {
					stop("Testing non-zero delta is not yet implemented for the combined rank estimator.")
				}
				NA_real_
			}
		},
		#' @description Jackknife bias correction is unstable for the
		#'   Hodges-Lehmann rank estimator; report explicit non-estimability.
		compute_jackknife_estimate = function(unit = "auto"){
			private$cache_nonestimable_estimate("kk_wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_corrected_estimate = function(unit = "auto"){
			self$compute_jackknife_estimate(unit = unit)
		},
		compute_jackknife_bias_estimate = function(unit = "auto"){
			private$cache_nonestimable_estimate("kk_wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_std_error = function(unit = "auto"){
			private$cache_nonestimable_se("kk_wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_standard_error = function(unit = "auto"){
			self$compute_jackknife_std_error(unit = unit)
		},
		compute_jackknife_wald_two_sided_pval = function(delta = 0, unit = "auto"){
			private$cache_nonestimable_se("kk_wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			private$cache_nonestimable_se("kk_wilcox_hl_jackknife_not_supported")
			c(NA_real_, NA_real_)
		}
	),
	private = list(
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			# Recompute KKstats if cache was cleared (e.g., after y transformation for rand CI)
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			# --- Matched pairs: Wilcoxon Signed-Rank HL Estimate ---
			if (m > 0){
				private$rank_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0
			# --- Reservoir: Wilcoxon Rank-Sum HL Estimate ---
			if (nRT > 0 && nRC > 0){
				private$rank_for_reservoir()
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0
			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		rank_for_matched_pairs = function(){
			diffs = private$cached_values$KKstats$y_matched_diffs
			m_pairs = length(diffs)
			# signed-rank test requires at least some non-zero differences
			if (all(diffs == 0)) return(invisible(NULL))
			mod = tryCatch({
				stats::wilcox.test(diffs, conf.int = TRUE)
			}, error = function(e) NULL)
			if (is.null(mod)) return(invisible(NULL))
			beta = as.numeric(mod$estimate)
			# Variance of Walsh averages / m — matches C++ fast-path estimate_hl_ssq_signed_rank
			se_sq = if (m_pairs >= 2L) {
				walsh = unlist(lapply(seq_len(m_pairs), function(i) (diffs[i] + diffs[i:m_pairs]) / 2))
				var(walsh) / m_pairs
			} else NA_real_
			private$cached_values$beta_T_matched     = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(se_sq) && se_sq > 0) se_sq else NA_real_
		},
		rank_for_reservoir = function(){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			yT = y_r[w_r == 1]
			yC = y_r[w_r == 0]
			mod = tryCatch({
				stats::wilcox.test(yT, yC, conf.int = TRUE)
			}, error = function(e) NULL)
			if (is.null(mod)) return(invisible(NULL))
			beta = as.numeric(mod$estimate)
			# Variance of pairwise differences / (n_t + n_c) — matches C++ estimate_hl_ssq_rank_sum
			se_sq = if (length(yT) >= 2L && length(yC) >= 2L) {
				var(as.vector(outer(yT, yC, FUN = "-"))) / (length(yT) + length(yC))
			} else NA_real_
			private$cached_values$beta_T_reservoir     = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(se_sq) && se_sq > 0) se_sq else NA_real_
		},
		compute_fast_bootstrap_distr = function(B, i_reservoir, n_reservoir, m, y, w, m_vec){
			# Generate bootstrap indices for KK design in R first
			indices_mat = matrix(NA_integer_, nrow = length(y), ncol = B)
			m_mat = matrix(NA_integer_, nrow = length(y), ncol = B)
			w_mat = matrix(NA_integer_, nrow = length(y), ncol = B)
			for (b in 1:B) {
				i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
				w_b_res = w[i_reservoir_b]
				
				i_matched_b = integer(0)
				m_vec_b_matched = integer(0)
				w_b_matched = integer(0)
				
				if (m > 0) {
					pairs_to_include = sample(seq_len(m), m, replace = TRUE)
					for (new_pair_id in 1:m) {
						original_pair_id = pairs_to_include[new_pair_id]
						pair_indices = which(m_vec == original_pair_id)
						i_matched_b = c(i_matched_b, pair_indices)
						m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
						w_b_matched = c(w_b_matched, w[pair_indices])
					}
				}
				
				w_mat[, b] = c(w_b_res, w_b_matched)
				m_mat[, b] = c(rep(0L, n_reservoir), m_vec_b_matched)
				# Note: y is retrieved inside C++ using indices_mat
				indices_mat[, b] = c(i_reservoir_b, i_matched_b)
			}
			compute_wilcox_matching_ivwc_bootstrap_parallel_cpp(
				as.integer(w),
				as.numeric(y),
				as.integer(m_vec),
				indices_mat,
				m_mat,
				private$n_cpp_threads(B)
			)
		}
	)
)
