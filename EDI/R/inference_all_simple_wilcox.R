#' Wilcoxon Rank-Sum Inference
#'
#' Fits a Wilcoxon rank-sum (Hodges-Lehmann) inference for continuous responses
#' using the treatment indicator.
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
#' seq_des_inf = InferenceAllSimpleWilcox$new(seq_des)
#' seq_des_inf$compute_estimate()
#' }
#' @export
InferenceAllSimpleWilcox = R6::R6Class("InferenceAllSimpleWilcox",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj  A completed \code{DesignSeqOneByOne} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose      Whether to print progress messages. Default \code{FALSE}.
		#' @param max_resample_attempts Maximum number of times a single bootstrap replicate
		#'   may be redrawn when the drawn sample fails validity screening. If all attempts
		#'   fail the replicate is recorded as \code{NA}, silently reducing the effective \code{B}.
		#'   Must be a positive integer. Default \code{50L}.
		#' @param smart_cold_start_default Flag for consistent API.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, max_resample_attempts = 50L, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertCount(max_resample_attempts, positive = TRUE)
			}
			res_type = des_obj$get_response_type()
			if (should_run_asserts()) {
				if (res_type == "incidence"){
					stop(
						"Wilcoxon rank-sum inference is not implemented for incidence (binary) ",
						"responses: the Hodges-Lehmann estimator is degenerate (almost always 0) ",
						"on 0/1 data. Use InferenceAllSimpleMeanDiff or a clogit estimator instead."
					)
				}
			}
			if (should_run_asserts()) {
				assertResponseType(res_type, c("continuous", "count", "proportion", "survival", "ordinal"))
			}
			super$initialize(des_obj = des_obj, verbose = verbose, harden = TRUE, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$max_resample_attempts = max_resample_attempts
			if (private$any_censoring){
				stop(
					"Wilcoxon rank-sum inference does not support censored survival data. ",
					"Use InferenceSurvivalGehanWilcox for censored survival outcomes."
				)
			}
		},
		#' @description Returns the Hodges-Lehmann pseudo-median of all pairwise treatment-minus-control
		#' differences.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the weighted Hodges-Lehmann style bootstrap estimate.
		#' @param subject_or_block_weights Bootstrap weights at the subject/block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			beta = private$hl_point_estimate(private$y, private$w, row_weights)
			private$cached_values$beta_hat_T = beta
			if (!estimate_only) {
				y_vals = as.numeric(private$y); w_vals = as.integer(private$w); rw = as.numeric(row_weights)
				i_t = which(w_vals == 1L & is.finite(y_vals) & is.finite(rw) & rw > 0)
				i_c = which(w_vals == 0L & is.finite(y_vals) & is.finite(rw) & rw > 0)
				se = NA_real_
				if (length(i_t) >= 2L && length(i_c) >= 2L) {
					diffs = as.numeric(outer(y_vals[i_t], y_vals[i_c], "-"))
					wdiff = as.numeric(outer(rw[i_t], rw[i_c], "*"))
					ok = is.finite(diffs) & is.finite(wdiff) & wdiff > 0
					if (any(ok)) {
						diffs = diffs[ok]; wdiff = wdiff[ok]
						o = order(diffs); diffs = diffs[o]; wdiff = wdiff[o]
						cw = cumsum(wdiff) / sum(wdiff)
						q025 = diffs[which(cw >= 0.025)[1L]]
						q975 = diffs[which(cw >= 0.975)[1L]]
						if (is.finite(q025) && is.finite(q975) && q975 > q025)
							se = (q975 - q025) / (2 * 1.96)
					}
				}
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				private$cached_values$df = NA_real_
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
			}
			private$cached_values$beta_hat_T
		},
		#' @description Jackknife bias correction is unstable for the
		#'   Hodges-Lehmann estimator; report explicit non-estimability.
		compute_jackknife_estimate = function(unit = "auto"){
			private$cache_nonestimable_estimate("wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_corrected_estimate = function(unit = "auto"){
			self$compute_jackknife_estimate(unit = unit)
		},
		compute_jackknife_bias_estimate = function(unit = "auto"){
			private$cache_nonestimable_estimate("wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_std_error = function(unit = "auto"){
			private$cache_nonestimable_se("wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_standard_error = function(unit = "auto"){
			self$compute_jackknife_std_error(unit = unit)
		},
		compute_jackknife_wald_two_sided_pval = function(delta = 0, unit = "auto"){
			private$cache_nonestimable_se("wilcox_hl_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			private$cache_nonestimable_se("wilcox_hl_jackknife_not_supported")
			c(NA_real_, NA_real_)
		}
	),
	private = list(
		max_resample_attempts = 50L,
		hl_point_estimate = function(y_vals, w_vals, row_weights = NULL){
			if (is.null(row_weights)) {
				return(wilcox_hl_point_estimate_cpp(as.integer(w_vals), as.numeric(y_vals)))
			}
			# Weighted implementation
			y_vals = as.numeric(y_vals)
			w_vals = as.integer(w_vals)
			row_weights = as.numeric(row_weights)
			i_t = which(w_vals == 1L & is.finite(y_vals) & is.finite(row_weights) & row_weights > 0)
			i_c = which(w_vals == 0L & is.finite(y_vals) & is.finite(row_weights) & row_weights > 0)
			if (length(i_t) == 0L || length(i_c) == 0L) return(NA_real_)
			diffs = as.numeric(outer(y_vals[i_t], y_vals[i_c], "-"))
			wdiff = as.numeric(outer(row_weights[i_t], row_weights[i_c], "*"))
			ok = is.finite(diffs) & is.finite(wdiff) & wdiff > 0
			if (!any(ok)) return(NA_real_)
			diffs = diffs[ok]
			wdiff = wdiff[ok]
			o = order(diffs)
			diffs = diffs[o]
			wdiff = wdiff[o]
			cw = cumsum(wdiff) / sum(wdiff)
			idx = which(cw >= 0.5)[1L]
			if (!is.finite(idx) || is.na(idx)) return(NA_real_)
			as.numeric(diffs[idx])
		},
		compute_fast_bootstrap_distr = function(B, ...) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			if (private$is_KK) return(NULL)
			args = list(...)
			n = args[[1]]; y = args[[2]]; dead = args[[3]]; w = args[[4]]
			indices_mat = matrix(-1L, nrow = n, ncol = B)
			for (b in seq_len(B)) {
				attempt = 1L
				repeat {
					i_b = sample_int_replace_cpp(n, n)
					w_b = w[i_b]
					if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) {
						indices_mat[, b] = i_b - 1L
						break
					}
					attempt = attempt + 1L
					if (attempt > private$max_resample_attempts) break
				}
			}
			compute_wilcox_hl_distr_parallel_cpp(as.numeric(y), as.integer(w), matrix(as.integer(indices_mat), nrow=n), private$n_cpp_threads(B))
		},
		get_standard_error = function(){
			if (is.null(private$cached_values$s_beta_hat_T)) private$shared()
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			NA_real_
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			res = compute_wilcox_hl_distr_parallel_cpp(
				w_mat = as.matrix(w_mat),
				y = as.numeric(y),
				delta = as.numeric(delta),
				transform_code = 0L,
				zero_one_logit_clamp = as.numeric(zero_one_logit_clamp),
				num_cores = private$n_cpp_threads(ncol(w_mat))
			)
			return(res)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			yT = private$y[private$w == 1]
			yC = private$y[private$w == 0]
			if (length(yT) == 0L || length(yC) == 0L){
				private$cache_nonestimable_estimate("wilcox_empty_treatment_arm")
				return(invisible(NULL))
			}
			mod = tryCatch(stats::wilcox.test(yT, yC, conf.int = TRUE), error = function(e) NULL)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("wilcox_fit_unavailable")
				return(invisible(NULL))
			}
			beta = private$hl_point_estimate(private$y, private$w)
			ci   = mod$conf.int
			se   = if (length(ci) == 2L) (ci[2] - ci[1]) / (2 * 1.96) else NA_real_
			private$cached_values$beta_hat_T   = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
			private$cached_values$s_beta_hat_T = if (length(se)   == 1L && is.finite(se) && se > 0) se else NA_real_

			beta_hl = private$cached_values$beta_hat_T
			private$cached_values$likelihood_test_context = list(
				X = cbind(1, private$w),
				j = 2L,
				full_fit = list(b = c(mean(private$y - beta_hl * private$w), beta_hl), vt = var(yT), vc = var(yC))
			)
		},
		supports_lik_ratio_param_bootstrap = function() FALSE,
		supports_likelihood_tests = function() FALSE,
		get_supported_testing_types_impl = function(){
			"wald"
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$b)
			vt = spec$full_fit$vt; vc = spec$full_fit$vc
			w = spec$X[, 2]; n = length(w)
			y_sim = numeric(n)
			y_sim[w == 1] = b_null[1] + b_null[2] + rnorm(sum(w == 1), 0, sqrt(vt))
			y_sim[w == 0] = b_null[1] + rnorm(sum(w == 0), 0, sqrt(vc))
			
			hl_sim = private$hl_point_estimate(y_sim, w)
			list(
				full_fit = list(b = c(mean(y_sim - hl_sim * w), hl_sim), vt = var(y_sim[w==1]), vc = var(y_sim[w==0])),
				fit_null = function(d, start = NULL){
					m_joint = mean(y_sim - w * d)
					list(b = c(m_joint, d), vt = var(y_sim[w==1]), vc = var(y_sim[w==0]))
				},
				neg_loglik = function(fit){
					mu = fit$b[1] + fit$b[2]*w
					sum((y_sim - mu)^2)
				}
			)
		},
		get_likelihood_test_spec = function(){
			NULL
		}
	)
)
