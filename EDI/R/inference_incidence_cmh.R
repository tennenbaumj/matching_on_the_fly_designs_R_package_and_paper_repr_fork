#' CMH Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a randomization-based standard error.
#'
#' Legacy inference class. This class is retained for backwards compatibility
#' and is not comprehensively tested by the package comprehensive-test harness.
#'
#' @details
#' Treatment assignments are encoded as \eqn{w_i \in \{-1, +1\}}.  For a balanced
#' design the treatment-effect estimator is \eqn{\hat\tau = (2/n)\,\mathbf{y}'\mathbf{w}},
#' and since \eqn{E_w[\mathbf{y}'\mathbf{w}] = 0} for any balanced randomization the
#' standard error is
#' \deqn{SE(\hat\tau) = \frac{2}{n}\sqrt{\frac{\sum_k (\mathbf{y}'\mathbf{w}_k)^2}{K}}}{SE = (2/n) * sqrt(sum(ytw^2) / K)}
#' where \eqn{K} draws \eqn{\mathbf{w}_1,\ldots,\mathbf{w}_K} come from the design's
#' reference distribution.  Centering at the known zero mean (rather than the sample mean)
#' makes the denominator \eqn{K} rather than \eqn{K-1}.
#'
#' For blocking designs the expectation is evaluated exactly:
#' \deqn{SE(\hat\tau) = \frac{2}{n}\sqrt{\sum_b \frac{n_{1b}\,n_{0b}}{n_B - 1}}}{SE = (2/n) * sqrt(sum_b n1b*n0b / (nB - 1))}
#' where \eqn{n_{1b}, n_{0b}} are the numbers of positive and negative responses in
#' block \eqn{b} and \eqn{n_B} is the (common) block size.  This equals
#' \eqn{2\sqrt{V_{\rm CMH}}} where \eqn{V_{\rm CMH}} is the CMH variance from
#' Azriel et al. (2026), Equation 3.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneRandomBlockSize$new(n = 20, response_type = 'incidence', strata_cols = 'x1')
#' for (i in 1:20) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = factor(rep(1:2, 10)[i], levels=1:2)))
#' }
#' seq_des$add_all_subject_responses(rbinom(20, 1, 0.5))
#' inf = InferenceIncidCMH$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceIncidCMH = R6::R6Class("InferenceIncidCMH",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description Computes an approximate confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			self$compute_estimate()
			private$get_standard_error()
			super$compute_asymp_confidence_interval(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			self$compute_estimate()
			private$get_standard_error()
			super$compute_asymp_two_sided_pval(delta)
		},
		#' @description Initialize CMH incidence inference.
		#' @param des_obj A completed design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param se_est_num_vectors For non-block designs, the number of randomization vectors
		#'   drawn from the design to estimate the standard error. Default \code{1000L}.
		#' @param verbose Logical. Whether to print progress messages.
		#' @return A new \code{InferenceIncidCMH} object.
		initialize = function(des_obj, model_formula = NULL, se_est_num_vectors = 5000L, verbose = FALSE){
			if (des_obj$is_blocking_design()) {
				if (des_obj$get_prob_T() != 0.5) {
					stop("InferenceIncidCMH requires even treatment allocation for blocking designs.")
				}
				block_ids = des_obj$get_block_ids()
				block_sizes = as.integer(table(block_ids))
				if (length(block_sizes) > 1L && any(block_sizes != block_sizes[1L])) {
					stop("InferenceIncidCMH requires equal block sizes for blocking designs.")
				}
			}
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
				assertCount(se_est_num_vectors, positive = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$se_est_num_vectors = as.integer(se_est_num_vectors)
		}
	),
	private = list(
		se_est_num_vectors = NULL,
		supports_lik_ratio_param_bootstrap = function() FALSE,
		supports_likelihood_tests = function() FALSE,
		get_supported_testing_types_impl = function(){
			"wald"
		},
		get_standard_error = function(){
			if (!is.null(private$cached_values$cmh_s_beta_hat_T)) {
				se = private$cached_values$cmh_s_beta_hat_T
				if (is.finite(se) && se > 0) return(se)
				private$cache_nonestimable_se("cmh_standard_error_unavailable")
				return(NA_real_)
			}
			if (private$des_obj$is_blocking_design()) {
				private$cached_values$cmh_s_beta_hat_T = compute_cmh_block_se_cpp(
					private$des_obj_priv_int$y,
					private$des_obj$get_block_ids(),
					private$des_obj_priv_int$n
				)
			} else {
				precomp = private$des_obj$get_cmh_se_w_mat()
				w_mat = if (!is.null(precomp)) precomp else private$des_obj$draw_ws_according_to_design(private$se_est_num_vectors)
				ytw      = drop(private$y %*% w_mat)
				# With {-1,+1} encoding: τ̂ = (2/n)*y'w, so SE[τ̂] = (2/n)*SD[ytw].
				# E[y·w] = 0 exactly for all balanced designs, so the unbiased variance
				# estimator uses K (not K-1) in the denominator.
				K        = length(ytw)
				private$cached_values$cmh_s_beta_hat_T = 2 / private$n * sqrt(max(0, sum(ytw^2) / K))
			}
			if (!is.finite(private$cached_values$cmh_s_beta_hat_T) || private$cached_values$cmh_s_beta_hat_T <= 0) {
				private$cached_values$cmh_s_beta_hat_T = NA_real_
				private$cache_nonestimable_se("cmh_standard_error_unavailable")
				return(NA_real_)
			}
			private$cached_values$s_beta_hat_T = private$cached_values$cmh_s_beta_hat_T
			private$cached_values$df = NA_real_
			private$cached_values$cmh_s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			NA_real_
		}
	)
)
