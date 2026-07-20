.inflate_kk_onelik_standard_error_with_jackknife = function(private_env, self_obj){
	if (!is.null(private_env$active_resampling_operation)) return(invisible(NULL))
	theta_hat = as.numeric(private_env$cached_values$beta_hat_T %||% NA_real_)[1L]
	if (!is.finite(theta_hat)) return(invisible(NULL))
	se_model = as.numeric(private_env$cached_values$s_beta_hat_T %||% NA_real_)[1L]
	jack_summary = tryCatch(private_env$compute_jackknife_summary(unit = "auto"), error = function(e) NULL)
	se_jack = as.numeric(jack_summary$std_error %||% NA_real_)[1L]
	if (!is.finite(se_jack) || se_jack <= 0) return(invisible(NULL))
	if (!is.finite(se_model) || se_jack > se_model * (1 + sqrt(.Machine$double.eps))) {
		private_env$cached_values$s_beta_hat_T_model = se_model
		private_env$cached_values$s_beta_hat_T = se_jack
		private_env$cached_values$s_beta_hat_T_source = "jackknife"
	}
	invisible(NULL)
}

.conservative_kk_onelik_pval = function(model_p, design_p){
	vals = as.numeric(c(model_p, design_p))
	vals = vals[is.finite(vals)]
	if (!length(vals)) return(NA_real_)
	min(1, max(0, max(vals)))
}

.conservative_kk_onelik_ci = function(model_ci, design_ci, alpha = 0.05){
	model_ci = as.numeric(model_ci)
	design_ci = as.numeric(design_ci)
	ok_model = length(model_ci) >= 2L && all(is.finite(model_ci[1:2]))
	ok_design = length(design_ci) >= 2L && all(is.finite(design_ci[1:2]))
	ci = c(NA_real_, NA_real_)
	names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
	if (!ok_model && !ok_design) return(ci)
	if (!ok_model) {
		ci[] = design_ci[1:2]
		return(ci)
	}
	if (!ok_design) {
		ci[] = model_ci[1:2]
		return(ci)
	}
	ci[] = c(min(model_ci[1L], design_ci[1L]), max(model_ci[2L], design_ci[2L]))
	ci
}

#' KK Hurdle Poisson IVWC Inference for Count Responses
#'
#' Internal base class for KK hurdle-Poisson inverse-variance weighted combined
#' inference. The matched-pair component is fit with a hurdle-Poisson mixed model
#' using pair random intercepts, and the reservoir component is fit with an
#' ordinary Poisson log-link regression. The reported treatment effect is on the
#' log-rate scale.
#'
#' @keywords internal
#' @noRd
InferenceCountKKHurdlePoissonIVWC = R6::R6Class("InferenceCountKKHurdlePoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with count responses.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB} for
		#'   the matched-pair component.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(use_rcpp)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or DesignFixedBinaryMatch.")
				}
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (inherits(des_obj, "DesignFixedBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			private$use_rcpp = use_rcpp
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !private$use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
		},
		#' @description Compute treatment estimate
		#' @param estimate_only If \code{TRUE}, skip standard-error calculations.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute asymp confidence interval
		#' @param alpha Confidence level.
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute asymp two sided pval for treatment effect
		#' @param delta The null treatment effect (default 0).
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
					stop("Testing non-zero delta is not yet implemented for this class.")
				}
				NA_real_
			}
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			private$shared_combined_bootstrap(subject_or_block_weights, estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B Integer. Number of bootstrap samples (default 501).
		#' @param show_progress Logical. Whether to show a progress bar.
		#' @param debug Logical. Whether to return diagnostics.
		#' @param bootstrap_type Character. Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		use_rcpp = TRUE,
		max_abs_reasonable_coef = 1e4,
		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the fixed-effect coefficient for "w" directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			X = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(NA_integer_, nrow(X))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)
			beta_m = NA_real_
			ssq_m = NA_real_
			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(X, matched_idx, m_vec, se = FALSE)
				beta_m = res_m$beta_hat
				ssq_m = res_m$se^2
			}
			m_ok = !is.na(beta_m) && is.finite(beta_m) && !is.na(ssq_m) && is.finite(ssq_m) && ssq_m > 0
			beta_r = NA_real_
			ssq_r = NA_real_
			if (length(reservoir_idx) > 1L && length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(X, reservoir_idx, estimate_only = estimate_only)
				beta_r = res_r$beta_hat
				ssq_r = res_r$ssq_hat
			}
			r_ok = !is.na(beta_r) && is.finite(beta_r) && !is.na(ssq_r) && is.finite(ssq_r) && ssq_r > 0
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		},
		compute_basic_match_data = function(){
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv     = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		build_model_matrix = function(){
			if (ncol(as.matrix(private$X)) > 0){
				X = private$create_design_matrix()
				full_names = c("(Intercept)", "w", if (ncol(X) > 2L) paste0("x", seq_len(ncol(X) - 2L)) else NULL)
				colnames(X) = full_names[seq_len(ncol(X))]
			} else {
				X = cbind(1, private$w)
				colnames(X) = c("(Intercept)", "w")
			}
			X
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			X = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(NA_integer_, nrow(X))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)
			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(X, matched_idx, m_vec, se = !estimate_only)
				private$cached_values$beta_T_matched = res_m$beta_hat
				if (!estimate_only) private$cached_values$ssq_beta_T_matched = res_m$se^2 else private$cached_values$ssq_beta_T_matched = 1.0
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
				!is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0
			if (length(reservoir_idx) > 1L &&
				length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(X, reservoir_idx, estimate_only = estimate_only)
				private$cached_values$beta_T_reservoir = res_r$beta_hat
				if (!estimate_only) private$cached_values$ssq_beta_T_reservoir = res_r$ssq_hat else private$cached_values$ssq_beta_T_reservoir = 1.0
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) &&
				!is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cache_nonestimable_estimate("kk_hurdle_poisson_ivwc_no_usable_component")
			}
			invisible(NULL)
		},
		build_glmm_formula = function(dat){
			fixed_terms = setdiff(colnames(dat), c("y", "pair_group"))
			rhs = paste(c(fixed_terms, "(1 | pair_group)"), collapse = " + ")
			stats::as.formula(paste("y ~", rhs))
		},
		fit_hurdle_for_matched_pairs = function(X, matched_idx, m_vec, se = TRUE){
			X_matched = X[matched_idx, , drop = FALSE]
			if (is.null(dim(X_matched)) || ncol(X_matched) < 2L) {
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			reduced = private$reduce_design_matrix_preserving_treatment(X_matched)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			if (private$use_rcpp) {
				res = private$fit_hurdle_for_matched_pairs_rcpp(
					X_fit = X_fit,
					y_fit = private$y[matched_idx],
					group_id = m_vec[matched_idx],
					j_treat = reduced$j_treat,
					se = se
				)
				if (is.finite(res$beta_hat) && (!se || (is.finite(res$se) && res$se > 0))) {
					return(res)
				}
			}
			private$fit_hurdle_for_matched_pairs_glmm_tmb(
				X_fit = X_fit,
				y_fit = private$y[matched_idx],
				group_id = m_vec[matched_idx],
				se = se
			)
		},
		fit_hurdle_for_matched_pairs_rcpp = function(X_fit, y_fit, group_id, j_treat, se = TRUE){
			fit = tryCatch(
				fast_hurdle_poisson_glmm_cpp(
					X = as.matrix(X_fit),
					y = as.numeric(y_fit),
					group_id = as.integer(group_id),
					j_T = as.integer(j_treat - 1L),
					estimate_only = !se,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) return(list(beta_hat = NA_real_, se = NA_real_))
			beta_hat = as.numeric(fit$b[j_treat])
			if (!is.finite(beta_hat) || abs(beta_hat) > private$max_abs_reasonable_coef) {
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			if (!se) return(list(beta_hat = beta_hat, se = 1.0))
			ssq = as.numeric(fit$ssq_b_T)
			if (!is.finite(ssq) || ssq <= 0) return(list(beta_hat = NA_real_, se = NA_real_))
			se_val = sqrt(ssq)
			if (!is.finite(se_val) || se_val <= 0 || se_val > private$max_abs_reasonable_coef) {
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			list(beta_hat = beta_hat, se = se_val)
		},
		fit_hurdle_for_matched_pairs_glmm_tmb = function(X_fit, y_fit, group_id, se = TRUE){
			if (!check_package_installed("glmmTMB")) return(list(beta_hat = NA_real_, se = NA_real_))
			pred_df = as.data.frame(X_fit[, -1, drop = FALSE])
			colnames(pred_df)[1] = "w"
			dat = data.frame(
				y = y_fit,
				pred_df,
				pair_group = factor(group_id)
			)
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)
			formula_cond = private$build_glmm_formula(dat)
			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = stats::as.formula(sub("^y ~ ", "~ ", deparse(formula_cond))),
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat,
						control = glmm_control,
						se = se
					)
				)),
				error = function(e) NULL
			)
			if (is.null(mod) && ncol(dat) > 3L){
				dat = dat[, c("y", "w", "pair_group"), drop = FALSE]
				mod = tryCatch(
					suppressWarnings(suppressMessages(
						glmmTMB::glmmTMB(
							y ~ w + (1 | pair_group),
							ziformula = ~ w + (1 | pair_group),
							family = glmmTMB::truncated_poisson(link = "log"),
							data = dat,
							control = glmm_control,
							se = se
						)
					)),
					error = function(e) NULL
				)
			}
			if (is.null(mod)) return(list(beta_hat = NA_real_, se = NA_real_))
			if (!se){
				beta = glmmTMB::fixef(mod)$cond
				if ("w" %in% names(beta)){
					return(list(beta_hat = as.numeric(beta["w"]), se = 1.0)) # Return dummy SE > 0
				}
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table))) return(list(beta_hat = NA_real_, se = NA_real_))
			beta_hat = as.numeric(coef_table["w", "Estimate"])
			se_val = as.numeric(coef_table["w", "Std. Error"])
			if (!is.finite(beta_hat) || !is.finite(se_val) || se_val <= 0) return(list(beta_hat = NA_real_, se = NA_real_))
			list(beta_hat = beta_hat, se = se_val)
		},
		fit_poisson_for_reservoir = function(X, reservoir_idx, estimate_only = FALSE){
			X_res = X[reservoir_idx, , drop = FALSE]
			if (is.null(dim(X_res)) || ncol(X_res) < 2L) {
				return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			}
			reduced = private$reduce_design_matrix_preserving_treatment(X_res)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			}
			mod = tryCatch({
				if (estimate_only) {
					fast_poisson_regression_cpp(X_fit, private$y[reservoir_idx])
				} else {
					fast_poisson_regression_with_var_cpp(X_fit, private$y[reservoir_idx], j = reduced$j_treat)
				}
			}, error = function(e) NULL)
			if (is.null(mod) || !isTRUE(mod$converged)) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			beta_hat = as.numeric(mod$b[reduced$j_treat])
			if (estimate_only) {
				if (!is.finite(beta_hat)) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
				return(list(beta_hat = beta_hat, ssq_hat = 1))
			}
			ssq_hat = as.numeric(mod$ssq_b_j)
			if (!is.finite(beta_hat) || !is.finite(ssq_hat) || ssq_hat <= 0) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			list(beta_hat = beta_hat, ssq_hat = ssq_hat)
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)))
)
#' KK Hurdle Poisson Combined-Likelihood Inference for Count Responses
#'
#' @export
InferenceCountKKHurdlePoissonOneLik = R6::R6Class("InferenceCountKKHurdlePoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with count responses.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB}.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param optimization_alg Optimization algorithm.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(use_rcpp)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$use_rcpp = use_rcpp
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !private$use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
			if (private$has_match_structure) {
				if (inherits(des_obj, "DesignFixedBinaryMatch")) {
					des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
				}
				private$m = des_obj$.__enclos_env__$private$m
				private$compute_basic_match_data()
			}
		},
		#' @description Compute treatment estimate
		#' @param estimate_only If \code{TRUE}, skip standard-error calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_hurdle(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the combined hurdle-Poisson estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			private$compute_weighted_combined_hurdle_estimate(row_weights, estimate_only = estimate_only)
		},
		#' @description Compute asymp confidence interval
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			switch(
				self$get_testing_type(),
				wald = {
					private$shared_combined_hurdle()
					if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
						return(self$compute_bootstrap_confidence_interval(alpha = alpha))
					}
					private$compute_z_or_t_ci_from_s_and_df(alpha)
				},
				score = self$compute_score_confidence_interval(alpha = alpha),
				lik_ratio = self$compute_lik_ratio_confidence_interval(alpha = alpha),
				gradient = self$compute_gradient_confidence_interval(alpha = alpha)
			)
		},
		#' @description Computes a design-conservative score confidence interval.
		#' @param alpha Significance level.
		compute_score_confidence_interval = function(alpha = 0.05){
			private$shared_combined_hurdle()
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			ci_design = private$compute_z_or_t_ci_from_s_and_df(alpha)
			ci_model = tryCatch(super$compute_score_confidence_interval(alpha = alpha), error = function(e) c(NA_real_, NA_real_))
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (all(is.finite(ci_design[1:2]))) private$clear_nonestimable_state()
			.conservative_kk_onelik_ci(ci_model, ci_design, alpha)
		},
		#' @description Computes a design-conservative likelihood-ratio confidence interval.
		#' @param alpha Significance level.
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			private$shared_combined_hurdle()
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			ci_design = private$compute_z_or_t_ci_from_s_and_df(alpha)
			ci_model = tryCatch(super$compute_lik_ratio_confidence_interval(alpha = alpha), error = function(e) c(NA_real_, NA_real_))
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (all(is.finite(ci_design[1:2]))) private$clear_nonestimable_state()
			.conservative_kk_onelik_ci(ci_model, ci_design, alpha)
		},
		#' @description Computes a design-conservative gradient confidence interval.
		#' @param alpha Significance level.
		compute_gradient_confidence_interval = function(alpha = 0.05){
			private$shared_combined_hurdle()
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			ci_design = private$compute_z_or_t_ci_from_s_and_df(alpha)
			ci_model = tryCatch(super$compute_gradient_confidence_interval(alpha = alpha), error = function(e) c(NA_real_, NA_real_))
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (all(is.finite(ci_design[1:2]))) private$clear_nonestimable_state()
			.conservative_kk_onelik_ci(ci_model, ci_design, alpha)
		},
		#' @description Compute asymp two sided pval
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			switch(
				self$get_testing_type(),
				wald = {
					private$shared_combined_hurdle()
					if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
						return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
					}
					private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
				},
				score = self$compute_score_two_sided_pval(delta = delta),
				lik_ratio = self$compute_lik_ratio_two_sided_pval(delta = delta),
				gradient = self$compute_gradient_two_sided_pval(delta = delta)
			)
		},
		#' @description Computes a design-conservative score p-value.
		#' @param delta Null treatment effect value.
		compute_score_two_sided_pval = function(delta = 0){
			private$shared_combined_hurdle()
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			p_design = private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			p_model = tryCatch(super$compute_score_two_sided_pval(delta = delta), error = function(e) NA_real_)
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (is.finite(p_design)) private$clear_nonestimable_state()
			.conservative_kk_onelik_pval(p_model, p_design)
		},
		#' @description Computes a design-conservative likelihood-ratio p-value.
		#' @param delta Null treatment effect value.
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			private$shared_combined_hurdle()
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			p_design = private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			p_model = tryCatch(super$compute_lik_ratio_two_sided_pval(delta = delta), error = function(e) NA_real_)
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (is.finite(p_design)) private$clear_nonestimable_state()
			.conservative_kk_onelik_pval(p_model, p_design)
		},
		#' @description Computes a design-conservative gradient p-value.
		#' @param delta Null treatment effect value.
		compute_gradient_two_sided_pval = function(delta = 0){
			private$shared_combined_hurdle()
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			p_design = private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			p_model = tryCatch(super$compute_gradient_two_sided_pval(delta = delta), error = function(e) NA_real_)
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (is.finite(p_design)) private$clear_nonestimable_state()
			.conservative_kk_onelik_pval(p_model, p_design)
		},
		#' @description Compute asymp confidence interval
		#' @param alpha Confidence level.
		compute_wald_confidence_interval = function(alpha = 0.05){
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute Wald two sided pval
		#' @param delta Null treatment effect value.
		compute_wald_two_sided_pval = function(delta = 0){
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B Integer. Number of bootstrap samples (default 501).
		#' @param show_progress Logical. Whether to show a progress bar.
		#' @param debug Logical. Whether to return diagnostics.
		#' @param bootstrap_type Character. Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		},
		supports_lik_ratio_param_bootstrap = function() TRUE
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		m = NULL,
		cached_mod = NULL,
		use_rcpp = TRUE,
		max_abs_reasonable_coef = 1e4,
		supports_likelihood_tests = function(){
			TRUE
		},
		get_supported_testing_types_impl = function(){
			c("wald", "score", "gradient", "lik_ratio")
		},
		warn_bootstrap_fallback_once = function(){
			if (!isTRUE(private$cached_values$warned_bootstrap_se_unavailable)) {
				private$cached_values$warned_bootstrap_se_unavailable = TRUE
				warning("KK hurdle-Poisson combined-likelihood: falling back to bootstrap because standard error is unavailable.")
			}
		},
		compute_basic_match_data = function(){
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},
		build_model_matrix = function(){
			if (ncol(as.matrix(private$X)) > 0){
				X = private$create_design_matrix()
				full_names = c("(Intercept)", "w", if (ncol(X) > 2L) paste0("x", seq_len(ncol(X) - 2L)) else NULL)
				colnames(X) = full_names[seq_len(ncol(X))]
			} else {
				X = cbind(1, private$w)
				colnames(X) = c("(Intercept)", "w")
			}
			as.matrix(X)
		},
		build_combined_hurdle_data = function(X_fit, j_treat){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, nrow(X_fit))
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)
			list(
				X_fit = X_fit,
				matched_idx = matched_idx,
				reservoir_idx = reservoir_idx,
				X_matched = X_fit[matched_idx, , drop = FALSE],
				y_matched = as.numeric(private$y[matched_idx]),
				group_id = as.integer(m_vec[matched_idx]),
				X_reservoir = X_fit[reservoir_idx, , drop = FALSE],
				y_reservoir = as.numeric(private$y[reservoir_idx]),
				j_treat = as.integer(j_treat)
			)
		},
		build_weighted_combined_hurdle_data = function(X_fit, j_treat, row_weights = NULL){
			dat = private$build_combined_hurdle_data(X_fit, j_treat)
			if (is.null(row_weights)) row_weights = rep(1, nrow(X_fit))
			row_weights = as.numeric(row_weights)
			row_weights[!is.finite(row_weights) | row_weights < 0] = 0
			dat$weights_matched = row_weights[dat$matched_idx]
			dat$weights_reservoir = row_weights[dat$reservoir_idx]
			dat
		},
		combined_hurdle_neg_loglik = function(params, dat){
			p = ncol(dat$X_fit)
			beta = as.numeric(params[seq_len(p)])
			total = 0
			w_m = dat$weights_matched %||% rep(1, length(dat$y_matched))
			if (length(dat$matched_idx) > 0L && sum(dat$y_matched > 0 & w_m > 0, na.rm = TRUE) > p) {
				total = total + if (is.null(dat$weights_matched)) {
					as.numeric(get_hurdle_poisson_glmm_neg_loglik_cpp(
						X = dat$X_matched,
						y = dat$y_matched,
						group_id = dat$group_id,
						params = as.numeric(params)
					))
				} else {
					as.numeric(get_hurdle_poisson_glmm_weighted_neg_loglik_cpp(
						X = dat$X_matched,
						y = dat$y_matched,
						group_id = dat$group_id,
						weights = w_m,
						params = as.numeric(params)
					))
				}
			}
			if (length(dat$reservoir_idx) > 0L) {
				w_r = dat$weights_reservoir %||% rep(1, length(dat$y_reservoir))
				eta_r = as.numeric(dat$X_reservoir %*% beta)
				total = total - sum(w_r * (dat$y_reservoir * eta_r - exp(pmin(eta_r, 20)) - lgamma(dat$y_reservoir + 1)))
			}
			as.numeric(total)
		},
		combined_hurdle_score = function(params, dat){
			p = ncol(dat$X_fit)
			score = numeric(p + 1L)
			w_m = dat$weights_matched %||% rep(1, length(dat$y_matched))
			if (length(dat$matched_idx) > 0L && sum(dat$y_matched > 0 & w_m > 0, na.rm = TRUE) > p) {
				score = score + if (is.null(dat$weights_matched)) {
					as.numeric(get_hurdle_poisson_glmm_score_cpp(
						X = dat$X_matched,
						y = dat$y_matched,
						group_id = dat$group_id,
						params = as.numeric(params)
					))
				} else {
					as.numeric(get_hurdle_poisson_glmm_weighted_score_cpp(
						X = dat$X_matched,
						y = dat$y_matched,
						group_id = dat$group_id,
						weights = w_m,
						params = as.numeric(params)
					))
				}
			}
			if (length(dat$reservoir_idx) > 0L) {
				w_r = dat$weights_reservoir %||% rep(1, length(dat$y_reservoir))
				score[seq_len(p)] = score[seq_len(p)] + as.numeric(get_poisson_regression_weighted_score_cpp(
					X = dat$X_reservoir,
					y = dat$y_reservoir,
					weights = w_r,
					beta = as.numeric(params[seq_len(p)])
				))
			}
			score
		},
		combined_hurdle_hessian = function(params, dat){
			p = ncol(dat$X_fit)
			H = matrix(0, nrow = p + 1L, ncol = p + 1L)
			w_m = dat$weights_matched %||% rep(1, length(dat$y_matched))
			if (length(dat$matched_idx) > 0L && sum(dat$y_matched > 0 & w_m > 0, na.rm = TRUE) > p) {
				H = H + if (is.null(dat$weights_matched)) {
					as.matrix(get_hurdle_poisson_glmm_hessian_cpp(
						X = dat$X_matched,
						y = dat$y_matched,
						group_id = dat$group_id,
						params = as.numeric(params)
					))
				} else {
					as.matrix(get_hurdle_poisson_glmm_weighted_hessian_cpp(
						X = dat$X_matched,
						y = dat$y_matched,
						group_id = dat$group_id,
						weights = w_m,
						params = as.numeric(params)
					))
				}
			}
			if (length(dat$reservoir_idx) > 0L) {
				w_r = dat$weights_reservoir %||% rep(1, length(dat$y_reservoir))
				H[seq_len(p), seq_len(p)] = H[seq_len(p), seq_len(p), drop = FALSE] + as.matrix(get_poisson_regression_weighted_hessian_cpp(
					X = dat$X_reservoir,
					weights = w_r,
					beta = as.numeric(params[seq_len(p)])
				))
			}
			H
		},
		information_inverse_diagonal_entry = function(information, j){
			info = as.matrix(information)
			if (!is.matrix(info) || nrow(info) != ncol(info) || any(!is.finite(info))) return(NA_real_)
			j = as.integer(j)
			if (length(j) != 1L || !is.finite(j) || j < 1L || j > nrow(info)) return(NA_real_)
			info = (info + t(info)) / 2
			vcov = tryCatch(solve(info), error = function(e) NULL)
			if (!is.null(vcov) && all(is.finite(vcov))) {
				vcov = (vcov + t(vcov)) / 2
				val = as.numeric(vcov[j, j])
				if (is.finite(val) && val >= 0) return(val)
			}

			eig = tryCatch(eigen(info, symmetric = TRUE), error = function(e) NULL)
			if (is.null(eig) || any(!is.finite(eig$values)) || any(!is.finite(eig$vectors))) return(NA_real_)
			scale = max(abs(eig$values), 1)
			tol = scale * sqrt(.Machine$double.eps)
			if (any(eig$values < -tol)) return(NA_real_)
			pos = eig$values > tol
			if (!any(pos)) return(NA_real_)

			e_j = numeric(nrow(info))
			e_j[j] = 1
			V_pos = eig$vectors[, pos, drop = FALSE]
			resid = e_j - as.numeric(V_pos %*% crossprod(V_pos, e_j))
			if (sqrt(sum(resid^2)) > 100 * sqrt(.Machine$double.eps)) return(NA_real_)

			val = sum((eig$vectors[j, pos]^2) / eig$values[pos])
			if (is.finite(val) && val >= 0) val else NA_real_
		},
		record_combined_hurdle_fit_summary = function(fit, X_model, X_fit, fallback_used = FALSE, fallback_reason = NULL){
			model_names = colnames(X_model)
			fit_names = colnames(X_fit)
			if (is.null(model_names)) model_names = paste0("V", seq_len(ncol(X_model)))
			if (is.null(fit_names)) fit_names = paste0("V", seq_len(ncol(X_fit)))
			beta_fit = as.numeric(fit$b %||% fit$params[seq_along(fit_names)])
			if (length(beta_fit) < length(fit_names)) return(invisible(NULL))
			beta_full = rep(NA_real_, length(model_names))
			names(beta_full) = model_names
			beta_full[match(fit_names, model_names)] = beta_fit[seq_along(fit_names)]
			
			se_full = rep(NA_real_, length(model_names))
			names(se_full) = model_names
			vcov_fit = tryCatch({
				info = as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information)
				solve(info)
			}, error = function(e) NULL)
			se_log_sigma = NA_real_
			if (!is.null(vcov_fit) && all(is.finite(vcov_fit))) {
				diag_v = diag(vcov_fit)
				if (length(diag_v) >= length(fit_names)) {
					se_fit = sqrt(pmax(0, diag_v[seq_along(fit_names)]))
					se_full[match(fit_names, model_names)] = se_fit
				}
				if (length(diag_v) >= length(fit_names) + 1L) {
					se_log_sigma = sqrt(max(0, as.numeric(diag_v[length(fit_names) + 1L])))
				}
			}
			
			log_sigma = as.numeric(fit$log_sigma %||% fit$params[length(fit$params)] %||% NA_real_)[1L]
			coef_all = c(beta_full, log_sigma = log_sigma)
			se_all = c(se_full, log_sigma = se_log_sigma)
			summary_table = matrix(NA_real_, nrow = length(coef_all), ncol = 4L)
			rownames(summary_table) = names(coef_all)
			colnames(summary_table) = c("Value", "Std. Error", "z value", "Pr(>|z|)")
			summary_table[, 1L] = coef_all
			summary_table[, 2L] = se_all
			ok = is.finite(coef_all) & is.finite(se_all) & se_all > 0
			summary_table[ok, 3L] = coef_all[ok] / se_all[ok]
			summary_table[ok, 4L] = 2 * stats::pnorm(-abs(summary_table[ok, 3L]))
			private$cached_values$full_coefficients = beta_full
			private$cached_values$summary_table = summary_table
			private$cached_values$model_fit_fallback = NULL
			if (isTRUE(fallback_used)) {
				private$cached_values$model_fit_fallback = list(
					used = TRUE,
					reason = fallback_reason %||% "full_kk_hurdle_poisson_onelik_fit_failed_to_converge",
					requested_model = "full covariate-adjusted KK combined hurdle-Poisson likelihood",
					fitted_model = "treatment-only KK combined hurdle-Poisson likelihood",
					omitted_conditional = setdiff(model_names, fit_names),
					family = "KK combined hurdle-Poisson"
				)
			}
			invisible(NULL)
		},
		fit_combined_hurdle = function(dat, estimate_only = FALSE, fixed_idx = NULL, fixed_values = NULL, warm_start_params = NULL){
			p = ncol(dat$X_fit)
			n_params = p + 1L
			start = as.numeric(warm_start_params %||% private$get_fit_warm_start_for_length("params", n_params))
			if (length(start) != n_params) {
				start = rep(0, n_params)
				start[n_params] = -3
				if (length(dat$reservoir_idx) > 0L) {
					w_r_start = dat$weights_reservoir %||% rep(1, length(dat$y_reservoir))
					pois_fit = tryCatch(
						fast_poisson_regression_weighted_cpp(
							X = dat$X_reservoir,
							y = dat$y_reservoir,
							weights = w_r_start,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (!is.null(pois_fit) && length(pois_fit$b) == p && all(is.finite(pois_fit$b))) {
						start[seq_len(p)] = as.numeric(pois_fit$b)
					}
				}
				if (length(dat$matched_idx) > 0L && sum(dat$y_matched > 0, na.rm = TRUE) > p) {
					match_fit = tryCatch(
						fast_hurdle_poisson_glmm_cpp(
							X = dat$X_matched,
							y = dat$y_matched,
							group_id = dat$group_id,
							j_T = as.integer(dat$j_treat - 1L),
							estimate_only = TRUE,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (!is.null(match_fit) && isTRUE(match_fit$converged) &&
					    length(match_fit$params) == n_params && all(is.finite(match_fit$params))) {
						start = as.numeric(match_fit$params)
					}
				}
			}
			if (is.null(fixed_idx)) fixed_idx = integer(0)
			if (is.null(fixed_values)) fixed_values = numeric(0)
			free_idx = setdiff(seq_len(n_params), fixed_idx)
			par_to_full = function(par_free){
				par = start
				if (length(fixed_idx) > 0L) par[fixed_idx] = fixed_values
				par[free_idx] = par_free
				par
			}
			obj = function(par_free){
				par = par_to_full(par_free)
				private$combined_hurdle_neg_loglik(par, dat)
			}
			gr = function(par_free){
				par = par_to_full(par_free)
				-private$combined_hurdle_score(par, dat)[free_idx]
			}
			opt = tryCatch(
				stats::optim(
					par = start[free_idx],
					fn = obj,
					gr = gr,
					method = if (identical(private$optimization_alg %||% "lbfgs", "newton_raphson")) "BFGS" else "BFGS",
					control = list(maxit = 300, reltol = 1e-8)
				),
				error = function(e) NULL
			)
			if (is.null(opt)) return(NULL)
			params = par_to_full(opt$par)
			score = private$combined_hurdle_score(params, dat)
			hessian = private$combined_hurdle_hessian(params, dat)
			info = -hessian
			ssq_b_T = NA_real_
			if (!estimate_only) {
				info_free = info[free_idx, free_idx, drop = FALSE]
				j_free = match(dat$j_treat, free_idx)
				ssq_b_T = private$information_inverse_diagonal_entry(info_free, j_free)
			}
			list(
				params = params,
				b = params[seq_len(p)],
				log_sigma = params[n_params],
				ssq_b_T = ssq_b_T,
				score = score,
				observed_information = info,
				information = info,
				fisher_information = info,
				hessian = hessian,
				neg_loglik = private$combined_hurdle_neg_loglik(params, dat),
				neg_ll = private$combined_hurdle_neg_loglik(params, dat),
				converged = is.finite(opt$value) && isTRUE(opt$convergence == 0L)
			)
		},
		shared_combined_hurdle = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$cached_values$likelihood_test_context = NULL
			private$cached_values$model_fit_fallback = NULL
			private$cached_mod = NULL
			X_full = private$build_model_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			if (!is.null(X_fit)) colnames(X_fit) = colnames(X_full)[reduced$keep]
			if (is.null(X_fit) || !is.finite(reduced$j_treat)) {
				private$cache_nonestimable_estimate("kk_hurdle_poisson_onelik_design_unusable")
				return(invisible(NULL))
			}
			X_model = X_fit
			dat = private$build_combined_hurdle_data(X_fit, reduced$j_treat)
			fit = private$fit_combined_hurdle(dat, estimate_only = estimate_only)
			fallback_used = FALSE
			fallback_reason = NULL
			if ((is.null(fit) || !isTRUE(fit$converged) || length(fit$b) < dat$j_treat || !is.finite(fit$b[dat$j_treat])) &&
			    ncol(X_fit) > 2L) {
				keep = sort(unique(c(1L, reduced$j_treat)))
				if (length(keep) == 2L) {
					X_fit_fallback = X_fit[, keep, drop = FALSE]
					j_treat_fallback = match(colnames(X_fit)[reduced$j_treat], colnames(X_fit_fallback))
					if (length(j_treat_fallback) == 1L && is.finite(j_treat_fallback)) {
						dat_fallback = private$build_combined_hurdle_data(X_fit_fallback, j_treat_fallback)
						fit_fallback = private$fit_combined_hurdle(dat_fallback, estimate_only = estimate_only)
						if (!is.null(fit_fallback) && isTRUE(fit_fallback$converged) &&
						    length(fit_fallback$b) >= dat_fallback$j_treat &&
						    is.finite(fit_fallback$b[dat_fallback$j_treat])) {
							X_fit = X_fit_fallback
							dat = dat_fallback
							fit = fit_fallback
							fallback_used = TRUE
							fallback_reason = "full_kk_hurdle_poisson_onelik_fit_failed_to_converge"
						}
					}
				}
			}
			if (is.null(fit) || !isTRUE(fit$converged) || length(fit$b) < dat$j_treat || !is.finite(fit$b[dat$j_treat])) {
				private$cache_nonestimable_estimate("kk_hurdle_poisson_onelik_fit_failed")
				return(invisible(NULL))
			}
			private$cached_mod = fit
			private$cached_values$beta_hat_T = as.numeric(fit$b[dat$j_treat])
			private$cached_values$df = Inf
			private$set_fit_warm_start(as.numeric(fit$params), "params", fisher = fit$fisher_information)
			private$cached_values$likelihood_test_context = dat
			private$record_combined_hurdle_fit_summary(
				fit = fit,
				X_model = X_model,
				X_fit = X_fit,
				fallback_used = fallback_used,
				fallback_reason = fallback_reason
			)
			if (!estimate_only) {
				se = sqrt(as.numeric(fit$ssq_b_T))
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef) se else NA_real_
				.inflate_kk_onelik_standard_error_with_jackknife(private, self)
			}
			private$clear_nonestimable_state()
			invisible(NULL)
		},
		compute_weighted_combined_hurdle_estimate = function(row_weights, estimate_only = FALSE){
			X_full = private$build_model_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			if (!is.null(X_fit)) colnames(X_fit) = colnames(X_full)[reduced$keep]
			if (is.null(X_fit) || !is.finite(reduced$j_treat)) {
				private$cache_nonestimable_estimate("kk_hurdle_poisson_onelik_weighted_design_unusable")
				return(NA_real_)
			}
			row_weights = as.numeric(row_weights)
			if (length(row_weights) != nrow(X_fit)) {
				private$cache_nonestimable_estimate("kk_hurdle_poisson_onelik_weighted_length_mismatch")
				return(NA_real_)
			}
			if (weights_are_effectively_constant(row_weights)) {
				return(as.numeric(self$compute_estimate(estimate_only = estimate_only))[1L])
			}
			dat = private$build_weighted_combined_hurdle_data(X_fit, reduced$j_treat, row_weights = row_weights)
			fit = private$fit_combined_hurdle(dat, estimate_only = estimate_only)
			if ((is.null(fit) || !isTRUE(fit$converged) || length(fit$b) < dat$j_treat || !is.finite(fit$b[dat$j_treat])) &&
			    ncol(X_fit) > 2L) {
				keep = sort(unique(c(1L, reduced$j_treat)))
				if (length(keep) == 2L) {
					X_fit_fallback = X_fit[, keep, drop = FALSE]
					j_treat_fallback = match(colnames(X_fit)[reduced$j_treat], colnames(X_fit_fallback))
					if (length(j_treat_fallback) == 1L && is.finite(j_treat_fallback)) {
						dat_fallback = private$build_weighted_combined_hurdle_data(X_fit_fallback, j_treat_fallback, row_weights = row_weights)
						fit_fallback = private$fit_combined_hurdle(dat_fallback, estimate_only = estimate_only)
						if (!is.null(fit_fallback) && isTRUE(fit_fallback$converged) &&
						    length(fit_fallback$b) >= dat_fallback$j_treat &&
						    is.finite(fit_fallback$b[dat_fallback$j_treat])) {
							dat = dat_fallback
							fit = fit_fallback
						}
					}
				}
			}
			if (is.null(fit) || !isTRUE(fit$converged) || length(fit$b) < dat$j_treat || !is.finite(fit$b[dat$j_treat])) {
				private$cache_nonestimable_estimate("kk_hurdle_poisson_onelik_weighted_fit_failed")
				return(NA_real_)
			}
			private$cached_mod = fit
			private$cached_values$beta_hat_T = as.numeric(fit$b[dat$j_treat])
			private$cached_values$df = Inf
			if (!estimate_only) {
				se = sqrt(as.numeric(fit$ssq_b_T))
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef) se else NA_real_
			}
			private$clear_nonestimable_state()
			private$cached_values$beta_hat_T
		},
		shared = function(estimate_only = FALSE){
			private$shared_combined_hurdle(estimate_only = estimate_only)
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			dat        = spec$dat
			params_null = as.numeric(null_fit$params)
			p          = ncol(dat$X_fit)
			beta       = params_null[seq_len(p)]
			log_sigma  = params_null[p + 1L]
			sigma      = exp(min(log_sigma, 8))
			if (!is.finite(sigma) || sigma < 0) return(NULL)

			y_matched_sim = dat$y_matched
			if (length(dat$matched_idx) > 0L) {
				G      = max(as.integer(dat$group_id))
				u_g    = rnorm(G, 0, sigma)
				eta_m  = as.numeric(dat$X_matched %*% beta) + u_g[as.integer(dat$group_id)]
				lam_m  = exp(pmin(eta_m, 20))
				pi_pos = plogis(eta_m)
				is_pos = as.logical(rbinom(length(eta_m), 1L, pi_pos))
				y_matched_sim = integer(length(eta_m))
				for (i in which(is_pos)){
					u = runif(1, exp(-lam_m[i]), 1)
					y_matched_sim[i] = max(qpois(u, lam_m[i]), 1L)
				}
			}

			y_reservoir_sim = dat$y_reservoir
			if (length(dat$reservoir_idx) > 0L) {
				lam_r           = exp(pmin(as.numeric(dat$X_reservoir %*% beta), 20))
				y_reservoir_sim = rpois(length(lam_r), lam_r)
			}

			dat_sim             = dat
			dat_sim$y_matched   = as.numeric(y_matched_sim)
			dat_sim$y_reservoir = as.numeric(y_reservoir_sim)

			j        = spec$j
			full_res = tryCatch(
				private$fit_combined_hurdle(dat_sim, estimate_only = FALSE),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged)) return(NULL)
			if (length(full_res$b) < j || !is.finite(full_res$b[j])) return(NULL)
			list(
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					private$fit_combined_hurdle(
						dat             = dat_sim,
						estimate_only   = FALSE,
						fixed_idx       = j,
						fixed_values    = d,
						warm_start_params = start %||% as.numeric(full_res$params)
					)
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% private$combined_hurdle_neg_loglik(as.numeric(fit$params), dat_sim))
				}
			)
		},
		get_likelihood_test_spec = function(){
			private$shared_combined_hurdle(estimate_only = FALSE)
			dat = private$cached_values$likelihood_test_context
			if (is.null(dat) || is.null(private$cached_mod)) return(NULL)
			j_treat = as.integer(dat$j_treat)
			list(
				dat = dat,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					private$fit_combined_hurdle(
						dat = dat,
						estimate_only = FALSE,
						fixed_idx = j_treat,
						fixed_values = delta,
						warm_start_params = start %||% as.numeric(private$cached_mod$params)
					)
				},
				extract_start = function(fit){
					as.numeric(fit$params)
				},
				score = function(fit){
					as.numeric(fit$score %||% private$combined_hurdle_score(as.numeric(fit$params), dat))
				},
				observed_information = function(fit){
					as.matrix(fit$observed_information %||% -private$combined_hurdle_hessian(as.numeric(fit$params), dat))
				},
				fisher_information = function(fit){
					as.matrix(fit$fisher_information %||% fit$observed_information %||% -private$combined_hurdle_hessian(as.numeric(fit$params), dat))
				},
				information = function(fit){
					as.matrix(fit$information %||% fit$fisher_information %||% fit$observed_information %||% -private$combined_hurdle_hessian(as.numeric(fit$params), dat))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% private$combined_hurdle_neg_loglik(as.numeric(fit$params), dat))
				}
			)
		}
	)))
)
#' Conditional-Poisson Inference for KK Designs with Combined Likelihood
#' @export
InferenceCountKKCondPoissonOneLik = R6::R6Class("InferenceCountKKCondPoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with count responses.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_cpoisson(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the combined conditional-Poisson estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$beta_hat_T = beta_hat_T
					private$cached_values$s_beta_hat_T = NA_real_
					return(private$cached_values$beta_hat_T)
				}
			}
			private$cached_values$beta_hat_T = private$compute_weighted_combined_estimate(row_weights)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			switch(
				self$get_testing_type(),
				wald = self$compute_wald_confidence_interval(alpha = alpha),
				score = self$compute_score_confidence_interval(alpha = alpha),
				lik_ratio = self$compute_lik_ratio_confidence_interval(alpha = alpha),
				gradient = self$compute_gradient_confidence_interval(alpha = alpha)
			)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			switch(
				self$get_testing_type(),
				wald = self$compute_wald_two_sided_pval(delta = delta),
				score = self$compute_score_two_sided_pval(delta = delta),
				lik_ratio = self$compute_lik_ratio_two_sided_pval(delta = delta),
				gradient = self$compute_gradient_two_sided_pval(delta = delta)
			)
		},
		#' @description Computes a Wald confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_wald_confidence_interval = function(alpha = 0.05){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes a Wald two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_wald_two_sided_pval = function(delta = 0){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes a design-adjusted score confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_score_confidence_interval = function(alpha = 0.05){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			ci_design = private$compute_z_or_t_ci_from_s_and_df(alpha)
			ci_model = tryCatch(super$compute_score_confidence_interval(alpha = alpha), error = function(e) c(NA_real_, NA_real_))
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (all(is.finite(ci_design[1:2]))) private$clear_nonestimable_state()
			.conservative_kk_onelik_ci(ci_model, ci_design, alpha)
		},
		#' @description Computes a design-adjusted likelihood-ratio confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			ci_design = private$compute_z_or_t_ci_from_s_and_df(alpha)
			ci_model = tryCatch(super$compute_lik_ratio_confidence_interval(alpha = alpha), error = function(e) c(NA_real_, NA_real_))
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (all(is.finite(ci_design[1:2]))) private$clear_nonestimable_state()
			.conservative_kk_onelik_ci(ci_model, ci_design, alpha)
		},
		#' @description Computes a design-adjusted gradient confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_gradient_confidence_interval = function(alpha = 0.05){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			ci_design = private$compute_z_or_t_ci_from_s_and_df(alpha)
			ci_model = tryCatch(super$compute_gradient_confidence_interval(alpha = alpha), error = function(e) c(NA_real_, NA_real_))
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (all(is.finite(ci_design[1:2]))) private$clear_nonestimable_state()
			.conservative_kk_onelik_ci(ci_model, ci_design, alpha)
		},
		#' @description Computes a design-adjusted score p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_score_two_sided_pval = function(delta = 0){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			p_design = private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			p_model = tryCatch(super$compute_score_two_sided_pval(delta = delta), error = function(e) NA_real_)
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (is.finite(p_design)) private$clear_nonestimable_state()
			.conservative_kk_onelik_pval(p_model, p_design)
		},
		#' @description Computes a design-adjusted likelihood-ratio p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			p_design = private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			p_model = tryCatch(super$compute_lik_ratio_two_sided_pval(delta = delta), error = function(e) NA_real_)
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (is.finite(p_design)) private$clear_nonestimable_state()
			.conservative_kk_onelik_pval(p_model, p_design)
		},
		#' @description Computes a design-adjusted gradient p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_gradient_two_sided_pval = function(delta = 0){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			beta_design = private$cached_values$beta_hat_T
			se_design = private$cached_values$s_beta_hat_T
			df_design = private$cached_values$df
			p_design = private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			p_model = tryCatch(super$compute_gradient_two_sided_pval(delta = delta), error = function(e) NA_real_)
			private$cached_values$beta_hat_T = beta_design
			private$cached_values$s_beta_hat_T = se_design
			private$cached_values$df = df_design
			if (is.finite(p_design)) private$clear_nonestimable_state()
			.conservative_kk_onelik_pval(p_model, p_design)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B Integer. Number of bootstrap samples (default 501).
		#' @param show_progress Logical. Whether to show a progress bar.
		#' @param debug Logical. Whether to return diagnostics.
		#' @param bootstrap_type Character. Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		},
		supports_lik_ratio_param_bootstrap = function() TRUE
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		cached_mod = NULL,
		max_abs_reasonable_coef = 1e4,
		get_supported_testing_types_impl = function(){
			c("wald", "score", "gradient", "lik_ratio")
		},
		compute_basic_match_data = function(){
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},
		build_model_matrix = function(){
			if (ncol(as.matrix(private$X)) > 0){
				X = private$create_design_matrix()
				full_names = c("(Intercept)", "w", if (ncol(X) > 2L) paste0("x", seq_len(ncol(X) - 2L)) else NULL)
				colnames(X) = full_names[seq_len(ncol(X))]
			} else {
				X = cbind(1, private$w)
				colnames(X) = c("(Intercept)", "w")
			}
			as.matrix(X)
		},
		build_combined_cpoisson_data = function(X_fit, j_treat){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, nrow(X_fit))
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			j_treat = as.integer(j_treat)
			cov_cols = setdiff(seq_len(ncol(X_fit)), c(1L, j_treat))
			pair_ids = sort(unique(m_vec[m_vec > 0L]))
			n_max = length(pair_ids)
			yT = numeric(n_max)
			nk = numeric(n_max)
			Xdiff_rows = if (length(cov_cols) > 0L) vector("list", n_max) else NULL
			pair_rows = vector("list", n_max)
			keep_pair = 0L
			for (pid in pair_ids){
				rows = which(m_vec == pid)
				if (length(rows) != 2L) next
				t_row = rows[private$w[rows] == 1]
				c_row = rows[private$w[rows] == 0]
				if (length(t_row) != 1L || length(c_row) != 1L) next
				keep_pair = keep_pair + 1L
				pair_rows[[keep_pair]] = c(t_row, c_row)
				yT[keep_pair] = as.numeric(private$y[t_row])
				nk[keep_pair] = as.numeric(private$y[t_row] + private$y[c_row])
				if (length(cov_cols) > 0L) {
					Xdiff_rows[[keep_pair]] = as.numeric(X_fit[t_row, cov_cols] - X_fit[c_row, cov_cols])
				}
			}
			if (keep_pair < n_max) {
				yT = yT[seq_len(keep_pair)]
				nk = nk[seq_len(keep_pair)]
				pair_rows = pair_rows[seq_len(keep_pair)]
				if (!is.null(Xdiff_rows)) Xdiff_rows = Xdiff_rows[seq_len(keep_pair)]
			}
			Xdiff = if (length(cov_cols) > 0L && keep_pair > 0L) {
				matrix(unlist(Xdiff_rows), nrow = keep_pair, ncol = length(cov_cols), byrow = TRUE)
			} else {
				matrix(numeric(0), nrow = 0L, ncol = length(cov_cols))
			}
			reservoir_idx = which(m_vec <= 0L)
			X_res = if (length(cov_cols) > 0L) X_fit[reservoir_idx, cov_cols, drop = FALSE] else matrix(nrow = length(reservoir_idx), ncol = 0L)
			list(
				yT_v = as.numeric(yT),
				n_k_v = as.numeric(nk),
				X_diff_v = as.matrix(Xdiff),
				y_r = as.numeric(private$y[reservoir_idx]),
				w_r = as.numeric(private$w[reservoir_idx]),
				X_r = as.matrix(X_res),
				j_treat = 2L,
				pair_rows = pair_rows,
				reservoir_idx = reservoir_idx
			)
		},
		set_failed_combined_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cache_nonestimable_estimate("kk_cpoisson_onelik_fit_failed")
		},
		reduce_combined_covariates = function(X_diff, X_r, w_r){
			n_p = nrow(as.matrix(X_diff))
			n_r = length(w_r)
			p = ncol(as.matrix(X_diff))
			pairs_part = if (n_p > 0L) {
				cbind(Intercept = rep(0, n_p), treatment = rep(1, n_p), as.matrix(X_diff))
			} else {
				matrix(0, nrow = 0, ncol = p + 2L)
			}
			reservoir_part = if (n_r > 0L) {
				cbind(Intercept = rep(1, n_r), treatment = w_r, as.matrix(X_r))
			} else {
				matrix(0, nrow = 0, ncol = p + 2L)
			}
			X_stack = rbind(pairs_part, reservoir_part)
			if (nrow(X_stack) == 0L) return(integer(0))
			qr_res = qr(X_stack)
			if (is.finite(qr_res$rank) && qr_res$rank < ncol(X_stack)){
				keep = qr_res$pivot[seq_len(qr_res$rank)]
				if (!(1L %in% keep)) keep = c(1L, keep)
				if (!(2L %in% keep)) keep = c(2L, keep)
				keep = sort(unique(keep))
				return(keep[keep > 2L] - 2L)
			}
			seq_len(p)
		},
		fit_combined_cpoisson = function(dat, estimate_only = FALSE, fixed_idx = NULL, fixed_values = NULL, warm_start_params = NULL){
			n_params = ncol(dat$X_diff_v) + 2L
			warm_fisher = private$get_fit_warm_start_fisher(n_params)
			tryCatch(
				fast_cpoisson_combined_with_var_cpp(
					yT_v = dat$yT_v,
					n_k_v = dat$n_k_v,
					X_diff_v = dat$X_diff_v,
					y_r = dat$y_r,
					w_r = dat$w_r,
					X_r = dat$X_r,
					fixed_idx = fixed_idx,
					fixed_values = fixed_values,
					warm_start_params = warm_start_params %||% private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = warm_fisher
				),
				error = function(e) NULL
			)
		},
		try_combined_fit = function(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v){
			n_params = ncol(X_diff_v) + 2L
			mod = tryCatch(
				fast_cpoisson_combined_with_var_cpp(
					yT_v = as.numeric(yT_v),
					n_k_v = as.numeric(n_k_v),
					X_diff_v = as.matrix(X_diff_v),
					y_r = as.numeric(y_r_v),
					w_r = as.numeric(w_r_v),
					X_r = as.matrix(X_r_v),
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
				),
				error = function(e) NULL
			)
			if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2L])) return(FALSE)
			if (!estimate_only) {
				ssq = mod$ssq_b_j
				if (is.null(ssq) || !is.finite(ssq) || ssq < 0) return(FALSE)
			}
			private$cached_mod = mod
			private$set_fit_warm_start(as.numeric(mod$params), "params", fisher = mod$fisher_information)
			private$cached_values$likelihood_test_context = list(
				yT_v = as.numeric(yT_v),
				n_k_v = as.numeric(n_k_v),
				X_diff_v = as.matrix(X_diff_v),
				y_r_v = as.numeric(y_r_v),
				w_r_v = as.numeric(w_r_v),
				X_r_v = as.matrix(X_r_v)
			)
			private$cached_values$beta_hat_T = as.numeric(mod$b[2L])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$df = NA_real_
			TRUE
		},
		try_pairs_only = function(estimate_only, yT_v, n_k_v, X_diff_v){
			if (length(yT_v) == 0L) return(FALSE)
			y_prop = yT_v / n_k_v
			X = if (ncol(X_diff_v) > 0L) cbind(1, X_diff_v) else matrix(1, nrow = length(yT_v), ncol = 1L)
			mod = tryCatch({
				res = fast_logistic_regression_weighted_cpp(X = X, y = y_prop, weights = n_k_v)
				list(b = res$b, ssq_b_1 = NA_real_, X_fit = X)
			}, error = function(e) NULL)
			if (is.null(mod) || length(mod$b) < 1L || !is.finite(mod$b[1L])) return(FALSE)
			private$cached_values$beta_hat_T = as.numeric(mod$b[1L])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			TRUE
		},
		try_reservoir_only = function(estimate_only, y_r_v, w_r_v, X_r_v){
			if (length(y_r_v) == 0L) return(FALSE)
			X_full = if (ncol(X_r_v) > 0L) cbind(1, w_r_v, X_r_v) else cbind(1, w_r_v)
			X_fit = X_full
			j_treat = 2L
			if (ncol(X_full) > 2L){
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (is.finite(r_full) && r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(2L %in% keep)) keep = c(2L, keep)
					keep = sort(unique(keep))
					X_fit = X_full[, keep, drop = FALSE]
					j_match = which(colnames(X_fit) == "w_r_v")
					if (length(j_match) > 0L) j_treat = j_match[1L]
				}
			}
			mod = tryCatch(
				fast_poisson_regression_with_var_cpp(X = X_fit, y = as.numeric(y_r_v), j = as.integer(j_treat)),
				error = function(e) NULL
			)
			if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
			if (!estimate_only) {
				ssq = mod$ssq_b_j
				if (is.null(ssq) || !is.finite(ssq) || ssq < 0) return(FALSE)
			}
			private$cached_values$beta_hat_T = as.numeric(mod$b[j_treat])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$df = NA_real_
			TRUE
		},
		weighted_cpoisson_neg_loglik = function(params, dat){
			p = ncol(dat$X_diff_v)
			beta_0 = params[1]
			beta_T = params[2]
			beta_x = if (p > 0L) params[seq.int(3L, p + 2L)] else numeric(0)
			ll = 0
			if (length(dat$pair_weights) > 0L) {
				eta_p = beta_T + if (p > 0L) as.numeric(dat$X_diff_v %*% beta_x) else 0
				ll = ll + sum(dat$pair_weights * (dat$yT_v * eta_p - dat$n_k_v * ifelse(eta_p > 0, eta_p + log1p(exp(-eta_p)), log1p(exp(eta_p)))))
			}
			if (length(dat$reservoir_weights) > 0L) {
				eta_r = beta_0 + beta_T * dat$w_r + if (p > 0L) as.numeric(dat$X_r %*% beta_x) else 0
				ll = ll + sum(dat$reservoir_weights * (dat$y_r * eta_r - exp(pmin(eta_r, 20)) - lgamma(dat$y_r + 1)))
			}
			-as.numeric(ll)
		},
		weighted_cpoisson_score = function(params, dat){
			p = ncol(dat$X_diff_v)
			score = numeric(p + 2L)
			beta_0 = params[1]
			beta_T = params[2]
			beta_x = if (p > 0L) params[seq.int(3L, p + 2L)] else numeric(0)
			if (length(dat$pair_weights) > 0L) {
				eta_p = beta_T + if (p > 0L) as.numeric(dat$X_diff_v %*% beta_x) else 0
				prob = plogis(eta_p)
				resid = dat$pair_weights * (dat$yT_v - dat$n_k_v * prob)
				score[2] = score[2] + sum(resid)
				if (p > 0L) score[seq.int(3L, p + 2L)] = score[seq.int(3L, p + 2L)] + as.numeric(crossprod(dat$X_diff_v, resid))
			}
			if (length(dat$reservoir_weights) > 0L) {
				eta_r = beta_0 + beta_T * dat$w_r + if (p > 0L) as.numeric(dat$X_r %*% beta_x) else 0
				mu_r = exp(pmin(eta_r, 20))
				resid = dat$reservoir_weights * (dat$y_r - mu_r)
				score[1] = score[1] + sum(resid)
				score[2] = score[2] + sum(resid * dat$w_r)
				if (p > 0L) score[seq.int(3L, p + 2L)] = score[seq.int(3L, p + 2L)] + as.numeric(crossprod(dat$X_r, resid))
			}
			score
		},
		compute_weighted_combined_estimate = function(row_weights){
			X_full = private$build_model_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			if (!is.null(X_fit)) colnames(X_fit) = colnames(X_full)[reduced$keep]
			if (is.null(X_fit) || !is.finite(reduced$j_treat)) return(NA_real_)
			dat = private$build_combined_cpoisson_data(X_fit, reduced$j_treat)
			dat$pair_weights = if (length(dat$pair_rows) > 0L) vapply(dat$pair_rows, function(rows) mean(row_weights[rows]), numeric(1)) else numeric(0)
			dat$reservoir_weights = as.numeric(row_weights[dat$reservoir_idx])
			n_params = ncol(dat$X_diff_v) + 2L
			start = private$get_fit_warm_start_for_length("params", n_params)
			if (length(start) != n_params) {
				full_fit = private$fit_combined_cpoisson(dat, estimate_only = TRUE)
				start = if (!is.null(full_fit) && length(full_fit$params) == n_params) as.numeric(full_fit$params) else rep(0, n_params)
			}
			opt = tryCatch(
				stats::optim(
					par = start,
					fn = function(par) private$weighted_cpoisson_neg_loglik(par, dat),
					gr = function(par) -private$weighted_cpoisson_score(par, dat),
					method = "BFGS",
					control = list(maxit = 300, reltol = 1e-8)
				),
				error = function(e) NULL
			)
			if (is.null(opt) || length(opt$par) < 2L || !is.finite(opt$par[2L])) return(NA_real_)
			private$set_fit_warm_start(as.numeric(opt$par), "params")
			as.numeric(opt$par[2L])
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		shared_combined_cpoisson = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$cached_values$likelihood_test_context = NULL
			private$cached_mod = NULL
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
				private$cached_values$combined_cov_keep = NULL
			}
			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			p = if (is.null(private$X)) 0L else ncol(as.matrix(private$X))
			has_reservoir = nRT > 0 && nRC > 0
			yT_v = numeric(0)
			n_k_v = numeric(0)
			X_diff_v = matrix(nrow = 0L, ncol = p)
			if (m > 0){
				yT = KKstats$yTs_matched
				yC = KKstats$yCs_matched
				n_k = yT + yC
				valid = which(n_k > 0)
				if (length(valid) > 0L){
					yT_v = yT[valid]
					n_k_v = n_k[valid]
					if (p > 0L) X_diff_v = as.matrix(KKstats$X_matched_diffs_full[valid, , drop = FALSE])
				}
			}
			has_pairs = length(yT_v) > 0L
			y_r_v = numeric(0)
			w_r_v = numeric(0)
			X_r_v = matrix(nrow = 0L, ncol = p)
			if (has_reservoir){
				y_r_v = KKstats$y_reservoir
				w_r_v = KKstats$w_reservoir
				if (p > 0L) X_r_v = as.matrix(KKstats$X_reservoir)
			}
			if (!has_pairs && !has_reservoir){
				private$set_failed_combined_cache()
				return(invisible(NULL))
			}
			if (p > 0L && is.null(private$cached_values$combined_cov_keep)){
				private$cached_values$combined_cov_keep = private$reduce_combined_covariates(X_diff_v, X_r_v, w_r_v)
			}
			if (p > 0L){
				keep = private$cached_values$combined_cov_keep
				if (length(keep) > 0L){
					X_diff_v = X_diff_v[, keep, drop = FALSE]
					X_r_v = X_r_v[, keep, drop = FALSE]
				} else {
					X_diff_v = matrix(nrow = nrow(X_diff_v), ncol = 0L)
					X_r_v = matrix(nrow = nrow(X_r_v), ncol = 0L)
				}
			}
			success = FALSE
			if (has_pairs && has_reservoir){
				success = private$try_combined_fit(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v)
			}
			if (!success){
				fallback_success = FALSE
				if (has_pairs){
					fallback_success = private$try_pairs_only(estimate_only, yT_v, n_k_v, X_diff_v)
				}
				if (!fallback_success && has_reservoir){
					fallback_success = private$try_reservoir_only(estimate_only, y_r_v, w_r_v, X_r_v)
				}
				if (!fallback_success){
					private$set_failed_combined_cache()
				} else {
					private$clear_nonestimable_state()
				}
			} else {
				private$clear_nonestimable_state()
			}
			if (!estimate_only) {
				.inflate_kk_onelik_standard_error_with_jackknife(private, self)
			}
			invisible(NULL)
		},
		shared = function(estimate_only = FALSE){
			private$shared_combined_cpoisson(estimate_only = estimate_only)
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			dat         = spec$dat
			params_null = as.numeric(null_fit$params %||% null_fit$b)
			p           = ncol(dat$X_diff_v)
			beta_0      = params_null[1L]
			beta_T      = params_null[2L]
			beta_x      = if (p > 0L) params_null[seq.int(3L, p + 2L)] else numeric(0)
			j           = spec$j

			yT_sim = dat$yT_v
			if (length(dat$yT_v) > 0L) {
				eta_pair = beta_T + if (p > 0L) as.numeric(dat$X_diff_v %*% beta_x) else 0
				prob     = plogis(eta_pair)
				yT_sim   = as.numeric(rbinom(length(dat$n_k_v), as.integer(dat$n_k_v), prob))
			}

			y_r_sim = dat$y_r
			if (length(dat$y_r) > 0L) {
				eta_r   = beta_0 + beta_T * dat$w_r + if (p > 0L) as.numeric(dat$X_r %*% beta_x) else 0
				y_r_sim = as.numeric(rpois(length(dat$y_r), exp(pmin(eta_r, 20))))
			}

			dat_sim       = dat
			dat_sim$yT_v  = yT_sim
			dat_sim$y_r   = y_r_sim

			full_res = tryCatch(
				private$fit_combined_cpoisson(dat_sim, estimate_only = FALSE),
				error = function(e) NULL
			)
			if (is.null(full_res) || length(full_res$b) < 2L || !is.finite(full_res$b[2L])) return(NULL)
			list(
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					private$fit_combined_cpoisson(
						dat           = dat_sim,
						estimate_only = FALSE,
						fixed_idx     = j,
						fixed_values  = d,
						warm_start_params = start %||% as.numeric(full_res$params)
					)
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik)
				}
			)
		},
		get_likelihood_test_spec = function(){
			private$shared_combined_cpoisson(estimate_only = FALSE)
			dat = private$cached_values$likelihood_test_context
			if (is.null(dat) || is.null(private$cached_mod)) return(NULL)
			list(
				dat = dat,
				j = as.integer(dat$j_treat),
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					private$fit_combined_cpoisson(
						dat = dat,
						estimate_only = FALSE,
						fixed_idx = dat$j_treat,
						fixed_values = delta,
						warm_start_params = start
					)
				},
				extract_start = function(fit){
					as.numeric(fit$params %||% fit$b)
				},
				score = function(fit){
					params = as.numeric(fit$params %||% fit$b)
					as.numeric(fit$score %||% get_cpoisson_combined_score_cpp(dat$yT_v, dat$n_k_v, dat$X_diff_v, dat$y_r, dat$w_r, dat$X_r, params))
				},
				observed_information = function(fit){
					params = as.numeric(fit$params %||% fit$b)
					as.matrix(fit$observed_information %||% -get_cpoisson_combined_hessian_cpp(dat$yT_v, dat$n_k_v, dat$X_diff_v, dat$y_r, dat$w_r, dat$X_r, params))
				},
				fisher_information = function(fit){
					params = as.numeric(fit$params %||% fit$b)
					as.matrix(fit$fisher_information %||% fit$observed_information %||% -get_cpoisson_combined_hessian_cpp(dat$yT_v, dat$n_k_v, dat$X_diff_v, dat$y_r, dat$w_r, dat$X_r, params))
				},
				information = function(fit){
					params = as.numeric(fit$params %||% fit$b)
					as.matrix(fit$information %||% fit$fisher_information %||% fit$observed_information %||% -get_cpoisson_combined_hessian_cpp(dat$yT_v, dat$n_k_v, dat$X_diff_v, dat$y_r, dat$w_r, dat$X_r, params))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik)
				}
			)
		}
	)))
)
