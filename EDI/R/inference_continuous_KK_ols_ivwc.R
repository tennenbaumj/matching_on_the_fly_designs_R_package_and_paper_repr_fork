#' OLS IVWC Compound Inference for KK Designs
#'
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using OLS regression for matched-pair differences
#' and reservoir outcomes, with the treatment indicator and, optionally, all
#' recorded covariates as predictors.
#' Note that warm starts are disabled for this class as OLS is a closed-form
#' estimator and does not benefit from initialization.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinKKOLSIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#'
#' \strong{Legacy class.} Not fully tested in \code{comprehensive_tests.R}.
#' @export
InferenceContinKKOLSIVWC = R6::R6Class("InferenceContinKKOLSIVWC",
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
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertFormula(model_formula, null.ok = TRUE)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$fit_warm_start_enabled = FALSE
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	))),
	private = list(
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			preserve = if (is.null(permutations$m_mat)) c("kk_ols_ivwc_matched_reduced_design", "kk_ols_ivwc_reservoir_reduced_design") else character()
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, preserve_cache_keys = preserve, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		reduce_design_matrix_once = function(X, j_treat, cache_key){
			cached = private$cached_values[[cache_key]]
			if (!is.null(cached)) return(cached)
			qr_X = qr(X)
			if (qr_X$rank < ncol(X)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_treat %in% keep)) keep[qr_X$rank] = j_treat
				keep = sort(unique(keep))
				X = X[, keep, drop = FALSE]
				j_treat = which(keep == j_treat)
			}
			cached = list(X = X, j_treat = j_treat)
			private$cached_values[[cache_key]] = cached
			cached
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (m > 0){
				private$ols_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m  = private$cached_values$ssq_beta_T_matched
			df_m   = private$cached_values$df_beta_T_matched
			m_ok   = !is.null(beta_m) && is.finite(beta_m) && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0
			if (nRT > 0 && nRC > 0){
				private$ols_for_reservoir(estimate_only = estimate_only)
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r  = private$cached_values$ssq_beta_T_reservoir
			df_r   = private$cached_values$df_beta_T_reservoir
			r_ok   = !is.null(beta_r) && is.finite(beta_r) && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				if (!estimate_only) {
					private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
					private$cached_values$df = private$satterthwaite_df(
						var_terms = c(w_star^2 * ssq_m, (1 - w_star)^2 * ssq_r),
						dfs = c(df_m, df_r)
					)
				}
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				if (!estimate_only) {
					private$cached_values$s_beta_hat_T = sqrt(ssq_m)
					private$cached_values$df = df_m
				}
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				if (!estimate_only) {
					private$cached_values$s_beta_hat_T = sqrt(ssq_r)
					private$cached_values$df = df_r
				}
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				if (!estimate_only) {
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$df = NA_real_
				}
			}
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		satterthwaite_df = function(var_terms, dfs){
			ok = is.finite(var_terms) & var_terms > 0 & is.finite(dfs) & dfs > 0
			if (!any(ok)) return(NA_real_)
			var_terms = var_terms[ok]
			dfs = dfs[ok]
			total_var = sum(var_terms)
			if (!is.finite(total_var) || total_var <= 0) return(NA_real_)
			total_var^2 / sum((var_terms^2) / dfs)
		},
		fit_ols_with_treatment = function(X, y, j_treat, estimate_only = FALSE){
			if (nrow(X) <= ncol(X)) return(NULL)
			
			fit = tryCatch(stats::lm.fit(X, y), error = function(e) NULL)
			if (is.null(fit) || length(stats::coef(fit)) < j_treat || !is.finite(stats::coef(fit)[j_treat])){
				return(NULL)
			}
			beta = as.numeric(stats::coef(fit)[j_treat])
			df = nrow(X) - ncol(X)
			if (estimate_only) return(list(beta = beta, ssq = NA_real_, df = df))
			if (df <= 0) return(NULL)
			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X, as.numeric(y), as.numeric(stats::coef(fit)), j_treat),
				error = function(e) NULL
			)
			if (is.null(post_fit)) return(NULL)
			se = as.numeric(post_fit$std_err[j_treat])
			if (!is.finite(se) || se <= 0) return(NULL)
			list(beta = beta, ssq = se^2, df = df)
		},
		ols_for_matched_pairs = function(estimate_only = FALSE){
			yd = private$cached_values$KKstats$y_matched_diffs
			m  = length(yd)
			if (ncol(as.matrix(private$X)) > 0){
				Xd = as.matrix(private$cached_values$KKstats$X_matched_diffs)
				X = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				reduced = private$reduce_design_matrix_once(
					X,
					1L,
					cache_key = "kk_ols_ivwc_matched_reduced_design"
				)
				X = reduced$X
			} else {
				X = matrix(1, nrow = m, ncol = 1L)
			}
			fit = private$fit_ols_with_treatment(X, yd, 1L, estimate_only = estimate_only)
			if (is.null(fit)) {
				private$cached_values$beta_T_matched     = if (m >= 1) mean(yd) else NA_real_
				private$cached_values$ssq_beta_T_matched = if (m >= 2) var(yd) / m else NA_real_
				private$cached_values$df_beta_T_matched  = if (m >= 2) m - 1 else NA_real_
			} else {
				private$cached_values$beta_T_matched     = fit$beta
				private$cached_values$ssq_beta_T_matched = fit$ssq
				private$cached_values$df_beta_T_matched  = fit$df
			}
		},
		ols_for_reservoir = function(estimate_only = FALSE){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			X_r = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L
			if (ncol(as.matrix(private$X)) > 0){
				X_full = cbind(1, w_r, X_r)
				reduced = private$reduce_design_matrix_once(
					X_full,
					j_treat,
					cache_key = "kk_ols_ivwc_reservoir_reduced_design"
				)
				X_full = reduced$X
				j_treat = reduced$j_treat
			} else {
				X_full = cbind(1, w_r)
			}
			fit = private$fit_ols_with_treatment(X_full, y_r, j_treat, estimate_only = estimate_only)
			if (is.null(fit)) {
				private$cached_values$beta_T_reservoir     = NA_real_
				private$cached_values$ssq_beta_T_reservoir = NA_real_
				private$cached_values$df_beta_T_reservoir  = NA_real_
			} else {
				private$cached_values$beta_T_reservoir     = fit$beta
				private$cached_values$ssq_beta_T_reservoir = fit$ssq
				private$cached_values$df_beta_T_reservoir  = fit$df
			}
		}
	)
)
