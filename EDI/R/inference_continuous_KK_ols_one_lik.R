#' OLS Combined-Likelihood Inference for KK Designs
#'
#' Fits a single stacked OLS regression over matched-pair differences and
#' reservoir observations for KK matching-on-the-fly designs with continuous
#' responses, using the treatment indicator and, optionally, all recorded
#' covariates as predictors.
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
#' inf = InferenceContinKKOLSOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinKKOLSOneLik = R6::R6Class("InferenceContinKKOLSOneLik",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompound,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertFormula(model_formula, null.ok = TRUE)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula)
			private$fit_warm_start_enabled = FALSE
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Returns the combined OLS estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$fit_combined(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute the treatment estimate with bootstrap weights.
		#' @param subject_or_block_weights Numeric vector. Row weights for bootstrap.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		#' @return The treatment estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				self$compute_estimate(estimate_only = estimate_only)
				beta_hat_T = as.numeric(private$cached_values$beta_hat_T)[1L]
				if (is.finite(beta_hat_T)) {
					return(private$cached_values$beta_hat_T)
				}
			}
			result = private$fit_weighted_combined(row_weights, estimate_only = estimate_only)
			if (is.list(result)) {
				private$cached_values$beta_hat_T = result$beta
				private$cached_values$s_beta_hat_T = result$se
			} else {
				private$cached_values$beta_hat_T = result
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$beta_hat_T
		},
		#' @description Computes the approximate confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 -
		#'   \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			switch(
				private$testing_type,
				wald = {
					private$fit_combined()
					if (should_run_asserts()) {
						private$assert_finite_se()
					}
					private$compute_z_or_t_ci_from_s_and_df(alpha)
				},
				score = private$invert_test_pval_confidence_interval(alpha),
				gradient = private$invert_gradient_ci_uniroot(alpha),
				lik_ratio = private$invert_lik_ratio_ci_newton(alpha)
			)
		},
		#' @description Computes the approximate p-value.
		#' @param delta The null difference to test against. For any treatment effect at all this
		#'   is set to zero (the default).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			switch(
				private$testing_type,
				wald = {
					private$fit_combined()
					if (should_run_asserts()) {
						private$assert_finite_se()
					}
					private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
				},
				score = private$compute_score_two_sided_pval_impl(delta),
				gradient = private$compute_gradient_two_sided_pval_impl(delta),
				lik_ratio = private$compute_lik_ratio_two_sided_pval_impl(delta)
			)
		},
		#' @description Returns the log-likelihood, gradient, and Hessian at the current estimate.
		#' @return A list with \code{loglik}, \code{gradient}, and \code{hessian}.
		get_likelihood_components = function(){
			private$fit_combined()
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) return(NULL)
			list(
				loglik = -spec$neg_loglik(spec$full_fit),
				gradient = spec$score(spec$full_fit),
				hessian = -spec$observed_information(spec$full_fit)
			)
		}
	),
	private = list(
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			preserve = if (is.null(permutations$m_mat)) c("kk_ols_combined_reduced_design") else character()
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp, preserve_cache_keys = preserve)
		},
		get_standard_error = function(){
			private$fit_combined(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T %||% NA_real_
		},
		get_degrees_of_freedom = function(){
			private$fit_combined(estimate_only = FALSE)
			private$cached_values$df %||% NA_real_
		},
		supports_likelihood_tests = function() TRUE,
		supports_lik_ratio_param_bootstrap = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null   = as.numeric(null_fit$b)
			sig2     = spec$full_fit$sigma2_hat
			n        = nrow(spec$X)
			X_fit    = spec$X
			j        = spec$j
			mu       = as.numeric(X_fit %*% b_null)
			y_sim    = private$simulate_param_boot_gaussian_y(mu, sig2)
			if (is.null(y_sim)) return(NULL)
			if (n <= ncol(X_fit)) return(NULL)
			lm_boot  = tryCatch(stats::lm.fit(X_fit, y_sim), error = function(e) NULL)
			if (is.null(lm_boot)) return(NULL)
			b_boot   = as.numeric(stats::coef(lm_boot))
			if (length(b_boot) < j || !is.finite(b_boot[j])) return(NULL)
			full_fit_boot = list(b = b_boot)
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					tryCatch(
						fast_ols_with_var_cpp(
							X = X_fit, y = y_sim, j = j,
							fixed_idx = j, fixed_values = d
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit){
					rss = sum((y_sim - X_fit %*% as.numeric(fit$b))^2)
					0.5 * rss / sig2
				}
			)
		},
		get_supported_testing_types_impl = function() c("wald", "score", "gradient", "lik_ratio"),
		get_score_test_information_matrix = function(spec, fit){
			tryCatch(spec$fisher_information(fit), error = function(e) NULL)
		},
		compute_score_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		get_likelihood_test_spec = function(){
			private$fit_combined()
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = ctx$y
			j_treat = as.integer(ctx$j_treat)
			sig2 = private$cached_mod$sigma2_hat
			
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					fast_ols_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						fixed_idx = j_treat,
						fixed_values = delta
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					as.numeric(t(X_fit) %*% (y - X_fit %*% as.numeric(fit$b)) / sig2)
				},
				observed_information = function(fit){
					(t(X_fit) %*% X_fit) / sig2
				},
				fisher_information = function(fit){
					(t(X_fit) %*% X_fit) / sig2
				},
				information = function(fit){
					(t(X_fit) %*% X_fit) / sig2
				},
				neg_loglik = function(fit){
					rss = sum((y - X_fit %*% as.numeric(fit$b))^2)
					# For the purposes of testing beta, we can treat sigma2 as fixed from the full model
					# to stay consistent with the standard score/LR test behavior in this package's GLM infrastructure.
					0.5 * rss / sig2
				}
			)
		},
		# Copied from InferenceAsympLikStdModCache to avoid multiple inheritance issues 
		# while still using the shared infrastructure.
		compute_likelihood_test_two_sided_pval = function(delta, testing_type){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) return(NA_real_)
			private$get_memoized_likelihood_test_pval(
				delta = delta,
				testing_type = testing_type,
				spec = spec,
				warm_cache_key = paste0("likelihood_test:", testing_type)
			)
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
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		fit_ols = function(X, y, j_treat, estimate_only = FALSE){
			if (nrow(X) <= ncol(X)) return(NULL)
			fit = tryCatch(stats::lm.fit(X, y), error = function(e) NULL)
			if (is.null(fit) || length(stats::coef(fit)) < j_treat || !is.finite(stats::coef(fit)[j_treat])){
				return(NULL)
			}
			beta_vec = as.numeric(stats::coef(fit))
			beta = beta_vec[j_treat]
			df = nrow(X) - ncol(X)
			if (estimate_only) return(list(b = beta_vec, beta = beta, se = NA_real_, mod = NULL, df = df))
			if (df <= 0) return(NULL)
			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X, as.numeric(y), beta_vec, j_treat),
				error = function(e) NULL
			)
			if (is.null(post_fit)) return(NULL)
			se = as.numeric(post_fit$std_err[j_treat])
			if (!is.finite(se) || se <= 0) return(NULL)
			res = stats::residuals(fit)
			rss = sum(res^2)
			sig2 = rss / df
			list(
				b = beta_vec,
				beta = beta,
				se = se,
				mod = fit,
				sigma2_hat = sig2,
				df = df,
				vcov = post_fit$vcov
			)
		},
		fit_combined = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			nR  = nRT + nRC
			if (m > 0 && nRT > 0 && nRC > 0){
				if (ncol(as.matrix(private$X)) > 0){
					Xd_full = as.matrix(KKstats$X_matched_diffs_full)
					p = ncol(Xd_full)
				} else {
					Xd_full = matrix(nrow = m, ncol = 0L)
					p = 0L
				}
				X_comb = rbind(
					if (p > 0L) cbind(0, 1, Xd_full) else cbind(numeric(m), rep(1.0, m)),
					if (p > 0L) cbind(rep(1, nR), KKstats$w_reservoir, as.matrix(KKstats$X_reservoir)) else cbind(rep(1, nR), KKstats$w_reservoir)
				)
				y_comb = c(KKstats$y_matched_diffs, KKstats$y_reservoir)
				j_treat = 2L
			} else if (m > 0){
				Xd = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_matched_diffs) else matrix(nrow = m, ncol = 0L)
				X_comb = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				y_comb = KKstats$y_matched_diffs
				j_treat = 1L
			} else if (nRT > 0 && nRC > 0){
				X_r = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_reservoir) else matrix(nrow = nR, ncol = 0L)
				X_comb = if (ncol(X_r) > 0L) cbind(rep(1, nR), KKstats$w_reservoir, X_r) else cbind(rep(1, nR), KKstats$w_reservoir)
				y_comb = KKstats$y_reservoir
				j_treat = 2L
			} else {
				private$cache_nonestimable_estimate("no_usable_matched_or_reservoir_data")
				return(invisible(NULL))
			}
			reduced = private$reduce_design_matrix_once(
				X_comb,
				j_treat,
				cache_key = "kk_ols_combined_reduced_design"
			)
			X_comb = reduced$X
			j_treat = reduced$j_treat
			fit = private$fit_ols(X_comb, y_comb, j_treat, estimate_only = estimate_only)
			if (is.null(fit)){
				private$cache_nonestimable_estimate("ols_fit_unavailable")
				private$cached_values$likelihood_test_context = NULL
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T   = fit$beta
			if (!estimate_only) {
				private$cached_values$s_beta_hat_T = fit$se
				private$cached_values$df = fit$df
				private$cached_mod = fit
				private$cached_values$likelihood_test_context = list(
					X = X_comb,
					y = y_comb,
					j_treat = j_treat
				)
				if (!is.null(fit$mod)) {
					private$cached_values$full_coefficients = stats::coef(fit$mod)
					private$cached_values$full_vcov = fit$vcov
				}
			}
			invisible(NULL)
		},
		fit_weighted_combined = function(row_weights, estimate_only = TRUE){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}
			if (is.null(KKstats)) return(NA_real_)
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			nR  = nRT + nRC
			kk_w = kk_pair_and_reservoir_bootstrap_weights(private, row_weights)
			if (m > 0 && nRT > 0 && nRC > 0){
				if (ncol(as.matrix(private$X)) > 0){
					Xd_full = as.matrix(KKstats$X_matched_diffs_full)
					p = ncol(Xd_full)
				} else {
					Xd_full = matrix(nrow = m, ncol = 0L)
					p = 0L
				}
				X_comb = rbind(
					if (p > 0L) cbind(0, 1, Xd_full) else cbind(numeric(m), rep(1.0, m)),
					if (p > 0L) cbind(rep(1, nR), KKstats$w_reservoir, as.matrix(KKstats$X_reservoir)) else cbind(rep(1, nR), KKstats$w_reservoir)
				)
				y_comb = c(KKstats$y_matched_diffs, KKstats$y_reservoir)
				w_comb = c(kk_w$pair_weights, kk_w$reservoir_weights)
				j_treat = 2L
			} else if (m > 0){
				Xd = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_matched_diffs) else matrix(nrow = m, ncol = 0L)
				X_comb = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				y_comb = KKstats$y_matched_diffs
				w_comb = kk_w$pair_weights
				j_treat = 1L
			} else if (nRT > 0 && nRC > 0){
				X_r = if (ncol(as.matrix(private$X)) > 0) as.matrix(KKstats$X_reservoir) else matrix(nrow = nR, ncol = 0L)
				X_comb = if (ncol(X_r) > 0L) cbind(rep(1, nR), KKstats$w_reservoir, X_r) else cbind(rep(1, nR), KKstats$w_reservoir)
				y_comb = KKstats$y_reservoir
				w_comb = kk_w$reservoir_weights
				j_treat = 2L
			} else {
				return(NA_real_)
			}
			ok = is.finite(w_comb) & w_comb > 0 & is.finite(y_comb)
			if (!any(ok)) return(NA_real_)
			X_comb = X_comb[ok, , drop = FALSE]
			y_comb = y_comb[ok]
			w_comb = w_comb[ok]
			reduced = private$reduce_design_matrix_once(X_comb, j_treat, cache_key = "kk_ols_combined_reduced_design_weighted")
			X_comb = reduced$X
			j_treat = reduced$j_treat
			fit = tryCatch(stats::lm.wfit(x = X_comb, y = y_comb, w = w_comb), error = function(e) NULL)
			coef_vec = if (is.null(fit)) NULL else as.numeric(fit$coefficients)
			if (is.null(coef_vec) || length(coef_vec) < j_treat || !is.finite(coef_vec[j_treat])) {
				return(if (estimate_only) NA_real_ else list(beta = NA_real_, se = NA_real_))
			}
			beta = coef_vec[j_treat]
			if (estimate_only) return(beta)
			df = nrow(X_comb) - ncol(X_comb)
			se = NA_real_
			if (df > 0) {
				sw = sqrt(w_comb)
				post_fit = tryCatch(
					ols_hc2_post_fit_cpp(X_comb * sw, as.numeric(y_comb * sw), coef_vec, j_treat),
					error = function(e) NULL
				)
				if (!is.null(post_fit)) {
					se_val = as.numeric(post_fit$std_err[j_treat])
					if (is.finite(se_val) && se_val > 0) se = se_val
				}
			}
			list(beta = beta, se = se)
		}
	)
)
