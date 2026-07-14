#' Lin (2013) Covariate-Adjusted OLS Inference for Continuous Responses
#'
#' Fits the Lin (2013) covariate-adjusted linear estimator for continuous responses.
#' The working model includes an intercept, treatment indicator and, optionally,
#' centered covariates and treatment-by-centered-covariate interactions.
#' Inference uses HC2 heteroskedasticity-robust standard errors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinLin$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinLin = R6::R6Class("InferenceContinLin",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = list(
		#' @description Initialize a Lin (2013) inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param harden Flag for consistent API.
		#' @param smart_cold_start_default Flag for consistent API.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, harden = TRUE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj = des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes Lin's estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_full = private$build_lin_design_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			keep = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
			if (is.null(X_fit) || !is.finite(j_treat) || !any(keep)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			w_eff = as.numeric(row_weights[keep])
			X_sub = X_fit[keep, , drop = FALSE]
			fit = tryCatch(
				stats::lm.wfit(x = X_sub, y = as.numeric(private$y[keep]), w = w_eff),
				error = function(e) NULL
			)
			coef_hat = if (!is.null(fit)) as.numeric(stats::coef(fit)) else numeric(0)
			if (length(coef_hat) < j_treat || !is.finite(coef_hat[j_treat])) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = coef_hat[j_treat]
			if (!estimate_only) {
				df = sum(keep) - length(coef_hat)
				se = if (df > 0L) {
					sigma2_w = sum(w_eff * fit$residuals^2) / df
					var_j = tryCatch(solve(crossprod(X_sub * sqrt(w_eff)))[j_treat, j_treat], error = function(e) NA_real_)
					if (is.finite(var_j) && var_j > 0) sqrt(sigma2_w * var_j) else NA_real_
				} else NA_real_
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				private$cached_values$df = df
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
			}
			private$cached_values$beta_hat_T
		},
		#' @description Computes a 1 - \code{alpha} confidence interval using HC2 robust standard error.
		#' @param alpha The confidence level. The default is 0.05.
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
		#' @description Computes a two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect. Defaults to 0.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		get_standard_error = function(){
			if (is.null(private$cached_values$s_beta_hat_T)) private$shared()
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			if (is.null(private$cached_values$df)) private$shared()
			private$cached_values$df
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		# Affine decomposition of the BRT null draws for the closed-form CI. The Lin estimator
		# is OLS on [1, w_fresh, Xc_b, w_fresh*Xc_b]. Under the additive sharp-null shift,
		# t0_b(delta) = A_b + delta * c_b where A_b = coef_2(lm(y ~ M_b)) and
		# c_b = 1 - g_b, with g_b = coef_2(lm(w_obs ~ M_b)). The slope is 1 - g_b because
		# h_b = coef_2(lm(w_fresh ~ M_b)) = 1 for any full-rank M_b that contains w_fresh as
		# its second column (algebraic identity: (M'M)^{-1}M'M[:,2] = e_2). Both A_b and g_b
		# are recovered from one QR per draw with a two-column response.
		compute_rand_bootstrap_ci_affine_coefs = function(rand_bootstrap_draws){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			n = as.integer(private$n)
			B = length(rand_bootstrap_draws)
			if (B == 0L) return(NULL)
			y_raw = as.numeric(private$y)
			w_obs = as.numeric(private$w)
			Xc_info = private$get_centered_covariates()
			has_covariates = !is.null(Xc_info)
			A = rep(NA_real_, B)
			cc = rep(NA_real_, B)
			for (b in seq_len(B)) {
				draw = rand_bootstrap_draws[[b]]
				if (is.null(draw$w_b) || length(draw$i_b) != n || length(draw$w_b) != n) return(NULL)
				w_f = as.numeric(draw$w_b)
				n_T = sum(w_f == 1)
				if (n_T == 0L || n_T == n) next
				M_b = if (!has_covariates) {
					cbind(1, w_f)
				} else {
					Xc_b = Xc_info$Xc[draw$i_b, , drop = FALSE]
					cbind(1, w_f, Xc_b, w_f * Xc_b)
				}
				coefs = tryCatch(
					qr.coef(qr(M_b), cbind(y_raw[draw$i_b], w_obs[draw$i_b])),
					error = function(e) NULL
				)
				if (is.null(coefs) || anyNA(coefs[2L, ])) next
				A[b] = as.numeric(coefs[2L, 1L])
				cc[b] = 1 - as.numeric(coefs[2L, 2L])
			}
			list(A = A, c = cc)
		},
		build_lin_design_matrix = function(){
			Xc_info = private$get_centered_covariates()
			if (is.null(Xc_info)){
				X_lin = cbind(1, private$w)
				colnames(X_lin) = c("(Intercept)", "treatment")
				return(X_lin)
			}
			Xc = Xc_info$Xc
			X_int = Xc * private$w
			colnames(X_int) = paste0("treatment:", colnames(Xc))
			X_lin = cbind(1, private$w, Xc, X_int)
			colnames(X_lin)[1:2] = c("(Intercept)", "treatment")
			X_lin
		},
		get_centered_covariates = function(){
			des_priv = private$des_obj_priv_int
			# Design-level cache: NULL = not computed, list() = computed with p=0, list(Xc=...) = computed
			if (!is.null(des_priv$lin_centered_covariates)) {
				cached = des_priv$lin_centered_covariates
				return(if (length(cached) == 0L) NULL else cached)
			}
			X = as.matrix(private$X)
			p = ncol(X)
			if (p == 0L){
				des_priv$lin_centered_covariates = list()  # sentinel: computed, no covariates
				return(NULL)
			}
			if (is.null(colnames(X))){
				colnames(X) = paste0("x", seq_len(p))
			}
			Xc = scale(X, center = TRUE, scale = FALSE)
			Xc = as.matrix(Xc)
			colnames(Xc) = colnames(X)
			result = list(Xc = Xc)
			des_priv$lin_centered_covariates = result
			result
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && isTRUE(private$cached_values$lin_estimate_only_complete)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$lin_full_complete)) return(invisible(NULL))
			X_full = private$build_lin_design_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$cache_nonestimable_estimate("linear_model_design_unusable")
				if (estimate_only) private$cached_values$lin_estimate_only_complete = TRUE
				private$cached_values$lin_full_complete = TRUE
				return(invisible(NULL))
			}
			mod = stats::lm.fit(X_fit, private$y)
			coef_hat = as.numeric(mod$coefficients)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				private$cache_nonestimable_estimate("linear_model_coefficients_unavailable")
				if (estimate_only){
					private$cached_values$lin_estimate_only_complete = TRUE
				} else {
					private$cached_values$lin_full_complete = TRUE
				}
				return(invisible(NULL))
			}
			beta_hat = coef_hat[j_treat]
			private$cached_values$beta_hat_T = beta_hat
			if (estimate_only){
				private$cached_values$lin_estimate_only_complete = TRUE
				return(invisible(NULL))
			}
			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X_fit, as.numeric(private$y), coef_hat, j_treat),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				private$cache_nonestimable_estimate("linear_model_post_fit_unavailable")
				private$cached_values$lin_full_complete = TRUE
				return(invisible(NULL))
			}
			coef_names = colnames(X_fit)
			vcov_hc2 = post_fit$vcov
			ssq_hat = post_fit$ssq_hat
			if (!is.finite(beta_hat) || !is.finite(ssq_hat) || ssq_hat < 0){
				beta_hat = NA_real_
				ssq_hat = NA_real_
			}
			names(coef_hat) = coef_names
			colnames(vcov_hc2) = rownames(vcov_hc2) = coef_names
			std_err = post_fit$std_err
			names(std_err) = coef_names
			z_vals = post_fit$z_vals
			names(z_vals) = coef_names
			private$cached_values$beta_hat_T = beta_hat
			private$cached_values$s_beta_hat_T = if (is.finite(ssq_hat)) sqrt(ssq_hat) else NA_real_
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
			private$cached_values$full_coefficients = coef_hat
			private$cached_values$full_vcov = vcov_hc2
			private$cached_values$lin_estimate_only_complete = TRUE
			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)
			private$cached_values$summary_table = summary_table
			private$cached_values$lin_full_complete = TRUE

			private$cached_values$likelihood_test_context = list(
				X = X_fit,
				j_treat = j_treat,
				full_fit = list(b = coef_hat)
			)
			private$cached_mod = list(b = coef_hat)
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		supports_likelihood_tests = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			# Using OLS as generative model for Lin PB
			b_null = as.numeric(null_fit$b)
			y_orig = as.numeric(private$y)
			X_orig = spec$X
			rss = sum((y_orig - as.numeric(X_orig %*% as.numeric(spec$full_fit$b)))^2)
			sig2 = rss / (nrow(X_orig) - ncol(X_orig))
			sigma = sqrt(sig2)

			mu = as.numeric(spec$X %*% b_null)
			y_sim = mu + rnorm(length(mu), mean = 0, sd = sigma)

			X_fit = spec$X
			j = spec$j

			# full fit for simulated data
			mod_boot = stats::lm.fit(X_fit, y_sim)
			full_fit_boot = list(b = as.numeric(mod_boot$coefficients))

			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					y_null = y_sim - as.numeric(X_fit[, j] * d)
					res = stats::lm.fit(X_fit[, -j, drop = FALSE], y_null)
					b_full = numeric(ncol(X_fit))
					b_full[j] = d
					b_full[-j] = as.numeric(res$coefficients)
					list(b = b_full, rss = sum((y_null - as.numeric(X_fit[, -j, drop = FALSE] %*% res$coefficients))^2))
				},
				neg_loglik = function(fit){
					rss_fit = if (!is.null(fit$rss)) fit$rss else sum((y_sim - as.numeric(X_fit %*% as.numeric(fit$b)))^2)
					0.5 * rss_fit / sig2
				}
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			b_full = as.numeric(private$cached_mod$b)
			rss_full = sum((y - X_fit %*% b_full)^2)
			sig2 = rss_full / (length(y) - ncol(X_fit))
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					y_null = y - as.numeric(X_fit[, j_treat] * delta)
					res = stats::lm.fit(X_fit[, -j_treat, drop = FALSE], y_null)
					b_full = numeric(ncol(X_fit))
					b_full[j_treat] = delta
					b_full[-j_treat] = as.numeric(res$coefficients)
					list(b = b_full, rss = sum((y_null - as.numeric(X_fit[, -j_treat, drop = FALSE] %*% res$coefficients))^2))
				},
				score = function(fit){
					as.numeric(t(X_fit) %*% (y - X_fit %*% as.numeric(fit$b)) / sig2)
				},
				observed_information = function(fit) (t(X_fit) %*% X_fit) / sig2,
				fisher_information   = function(fit) (t(X_fit) %*% X_fit) / sig2,
				neg_loglik = function(fit){
					if (!is.null(fit$rss)) return(0.5 * fit$rss / sig2)
					0.5 * sum((y - as.numeric(X_fit %*% as.numeric(fit$b)))^2) / sig2
				}
			)
		}
	)
)
