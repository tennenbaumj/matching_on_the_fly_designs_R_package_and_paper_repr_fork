#' OLS Inference for Continuous Responses
#'
#' Fits an ordinary least squares regression for continuous responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#' Note that warm starts are disabled for this class as OLS is a closed-form
#' estimator and does not benefit from initialization.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinOLS$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinOLS = R6::R6Class("InferenceContinOLS",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = list(
		#' @description Initialize an OLS inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param max_resample_attempts Maximum number of times a single bootstrap replicate
		#'   may be redrawn when the drawn sample fails validity screening. Default \code{50L}.
		#' @param harden  		Whether to apply robustness measures.
		#' @param smart_cold_start_default Flag for consistent API.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, max_resample_attempts = 50L, harden = TRUE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertCount(max_resample_attempts, positive = TRUE)
			}
			super$initialize(des_obj = des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$fit_warm_start_enabled = FALSE
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			private$max_resample_attempts = max_resample_attempts
		},
		#' @description Computes the OLS estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			cv = private$cached_values
			if (!is.null(cv$beta_hat_T)) {
				if (estimate_only || !is.null(cv$s_beta_hat_T)) return(cv$beta_hat_T)
			}
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_full = private$create_design_matrix()
			keep = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
			if (!any(keep)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			w_eff = as.numeric(row_weights[keep])
			X_sub = X_full[keep, , drop = FALSE]
			fit = tryCatch(
				stats::lm.wfit(x = X_sub, y = as.numeric(private$y[keep]), w = w_eff),
				error = function(e) NULL
			)
			coef_hat = if (!is.null(fit)) as.numeric(stats::coef(fit)) else numeric(0)
			if (length(coef_hat) < 2L || !is.finite(coef_hat[2L])) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = coef_hat[2L]
			if (!estimate_only) {
				df = sum(keep) - length(coef_hat)
				se = if (df > 0L) {
					sigma2_w = sum(w_eff * fit$residuals^2) / df
					var_j = tryCatch(solve(crossprod(X_sub * sqrt(w_eff)))[2L, 2L], error = function(e) NA_real_)
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
		#' @description Computes an approximate confidence interval for the treatment effect.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			cv = private$cached_values
			if (is.null(cv$s_beta_hat_T)) {
				private$shared(estimate_only = FALSE)
				cv = private$cached_values
			}
			est = cv$beta_hat_T
			se = cv$s_beta_hat_T
			df = cv$df
			if (is.null(df)) df = NA_real_
			if (length(est) != 1L || !is.finite(est) || !is.finite(se) || se <= 0) return(c(NA_real_, NA_real_))
			mult = if (!is.finite(df)) stats::qnorm(1 - alpha / 2) else stats::qt(1 - alpha / 2, df = df)
			ci = c(est - mult * se, est + mult * se)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		#' @description Computes an approximate two-sided p-value for the treatment effect.
		#' @param delta The null difference to test against. Default is zero.
		compute_asymp_two_sided_pval = function(delta = 0){
			cv = private$cached_values
			if (is.null(cv$s_beta_hat_T)) {
				private$shared(estimate_only = FALSE)
				cv = private$cached_values
			}
			est = cv$beta_hat_T
			se = cv$s_beta_hat_T
			df = cv$df
			if (is.null(df)) df = NA_real_
			if (length(est) != 1L || !is.finite(est) || !is.finite(se) || se <= 0) return(NA_real_)
			val = (est - delta) / se
			if (is.finite(df)) 2 * stats::pt(-abs(val), df = df) else 2 * stats::pnorm(-abs(val))
		}
	),
	private = list(
		max_resample_attempts = NULL,
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# the OLS kernel only implements the additive sharp-null shift
			if (delta != 0 && !identical(transform_responses, "none")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_full = private$create_design_matrix()
			Xc = if (ncol(X_full) > 2L) {
				as.matrix(X_full[, -(1:2), drop = FALSE])
			} else {
				matrix(numeric(0), nrow = as.integer(private$n), ncol = 0L)
			}
			compute_rand_bootstrap_ols_parallel_cpp(
				as.numeric(y0_full), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
		# Affine decomposition of the BRT null draws for the closed-form CI. The statistic is
		# the treatment coefficient of OLS on [1, w_fresh, Xc]; shifting responses by
		# delta * (w_fresh - w_obs) moves that coefficient by delta * (1 - g_b) where g_b is
		# the w_fresh-coefficient of regressing w_obs on the same design matrix. Both A_b and
		# g_b come from one QR per draw with a two-column response.
		compute_rand_bootstrap_ci_affine_coefs = function(rand_bootstrap_draws){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			n = as.integer(private$n)
			B = length(rand_bootstrap_draws)
			if (B == 0L) return(NULL)
			y_raw = as.numeric(private$y)
			w_obs = as.numeric(private$w)
			# use the same (hardened) covariate columns as the full-data design matrix
			X_full = private$create_design_matrix()
			Xc = if (ncol(X_full) > 2L) as.matrix(X_full[, -(1:2), drop = FALSE]) else NULL
			A = rep(NA_real_, B)
			cc = rep(NA_real_, B)
			for (b in seq_len(B)) {
				draw = rand_bootstrap_draws[[b]]
				if (is.null(draw$w_b) || length(draw$i_b) != n || length(draw$w_b) != n) return(NULL)
				w_f = as.numeric(draw$w_b)
				n_T = sum(w_f == 1)
				if (n_T == 0L || n_T == n) next
				M_b = if (is.null(Xc)) {
					cbind(1, w_f)
				} else {
					cbind(1, w_f, Xc[draw$i_b, , drop = FALSE])
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
		get_standard_error = function(){
			if (is.null(private$cached_values$s_beta_hat_T)) private$shared()
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			if (is.null(private$cached_values$df)) private$shared()
			private$cached_values$df
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			X_full = private$create_design_matrix()
			
			if (!private$harden) {
				if (estimate_only) {
					res = fast_ols_cpp(X_full, private$y)
					private$cached_values$beta_hat_T = as.numeric(res$b[2])
				} else {
					res = fast_ols_with_var_cpp(X_full, private$y, j = 2L)
					private$cached_values$beta_hat_T = as.numeric(res$b[2])
					private$cached_values$s_beta_hat_T = if (is.finite(res$ssq_b_2)) sqrt(res$ssq_b_2) else NA_real_
					private$cached_values$df = nrow(X_full) - ncol(X_full)
					private$cached_values$likelihood_test_context = list(X = X_full, j_treat = 2L, full_fit = res)
					private$cached_mod = res
				}
				return(invisible(NULL))
			}

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit){
					if (estimate_only) {
						res = fast_ols_cpp(X_fit, private$y)
						list(b = res$b, XtX = res$XtX, ssq_b_2 = NA_real_)
					} else {
						fast_ols_with_var_cpp(X_fit, private$y, j = 2L)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && length(mod$b) >= 2L && is.finite(mod$b[2])
				}
			)
			if (is.null(attempt$fit) || !is.finite(attempt$fit$b[2])){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[2])
			if (estimate_only) return(invisible(NULL))

			private$cached_values$s_beta_hat_T = if (is.finite(attempt$fit$ssq_b_2)) sqrt(attempt$fit$ssq_b_2) else NA_real_
			private$cached_values$df = nrow(attempt$X) - ncol(attempt$X)
			
			private$cached_values$likelihood_test_context = list(
				X = attempt$X,
				j_treat = 2L,
				full_fit = attempt$fit
			)
			private$cached_mod = attempt$fit
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		supports_likelihood_tests = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$b)
			# Residual variance estimate from full fit
			y_orig = as.numeric(private$y)
			X_orig = spec$X
			# Sigma_hat^2 = RSS / (n - p)
			rss = sum((y_orig - as.numeric(X_orig %*% as.numeric(spec$full_fit$b)))^2)
			sigma = sqrt(rss / (nrow(X_orig) - ncol(X_orig)))
			
			mu = as.numeric(spec$X %*% b_null)
			y_sim = mu + rnorm(length(mu), mean = 0, sd = sigma)
			
			X_fit = spec$X
			j = spec$j
			
			full_fit_boot = fast_ols_cpp(X_fit, y_sim)
			
			sig2 = spec$full_fit$sigma2_hat
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					y_null = y_sim - as.numeric(X_fit[, j] * d)
					res = fast_ols_cpp(X_fit[, -j, drop = FALSE], y_null)
					b_full = numeric(ncol(X_fit))
					b_full[j] = d
					b_full[-j] = as.numeric(res$b)
					list(b = b_full, rss = sum((y_null - as.numeric(X_fit[, -j, drop = FALSE] %*% res$b))^2))
				},
				neg_loglik = function(fit){
					rss = if (!is.null(fit$rss)) fit$rss else sum((y_sim - as.numeric(X_fit %*% as.numeric(fit$b)))^2)
					0.5 * rss / sig2
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
			sig2 = private$cached_mod$sigma2_hat
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					y_null = y - as.numeric(X_fit[, j_treat] * delta)
					res = fast_ols_cpp(X_fit[, -j_treat, drop = FALSE], y_null)
					b_full = numeric(ncol(X_fit))
					b_full[j_treat] = delta
					b_full[-j_treat] = as.numeric(res$b)
					list(b = b_full, rss = sum((y_null - as.numeric(X_fit[, -j_treat, drop = FALSE] %*% res$b))^2))
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
