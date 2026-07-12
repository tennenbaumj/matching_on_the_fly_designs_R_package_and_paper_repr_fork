# Internal survival::coxph.fit adapters used by Cox inference classes.
.fit_survival_coxph_kernel = function(X, y, dead, strata = NULL, offset = NULL, estimate_only = FALSE){
	X = as.matrix(X)
	y = as.numeric(y)
	dead = as.numeric(dead)
	if (is.null(colnames(X))) colnames(X) = paste0("x", seq_len(ncol(X)))
	strata_arg = if (is.null(strata)) NULL else as.integer(strata)
	offset_arg = if (is.null(offset)) NULL else as.numeric(offset)
	fit = tryCatch(
		suppressWarnings(survival::coxph.fit(
			x = X,
			y = survival::Surv(y, dead),
			strata = strata_arg,
			offset = offset_arg,
			init = NULL,
			control = survival::coxph.control(),
			weights = NULL,
			method = "breslow",
			rownames = as.character(seq_along(y)),
			resid = FALSE
		)),
		error = function(e) NULL
	)
	if (is.null(fit)) return(NULL)
	b = as.numeric(fit$coefficients %||% numeric(0))
	if (ncol(X) > 0L) {
		if (length(b) != ncol(X) || !all(is.finite(b))) return(NULL)
		names(b) = colnames(X)
	}
	if (estimate_only) {
		return(list(b = b, coefficients = b, vcov = NULL, var = NULL,
			neg_ll = NA_real_, neg_loglik = NA_real_, neg_log_lik = NA_real_,
			fisher_information = NULL, converged = TRUE))
	}
	vcov = if (ncol(X) > 0L) {
		v = as.matrix(fit$var)
		dimnames(v) = list(colnames(X), colnames(X))
		v
	} else {
		matrix(numeric(0), 0L, 0L)
	}
	neg_ll = -as.numeric(utils::tail(fit$loglik, 1L))
	if (!is.finite(neg_ll)) return(NULL)
	fisher = if (ncol(vcov) > 0L) {
		tryCatch(solve(vcov), error = function(e) NULL)
	} else {
		matrix(numeric(0), 0L, 0L)
	}
	list(
		b = b,
		coefficients = b,
		vcov = vcov,
		var = vcov,
		neg_ll = neg_ll,
		neg_loglik = neg_ll,
		neg_log_lik = neg_ll,
		fisher_information = fisher,
		converged = TRUE
	)
}

.fit_survival_coxph_fixed_kernel = function(X, y, dead, strata = NULL, fixed_idx = 1L, fixed_value = 0){
	X = as.matrix(X)
	p = ncol(X)
	if (p < 1L || fixed_idx < 1L || fixed_idx > p) return(NULL)
	if (is.null(colnames(X))) colnames(X) = paste0("x", seq_len(p))
	fixed_idx = as.integer(fixed_idx)
	X_fixed = X[, fixed_idx, drop = TRUE]
	free_idx = setdiff(seq_len(p), fixed_idx)
	X_free = X[, free_idx, drop = FALSE]
	offset = as.numeric(fixed_value) * as.numeric(X_fixed)
	fit = .fit_survival_coxph_kernel(X_free, y, dead, strata = strata, offset = offset)
	if (is.null(fit)) return(NULL)
	b = numeric(p)
	b[fixed_idx] = as.numeric(fixed_value)
	if (length(free_idx) > 0L) b[free_idx] = as.numeric(fit$b)
	names(b) = colnames(X)
	fisher = if (!is.null(fit$fisher_information) && length(free_idx) > 0L) {
		full = matrix(0, p, p, dimnames = list(colnames(X), colnames(X)))
		full[free_idx, free_idx] = fit$fisher_information
		full
	} else {
		NULL
	}
	list(
		b = b,
		coefficients = b,
		vcov = NULL,
		var = NULL,
		neg_ll = fit$neg_ll,
		neg_loglik = fit$neg_loglik,
		neg_log_lik = fit$neg_log_lik,
		fisher_information = fisher,
		converged = TRUE
	)
}

.cox_neg_loglik_breslow_r = function(X, y, dead, beta, strata = NULL){
	X = as.matrix(X)
	y = as.numeric(y)
	dead = as.numeric(dead)
	beta = as.numeric(beta)
	if (length(beta) != ncol(X)) return(NA_real_)
	strata_id = if (is.null(strata)) rep.int(1L, length(y)) else as.integer(strata)
	eta = as.numeric(X %*% beta)
	if (!all(is.finite(eta))) return(NA_real_)
	ll = 0
	for (s in unique(strata_id)) {
		idx = which(strata_id == s)
		if (length(idx) == 0L) next
		event_times = sort(unique(y[idx][dead[idx] == 1L]))
		for (tt in event_times) {
			event_idx = idx[dead[idx] == 1L & y[idx] == tt]
			risk_idx = idx[y[idx] >= tt]
			d = length(event_idx)
			if (d == 0L || length(risk_idx) == 0L) next
			den = sum(exp(eta[risk_idx]))
			if (!is.finite(den) || den <= 0) return(NA_real_)
			ll = ll + sum(eta[event_idx]) - d * log(den)
		}
	}
	-as.numeric(ll)
}

.cox_score_breslow_fd_r = function(X, y, dead, beta, strata = NULL){
	beta = as.numeric(beta)
	p = length(beta)
	score = numeric(p)
	for (j in seq_len(p)) {
		h = 1e-5 * max(1, abs(beta[j]))
		b_hi = beta; b_lo = beta
		b_hi[j] = b_hi[j] + h
		b_lo[j] = b_lo[j] - h
		nll_hi = .cox_neg_loglik_breslow_r(X, y, dead, b_hi, strata = strata)
		nll_lo = .cox_neg_loglik_breslow_r(X, y, dead, b_lo, strata = strata)
		if (!is.finite(nll_hi) || !is.finite(nll_lo)) return(rep(NA_real_, p))
		score[j] = -(nll_hi - nll_lo) / (2 * h)
	}
	score
}

.cox_information_breslow_fd_r = function(X, y, dead, beta, strata = NULL){
	beta = as.numeric(beta)
	p = length(beta)
	info = matrix(NA_real_, p, p)
	nll0 = .cox_neg_loglik_breslow_r(X, y, dead, beta, strata = strata)
	if (!is.finite(nll0)) return(info)
	h = 1e-4 * pmax(1, abs(beta))
	for (i in seq_len(p)) {
		for (j in i:p) {
			if (i == j) {
				b_hi = beta; b_lo = beta
				b_hi[i] = b_hi[i] + h[i]
				b_lo[i] = b_lo[i] - h[i]
				f_hi = .cox_neg_loglik_breslow_r(X, y, dead, b_hi, strata = strata)
				f_lo = .cox_neg_loglik_breslow_r(X, y, dead, b_lo, strata = strata)
				val = (f_hi - 2 * nll0 + f_lo) / (h[i]^2)
			} else {
				b_pp = beta; b_pm = beta; b_mp = beta; b_mm = beta
				b_pp[i] = b_pp[i] + h[i]; b_pp[j] = b_pp[j] + h[j]
				b_pm[i] = b_pm[i] + h[i]; b_pm[j] = b_pm[j] - h[j]
				b_mp[i] = b_mp[i] - h[i]; b_mp[j] = b_mp[j] + h[j]
				b_mm[i] = b_mm[i] - h[i]; b_mm[j] = b_mm[j] - h[j]
				f_pp = .cox_neg_loglik_breslow_r(X, y, dead, b_pp, strata = strata)
				f_pm = .cox_neg_loglik_breslow_r(X, y, dead, b_pm, strata = strata)
				f_mp = .cox_neg_loglik_breslow_r(X, y, dead, b_mp, strata = strata)
				f_mm = .cox_neg_loglik_breslow_r(X, y, dead, b_mm, strata = strata)
				val = (f_pp - f_pm - f_mp + f_mm) / (4 * h[i] * h[j])
			}
			info[i, j] = val
			info[j, i] = val
		}
	}
	info
}

#' Cox Proportional Hazards Regression Inference for Survival Responses
#'
#' Fits a Cox proportional hazards regression for survival responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalCoxPHRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalCoxPHRegr = R6::R6Class("InferenceSurvivalCoxPHRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(

		#' @description Initialize a Cox PH inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), enable internal Rcpp score/information helpers for likelihood inference.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			
			
			private$use_rcpp = use_rcpp
		},
		#' @description Computes the Cox PH estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			super$compute_estimate(estimate_only = estimate_only)
		},
		#' @description Recomputes the Cox PH treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) return(beta_hat_T)
			}
			X_cov = private$get_X()
			X_fit = if (!is.null(X_cov) && ncol(X_cov) > 0) cbind(treatment = private$w, X_cov) else matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			fit = weighted_cox_bootstrap_surrogate_fit(
				private$y, private$dead, X_fit, row_weights,
				warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(X_fit)) %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit))
			)
			beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			if (!is.finite(beta_hat_T) || private$cox_coefficients_extreme(beta_hat_T)) {
				private$cache_nonestimable_estimate("coxph_weighted_extreme_coefficients")
				beta_hat_T = NA_real_
			}
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		use_rcpp = TRUE,
		cox_extreme_coef_threshold = 20,
		cox_X_fit_cache = NULL,
		cox_data_cache = NULL,
		cox_w_cache = NULL,
		cox_coefficients_extreme = function(coefs){
			coefs = as.numeric(coefs)
			any(!is.finite(coefs)) || any(abs(coefs) > private$cox_extreme_coef_threshold)
		},
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$b)
			if (!all(is.finite(b_null))) return(NULL)
			X_fit = spec$X
			y_obs = as.numeric(private$y)
			dead_obs = as.numeric(private$dead)
			breslow = .breslow_hazard(y_obs, dead_obs, X_fit, b_null)
			if (length(breslow$times) == 0L) return(NULL)
			sim = .cox_simulate_from_breslow(breslow, y_obs, dead_obs, X_fit, b_null)
			y_sim = sim$y_sim; dead_sim = sim$dead_sim
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			j = spec$j
			full_res = .fit_survival_coxph_kernel(X_fit, y_sim, dead_sim)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$coefficients[j])) return(NULL)
			full_fit_boot = list(b = as.numeric(full_res$coefficients), neg_loglik = as.numeric(full_res$neg_ll))
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = .fit_survival_coxph_fixed_kernel(X_fit, y_sim, dead_sim, fixed_idx = j, fixed_value = d)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), neg_loglik = as.numeric(res$neg_ll))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X = ctx$X
			y = as.numeric(private$y)
			dead = as.numeric(private$dead)
			full_b_cox = as.numeric(private$cached_mod$b[-1])  # strip Cox no-intercept prefix
			full_fit = list(b = full_b_cox, neg_loglik = ctx$full_neg_loglik)
			list(
				X = X, y = y, j = 1L,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					res = .fit_survival_coxph_fixed_kernel(X, y, dead, fixed_idx = 1L, fixed_value = delta)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), neg_loglik = as.numeric(res$neg_ll), fisher_information = res$fisher_information)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_coxph_score_cpp(X, y, dead, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_coxph_hessian_cpp(X, y, dead, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					fit$fisher_information %||% -get_coxph_hessian_cpp(X, y, dead, as.numeric(fit$b))
				},
				information = function(fit){
					fit$information %||% fit$fisher_information %||% -get_coxph_hessian_cpp(X, y, dead, as.numeric(fit$b))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			if (is.null(private$cox_X_fit_cache) || is.null(private$cox_data_cache) || !identical(private$w, private$cox_w_cache)) {
				X_cov = private$get_X()
				private$cox_X_fit_cache = if (!is.null(X_cov) && ncol(X_cov) > 0){
					cbind(treatment = private$w, X_cov)
				} else {
					matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
				}
				if (private$harden && ncol(private$cox_X_fit_cache) > 1L) {
					orig_names = colnames(private$cox_X_fit_cache)
					reduced = qr_reduce_preserve_cols_cpp(as.matrix(private$cox_X_fit_cache), 1L)
					private$cox_X_fit_cache = as.matrix(reduced$X_reduced)
					colnames(private$cox_X_fit_cache) = orig_names[as.integer(reduced$keep)]
				}
				private$cox_data_cache = build_cox_data_cache_cpp(private$cox_X_fit_cache, private$y, private$dead)
				private$cox_w_cache = private$w
			}
			X_fit = private$cox_X_fit_cache

			if (private$use_rcpp) {
				fit = tryCatch(
					fast_coxph_regression_prebuilt_cpp(
						private$cox_data_cache,
						estimate_only = estimate_only,
						warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
						smart_cold_start = private$smart_cold_start_default %||% TRUE,
						optimization_alg = "newton_raphson"
					),
					error = function(e) NULL
				)
				
				if (is.null(fit) || !isTRUE(fit$converged)) {
					# Fallback to R if C++ fails
					fit = .fit_survival_coxph_kernel(X_fit, private$y, private$dead, estimate_only = estimate_only)
				}
				
				if (is.null(fit)) {
					private$cached_values$likelihood_test_context = NULL
					return(list(b = rep(NA_real_, ncol(X_fit) + 1L), vcov = matrix(NA_real_, ncol(X_fit) + 1L, ncol(X_fit) + 1L)))
				}
				
				private$cached_mod = fit
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					full_neg_loglik = fit$neg_ll %||% fit$neg_log_lik
				)
				
				coefs = as.numeric(fit$coefficients %||% fit$b)
				if (private$cox_coefficients_extreme(coefs)) {
					private$cache_nonestimable_estimate("coxph_extreme_coefficients")
					private$cached_values$likelihood_test_context = NULL
					return(list(
						beta_hat_T = NA_real_,
						ssq_b_2 = NA_real_,
						b = rep(NA_real_, ncol(X_fit) + 1L),
						params = rep(NA_real_, ncol(X_fit)),
						neg_log_lik = NA_real_,
						fisher_information = NULL,
						vcov = if (estimate_only) NULL else matrix(NA_real_, ncol(X_fit) + 1L, ncol(X_fit) + 1L)
					))
				}
				return(list(
					beta_hat_T = coefs[1L],
					ssq_b_2 = if (estimate_only) NA_real_ else fit$vcov[1, 1],
					b = c(0, coefs),
					params = coefs,
					neg_log_lik = as.numeric(fit$neg_ll %||% fit$neg_log_lik),
					fisher_information = fit$fisher_information,
					vcov = if (estimate_only) NULL else {
						v = matrix(0, ncol(X_fit) + 1, ncol(X_fit) + 1)
						v[2:(ncol(X_fit) + 1), 2:(ncol(X_fit) + 1)] = fit$vcov
						v
					}
				))
			}
			surv_obj = survival::Surv(private$y, private$dead)
			tryCatch({
				coxph_mod = suppressWarnings(survival::coxph(surv_obj ~ X_fit))
				if (estimate_only) {
					coefs = stats::coef(coxph_mod)
					if (private$cox_coefficients_extreme(coefs)) {
						private$cache_nonestimable_estimate("coxph_extreme_coefficients")
						return(list(beta_hat_T = NA_real_, b = rep(NA_real_, ncol(X_fit) + 1L), ssq_b_2 = NA_real_, vcov = NULL))
					}
					list(
						beta_hat_T = as.numeric(coefs[1]),
						b = c(0, coefs),
						ssq_b_2 = NA_real_,
						vcov = NULL
					)
				} else {
					coefs = stats::coef(coxph_mod)
					if (private$cox_coefficients_extreme(coefs)) {
						private$cache_nonestimable_estimate("coxph_extreme_coefficients")
						return(list(
							beta_hat_T = NA_real_,
							ssq_b_2 = NA_real_,
							b = rep(NA_real_, ncol(X_fit) + 1L),
							vcov = matrix(NA_real_, ncol(X_fit) + 1L, ncol(X_fit) + 1L),
							neg_log_lik = NA_real_
						))
					}
					vcov_mat = stats::vcov(coxph_mod)
					v = matrix(0, ncol(X_fit) + 1, ncol(X_fit) + 1)
					v[2:(ncol(X_fit) + 1), 2:(ncol(X_fit) + 1)] = vcov_mat
					list(
						beta_hat_T = as.numeric(coefs[1]),
						ssq_b_2 = as.numeric(vcov_mat[1, 1]),
						b = c(0, coefs),
						vcov = v,
						neg_log_lik = as.numeric(-stats::logLik(coxph_mod))
					)
				}
			}, error = function(e){
				list(
					b = rep(NA_real_, ncol(X_fit) + 1),
					vcov = matrix(NA_real_, ncol(X_fit) + 1, ncol(X_fit) + 1)
				)
			})
		}
	)
)
