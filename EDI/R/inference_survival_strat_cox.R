#' Stratified Cox PH Inference for Survival Responses
#'
#' Fits an auto-stratified Cox PH regression. Stratification variables are chosen
#' automatically from the recorded low-cardinality covariates. If no suitable
#' stratification covariates are found, the fit falls back to the corresponding
#' standard Cox PH model.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalStratCoxPHRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceSurvivalStratCoxPHRegr = R6::R6Class("InferenceSurvivalStratCoxPHRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   covariates from the design object are included. Use \code{~ 1} for univariate.
		#' @param use_rcpp Logical. If \code{TRUE} (default), enable internal Rcpp score/information helpers for likelihood inference.
		#'   Cox optimization uses \pkg{survival::coxph.fit}.
		#' @param optimization_alg Optimization algorithm: \code{"newton_raphson"} (default) or \code{"lbfgs"}.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = "lbfgs", verbose = FALSE, smart_cold_start_default = NULL) {
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFlag(use_rcpp)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$use_rcpp = use_rcpp
		},
		#' @description Recomputes the stratified Cox PH treatment estimate under
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
			X_cov = private$reduce_covariates_preserving_treatment(X_cov)
			X_fit = if (!is.null(X_cov) && ncol(X_cov) > 0) cbind(treatment = private$w, X_cov) else matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			strata_info = private$compute_strata_info(X_cov)
			use_strata = isTRUE(strata_info$num_strata > 1L)
			fit = weighted_cox_bootstrap_surrogate_fit(
				private$y, private$dead, X_fit, row_weights,
				strata = if (use_strata) strata_info$strata_id else NULL,
				warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(X_fit)) %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit))
			)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Compute confidence interval rand
		#' @param alpha The significance level (default 0.05).
		#' @param r Number of vectors to draw.
		#' @param pval_epsilon The bisection convergence tolerance.
		#' @param show_progress Whether to show a progress bar.
		#' @param ci_search_control Unused.
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			stop("Randomization confidence intervals are not supported for stratified Cox PH models because the estimator units (Log-Hazard Ratio) are inconsistent with the randomization test's required transformed scale (Log-Time Ratio / AFT effect).")
		}
	),
		private = list(
		use_rcpp = TRUE,
		strat_cox_X_linear_cache = NULL,
		strat_cox_strata_info_cache = NULL,
		strat_cox_informative_rows_cache = NULL,
		strat_cox_w_cache = NULL,
		strat_cox_data_cache = NULL,
		strat_cox_X_fit_cache = NULL,
		strat_cox_y_cache = NULL,
		strat_cox_dead_cache = NULL,
		strat_cox_strata_sub_cache = NULL,
		cached_mod = NULL,
		coxph_control = NULL,
		get_complexity_tier = function() "light",
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T) && is.finite(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			mod = private$generate_mod(estimate_only = estimate_only)
			private$cached_values$beta_hat_T = mod$beta_hat_T %||% as.numeric(mod$b[2])
			if (!is.finite(private$cached_values$beta_hat_T)) {
				private$cache_nonestimable_estimate("strat_cox_fit_unavailable")
				return(invisible(NULL))
			}
			if (abs(private$cached_values$beta_hat_T) > 0.5) {
				private$cache_nonestimable_estimate("strat_cox_extreme_estimate")
				return(invisible(NULL))
			}
			if (estimate_only) return(invisible(NULL))
			se = if (is.finite(mod$ssq_b_2 %||% NA_real_) && mod$ssq_b_2 > 0) sqrt(mod$ssq_b_2) else NA_real_
			private$cached_values$s_beta_hat_T = se
			private$cached_values$df = NA_real_
			if (!is.finite(se)) {
				private$cache_nonestimable_se("strat_cox_standard_error_unavailable")
			}
		},
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$coefficients %||% null_fit$b)
			if (!all(is.finite(b_null))) return(NULL)
			X_fit = spec$X
			y_obs = as.numeric(spec$y)
			dead_obs = as.numeric(spec$dead)
			strata = as.integer(spec$strata)
			stratified = isTRUE(spec$stratified)
			j = spec$j
			if (stratified && !is.null(strata)){
				sim = .cox_simulate_stratified(y_obs, dead_obs, X_fit, b_null, strata)
			} else {
				breslow = .breslow_hazard(y_obs, dead_obs, X_fit, b_null)
				if (length(breslow$times) == 0L) return(NULL)
				sim = .cox_simulate_from_breslow(breslow, y_obs, dead_obs, X_fit, b_null)
			}
			y_sim = sim$y_sim; dead_sim = sim$dead_sim
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			strata_arg = if (stratified && !is.null(strata)) strata else NULL
			full_res = tryCatch(
				fast_stratified_coxph_regression_cpp(X_fit, y_sim, dead_sim, strata = as.integer(strata_arg %||% rep.int(1L, length(y_sim))), estimate_only = FALSE),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$coefficients[j])) {
				full_res = .fit_survival_coxph_kernel(X_fit, y_sim, dead_sim, strata = strata_arg)
			}
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$coefficients[j])) return(NULL)
			full_fit_boot = list(b = as.numeric(full_res$coefficients), neg_loglik = as.numeric(full_res$neg_ll))
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = .fit_survival_coxph_fixed_kernel(X_fit, y_sim, dead_sim, strata = strata_arg, fixed_idx = j, fixed_value = d)
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
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			dead = as.numeric(ctx$dead)
			strata = as.integer(ctx$strata)
			stratified = isTRUE(ctx$stratified)
			j_treat = as.integer(ctx$j_treat %||% 1L)
			list(
				X = X_fit,
				y = y,
				dead = dead,
				strata = strata,
				stratified = stratified,
				j = j_treat,
				full_fit = list(
					b = as.numeric(private$cached_mod$coefficients %||% private$cached_mod$b),
					neg_loglik = as.numeric(ctx$full_neg_loglik)
				),
				fit_null = function(delta, start = NULL){
					strata_arg = if (stratified) strata else NULL
					.fit_survival_coxph_fixed_kernel(X_fit, y, dead, strata = strata_arg, fixed_idx = j_treat, fixed_value = delta)
				},
				extract_start = function(fit){
					as.numeric(fit$coefficients %||% fit$b)
				},
				score = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						get_stratified_coxph_score_cpp(X_fit, y, dead, as.integer(strata), beta)
					} else {
						get_coxph_score_cpp(X_fit, y, dead, beta)
					}
				},
				observed_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						get_stratified_coxph_hessian_cpp(X_fit, y, dead, as.integer(strata), beta)
					} else {
						get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				fisher_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					fit$fisher_information %||% if (stratified) {
						get_stratified_coxph_hessian_cpp(X_fit, y, dead, as.integer(strata), beta)
					} else {
						get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					fit$information %||% fit$fisher_information %||% if (stratified) {
						get_stratified_coxph_hessian_cpp(X_fit, y, dead, as.integer(strata), beta)
					} else {
						get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_ll %||% fit$neg_loglik %||% fit$neg_log_lik)
				}
			)
		},
		compute_strata_info = function(X_full) {
			n = length(private$y)
			if (is.null(X_full) || ncol(X_full) == 0){
				return(list(strata_id = rep.int(1L, n), selected_cols = integer(0), num_strata = 1L))
			}
			info = tryCatch(
				compute_survival_strata_ids_cpp(as.matrix(X_full)),
				error = function(e) NULL
			)
			if (is.null(info)){
				return(list(strata_id = rep.int(1L, n), selected_cols = integer(0), num_strata = 1L))
			}
			list(
				strata_id = as.integer(info$strata_id),
				selected_cols = as.integer(info$selected_cols),
				num_strata = as.integer(info$num_strata)
			)
		},
		reduce_covariates_preserving_treatment = function(X_covars){
			if (is.null(X_covars) || ncol(X_covars) == 0){
				return(matrix(numeric(0), nrow = length(private$y), ncol = 0))
			}
			X_covars = as.matrix(X_covars)
			if (ncol(X_covars) == 0){
				return(matrix(numeric(0), nrow = nrow(X_covars), ncol = 0))
			}
			full_design = cbind(w = private$w, X_covars)
			reduced = drop_linearly_dependent_cols(full_design)
			X_keep = reduced$M
			if (ncol(X_keep) == 0){
				return(matrix(numeric(0), nrow = nrow(full_design), ncol = 0))
			}
			if (!("w" %in% colnames(X_keep))){
				return(matrix(numeric(0), nrow = nrow(full_design), ncol = 0))
			}
			X_keep[, colnames(X_keep) != "w", drop = FALSE]
		},
		get_informative_rows = function(strata_id){
			if (length(strata_id) != length(private$y)) return(integer(0))
			good = rep(FALSE, length(strata_id))
			for (s in unique(strata_id)){
				i_s = which(strata_id == s)
				if (length(i_s) < 2) next
				if (length(unique(private$w[i_s])) < 2) next
				if (!any(private$dead[i_s] == 1, na.rm = TRUE)) next
				good[i_s] = TRUE
			}
			which(good)
		},
		fit_cox_with_formula = function(dat, formula_str){
			tryCatch(
				suppressWarnings(survival::coxph(stats::as.formula(formula_str), data = dat)),
				error = function(e) NULL
			)
		},
		format_mod_output = function(mod){
			if (is.null(mod)){
				return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_, neg_log_lik = NA_real_))
			}
			coef_w = tryCatch(as.numeric(stats::coef(mod)["w"]), error = function(e) NA_real_)
			if (!is.finite(coef_w) || abs(coef_w) > 0.5) {
				return(list(beta_hat_T = NA_real_, b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_, neg_log_lik = tryCatch(as.numeric(-stats::logLik(mod)), error = function(e) NA_real_)))
			}
			ssq_w = tryCatch(as.numeric(stats::vcov(mod)["w", "w"]), error = function(e) NA_real_)
			list(
				beta_hat_T = coef_w,
				b = c(0, coef_w),
				ssq_b_2 = if (is.finite(ssq_w) && ssq_w > 0) ssq_w else NA_real_,
				neg_log_lik = tryCatch(as.numeric(-stats::logLik(mod)), error = function(e) NA_real_)
			)
		},
		# Build (y, dead, X_mat) for a given row set.  X_mat = cbind(w, X_linear).
		build_rcpp_inputs = function(rows, X_linear){
			y_r    = as.numeric(private$y[rows])
			dead_r = as.numeric(private$dead[rows])
			w_r    = private$w[rows]
			if (ncol(X_linear) > 0){
				X_mat = cbind(w = w_r, X_linear[rows, , drop = FALSE])
			} else {
				X_mat = matrix(w_r, ncol = 1L, dimnames = list(NULL, "w"))
			}
			list(
				y = y_r,
				dead = dead_r,
				X = as.matrix(X_mat),
				surv_y = survival::Surv(y_r, dead_r),
				rownames = as.character(seq_along(y_r))
			)
		},
		fit_coxph_estimate_only_fast = function(X, surv_y, strata = NULL, rownames = NULL){
			if (is.null(private$coxph_control)) {
				private$coxph_control = survival::coxph.control()
			}
			fit = tryCatch(
				survival::coxph.fit(
					x = X,
					y = surv_y,
					strata = if (is.null(strata)) NULL else as.integer(strata),
					offset = NULL,
					init = NULL,
					control = private$coxph_control,
					weights = NULL,
					method = "breslow",
					rownames = rownames %||% as.character(seq_len(nrow(X))),
					resid = FALSE
				),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NULL)
			b = as.numeric(fit$coefficients %||% numeric(0))
			if (length(b) != ncol(X) || !all(is.finite(b))) return(NULL)
			names(b) = colnames(X)
			list(b = b, coefficients = b, vcov = NULL, var = NULL,
				neg_ll = -as.numeric(utils::tail(fit$loglik, 1L)), 
				neg_loglik = -as.numeric(utils::tail(fit$loglik, 1L)), 
				neg_log_lik = -as.numeric(utils::tail(fit$loglik, 1L)),
				fisher_information = NULL, converged = TRUE)
		},
		fit_rcpp_stratified = function(rows, X_linear, strata_id, estimate_only = FALSE){
			inp  = private$build_rcpp_inputs(rows, X_linear)
			strata_sub = as.integer(strata_id[rows])
			fit = tryCatch(
				fast_stratified_coxph_regression_cpp(
					inp$X, inp$y, inp$dead, strata = strata_sub,
					estimate_only = estimate_only,
					warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(inp$X)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(inp$X)),
					smart_cold_start = private$smart_cold_start_default %||% FALSE
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				# Fallback to R if C++ fails
				fit_r = if (estimate_only) {
					private$fit_coxph_estimate_only_fast(inp$X, inp$surv_y, strata = strata_sub, rownames = inp$rownames)
				} else {
					.fit_survival_coxph_kernel(inp$X, inp$y, inp$dead, strata = strata_sub, estimate_only = FALSE)
				}
				fit = fit_r
			}
			if (is.null(fit)) return(NULL)
			if (estimate_only) return(list(fit = fit, stratified = TRUE))
			list(fit = fit, X = inp$X, y = inp$y, dead = inp$dead, strata = strata_sub, stratified = TRUE)
		},

		fit_rcpp_unstratified = function(rows, X_linear, estimate_only = FALSE){
			inp = private$build_rcpp_inputs(rows, X_linear)
			fit = tryCatch(
				fast_coxph_regression_cpp(
					inp$X, inp$y, inp$dead,
					estimate_only = estimate_only,
					warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(inp$X)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(inp$X)),
					smart_cold_start = private$smart_cold_start_default %||% FALSE
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				fit_r = if (estimate_only) {
					private$fit_coxph_estimate_only_fast(inp$X, inp$surv_y, rownames = inp$rownames)
				} else {
					.fit_survival_coxph_kernel(inp$X, inp$y, inp$dead, estimate_only = FALSE)
				}
				fit = fit_r
			}
			if (is.null(fit)) return(NULL)
			if (estimate_only) return(list(fit = fit, stratified = FALSE))
			list(fit = fit, X = inp$X, y = inp$y, dead = inp$dead, strata = NULL, stratified = FALSE)
		},

		format_rcpp_output = function(fit){
			if (is.null(fit) || !isTRUE(fit$converged)) return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_, neg_log_lik = NA_real_))
			beta_w  = as.numeric(fit$coefficients %||% fit$b)[1L]
			if (!is.finite(beta_w) || abs(beta_w) > 0.5) {
				return(list(beta_hat_T = NA_real_, b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_, neg_log_lik = as.numeric(fit$neg_ll %||% fit$neg_log_lik), fisher_information = fit$fisher_information))
			}
			ssq_w   = if (!is.null(fit$vcov) && nrow(fit$vcov) >= 1L && is.finite(fit$vcov[1, 1]) && fit$vcov[1, 1] > 0)
				fit$vcov[1, 1] else NA_real_
			list(beta_hat_T = beta_w, b = c(0, beta_w), ssq_b_2 = ssq_w, neg_log_lik = as.numeric(fit$neg_ll %||% fit$neg_log_lik), fisher_information = fit$fisher_information)
		},
		generate_mod = function(estimate_only = FALSE){
			if (is.null(private$strat_cox_X_linear_cache) || is.null(private$strat_cox_data_cache) || !identical(private$w, private$strat_cox_w_cache)) {
				X_full      = private$X
				private$strat_cox_strata_info_cache = private$compute_strata_info(X_full)
				private$strat_cox_X_linear_cache = matrix(numeric(0), nrow = length(private$y), ncol = 0)
				if (ncol(as.matrix(private$X)) > 0 && !is.null(X_full) && ncol(X_full) > 0){
					keep_cols = setdiff(seq_len(ncol(X_full)), private$strat_cox_strata_info_cache$selected_cols)
					if (length(keep_cols) > 0){
						private$strat_cox_X_linear_cache = private$reduce_covariates_preserving_treatment(X_full[, keep_cols, drop = FALSE])
					}
				}
				private$strat_cox_informative_rows_cache = integer(0)
				if (!is.null(private$strat_cox_strata_info_cache$strata_id) && isTRUE(private$strat_cox_strata_info_cache$num_strata > 1L)){
					private$strat_cox_informative_rows_cache = private$get_informative_rows(private$strat_cox_strata_info_cache$strata_id)
				}
				# Build prebuilt C++ data cache to avoid per-call strata grouping and sort overhead
				X_lin_c = private$strat_cox_X_linear_cache
				si_c    = private$strat_cox_strata_info_cache
				inf_r   = private$strat_cox_informative_rows_cache
				if (length(inf_r) >= 4L) {
					w_s = private$w[inf_r]
					X_f = as.matrix(if (ncol(X_lin_c) > 0) cbind(w = w_s, X_lin_c[inf_r, , drop = FALSE]) else matrix(w_s, ncol = 1L, dimnames = list(NULL, "w")))
					y_s = as.numeric(private$y[inf_r])
					d_s = as.numeric(private$dead[inf_r])
					st_s = as.integer(si_c$strata_id[inf_r])
					private$strat_cox_X_fit_cache     = X_f
					private$strat_cox_y_cache         = y_s
					private$strat_cox_dead_cache      = d_s
					private$strat_cox_strata_sub_cache = st_s
					private$strat_cox_data_cache = tryCatch(build_stratified_cox_data_cache_cpp(X_f, y_s, d_s, st_s), error = function(e) NULL)
				} else {
					X_f = as.matrix(if (ncol(X_lin_c) > 0) cbind(w = private$w, X_lin_c) else matrix(private$w, ncol = 1L, dimnames = list(NULL, "w")))
					y_a = as.numeric(private$y)
					d_a = as.numeric(private$dead)
					private$strat_cox_X_fit_cache     = X_f
					private$strat_cox_y_cache         = y_a
					private$strat_cox_dead_cache      = d_a
					private$strat_cox_strata_sub_cache = NULL
					private$strat_cox_data_cache = tryCatch(build_cox_data_cache_cpp(X_f, y_a, d_a), error = function(e) NULL)
				}
				private$strat_cox_w_cache = private$w
			}
			
			X_linear = private$strat_cox_X_linear_cache
			strata_info = private$strat_cox_strata_info_cache
			informative_rows = private$strat_cox_informative_rows_cache
			
			if (length(informative_rows) >= 4){
				if (private$use_rcpp){
					p_fit = ncol(private$strat_cox_X_fit_cache)
					fit = NULL
					if (!is.null(private$strat_cox_data_cache)) {
						fit = tryCatch(
							fast_coxph_regression_prebuilt_cpp(
								private$strat_cox_data_cache,
								estimate_only = estimate_only,
								warm_start_beta = private$get_fit_warm_start_for_length("params", p_fit),
								warm_start_fisher_info = private$get_fit_warm_start_fisher(p_fit),
								smart_cold_start = private$smart_cold_start_default %||% FALSE,
								optimization_alg = "newton_raphson"
							),
							error = function(e) NULL
						)
					}
					if (!is.null(fit) && isTRUE(fit$converged)){
						private$cached_mod = fit
						private$set_fit_warm_start(as.numeric(fit$coefficients %||% fit$b), "params", fisher = fit$fisher_information, force_pd = TRUE)
						if (!estimate_only) {
							private$cached_values$likelihood_test_context = list(
								X = private$strat_cox_X_fit_cache, y = private$strat_cox_y_cache,
								dead = private$strat_cox_dead_cache, strata = private$strat_cox_strata_sub_cache,
								stratified = TRUE, j_treat = 1L,
								full_neg_loglik = fit$neg_ll %||% fit$neg_log_lik
							)
						}
						return(private$format_rcpp_output(fit))
					}
					# Fallback: old per-call path (C++ prebuilt failed)
					res = private$fit_rcpp_stratified(informative_rows, X_linear, strata_info$strata_id, estimate_only = estimate_only)
					if (!is.null(res) && isTRUE(res$fit$converged)){
						private$cached_mod = res$fit
						private$set_fit_warm_start(as.numeric(res$fit$coefficients %||% res$fit$b), "params", fisher = res$fit$fisher_information, force_pd = TRUE)
						if (!estimate_only) {
							private$cached_values$likelihood_test_context = list(
								X = res$X, y = res$y, dead = res$dead, strata = res$strata,
								stratified = TRUE, j_treat = 1L,
								full_neg_loglik = res$fit$neg_ll %||% res$fit$neg_log_lik
							)
						}
						return(private$format_rcpp_output(res$fit))
					}
					res = private$fit_rcpp_unstratified(informative_rows, X_linear, estimate_only = estimate_only)
					if (!is.null(res) && isTRUE(res$fit$converged)){
						private$cached_mod = res$fit
						private$set_fit_warm_start(as.numeric(res$fit$coefficients %||% res$fit$b), "params", fisher = res$fit$fisher_information, force_pd = TRUE)
						if (!estimate_only) {
							private$cached_values$likelihood_test_context = list(
								X = res$X, y = res$y, dead = res$dead, strata = res$strata,
								stratified = FALSE, j_treat = 1L,
								full_neg_loglik = res$fit$neg_ll %||% res$fit$neg_log_lik
							)
						}
						return(private$format_rcpp_output(res$fit))
					}
				} else {
					colnames(X_linear) = paste0("x", seq_len(ncol(X_linear)))
					dat_full = data.frame(y = private$y, dead = private$dead, w = private$w)
					if (ncol(X_linear) > 0) dat_full = cbind(dat_full, as.data.frame(X_linear))
					base_terms   = c("w", setdiff(colnames(dat_full), c("y", "dead", "w")))
					base_formula = paste("survival::Surv(y, dead) ~", paste(base_terms, collapse = " + "))
					dat_strat    = dat_full[informative_rows, , drop = FALSE]
					dat_strat$strata_id = factor(strata_info$strata_id[informative_rows])
					mod = private$fit_cox_with_formula(dat_strat, paste(base_formula, "+ strata(strata_id)"))
					if (!is.null(mod)){
						private$cached_mod = mod
						private$cached_values$likelihood_test_context = list(
							X = as.matrix(cbind(w = private$w[informative_rows], X_linear[informative_rows, , drop = FALSE])),
							y = private$y[informative_rows],
							dead = private$dead[informative_rows],
							strata = strata_info$strata_id[informative_rows],
							stratified = TRUE,
							j_treat = 1L,
							full_neg_loglik = as.numeric(-stats::logLik(mod))
						)
						return(private$format_mod_output(mod))
					}
					if (ncol(X_linear) > 0){
						mod = private$fit_cox_with_formula(
							dat_strat[, c("y", "dead", "w", "strata_id"), drop = FALSE],
							"survival::Surv(y, dead) ~ w + strata(strata_id)"
						)
						if (!is.null(mod)){
							private$cached_mod = mod
							private$cached_values$likelihood_test_context = list(
								X = as.matrix(cbind(w = private$w[informative_rows])),
								y = private$y[informative_rows],
								dead = private$dead[informative_rows],
								strata = strata_info$strata_id[informative_rows],
								stratified = TRUE,
								j_treat = 1L,
								full_neg_loglik = as.numeric(-stats::logLik(mod))
							)
							return(private$format_mod_output(mod))
						}
					}
				}
			}
			all_rows = seq_len(length(private$y))
			if (private$use_rcpp){
				p_fit = ncol(private$strat_cox_X_fit_cache)
				fit = NULL
				if (!is.null(private$strat_cox_data_cache)) {
					fit = tryCatch(
						fast_coxph_regression_prebuilt_cpp(
							private$strat_cox_data_cache,
							estimate_only = estimate_only,
							warm_start_beta = private$get_fit_warm_start_for_length("params", p_fit),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(p_fit),
							smart_cold_start = private$smart_cold_start_default %||% FALSE,
							optimization_alg = "newton_raphson"
						),
						error = function(e) NULL
					)
				}
				if (!is.null(fit) && isTRUE(fit$converged)){
					private$cached_mod = fit
					private$set_fit_warm_start(as.numeric(fit$coefficients %||% fit$b), "params", fisher = fit$fisher_information, force_pd = TRUE)
					if (!estimate_only) {
						private$cached_values$likelihood_test_context = list(
							X = private$strat_cox_X_fit_cache, y = private$strat_cox_y_cache,
							dead = private$strat_cox_dead_cache, strata = NULL,
							stratified = FALSE, j_treat = 1L,
							full_neg_loglik = fit$neg_ll %||% fit$neg_log_lik
						)
					}
					return(private$format_rcpp_output(fit))
				}
				res = private$fit_rcpp_unstratified(all_rows, X_linear, estimate_only = estimate_only)
				if (!is.null(res) && isTRUE(res$fit$converged)){
					private$cached_mod = res$fit
					private$set_fit_warm_start(as.numeric(res$fit$coefficients %||% res$fit$b), "params", fisher = res$fit$fisher_information, force_pd = TRUE)
					if (!estimate_only) {
						private$cached_values$likelihood_test_context = list(
							X = res$X, y = res$y, dead = res$dead, strata = NULL,
							stratified = FALSE, j_treat = 1L,
							full_neg_loglik = res$fit$neg_ll %||% res$fit$neg_log_lik
						)
					}
					return(private$format_rcpp_output(res$fit))
				}
			} else {
				colnames(X_linear) = paste0("x", seq_len(ncol(X_linear)))
				dat_full = data.frame(y = private$y, dead = private$dead, w = private$w)
				if (ncol(X_linear) > 0) dat_full = cbind(dat_full, as.data.frame(X_linear))
				base_terms   = c("w", setdiff(colnames(dat_full), c("y", "dead", "w")))
				base_formula = paste("survival::Surv(y, dead) ~", paste(base_terms, collapse = " + "))
				mod = private$fit_cox_with_formula(dat_full, base_formula)
				if (is.null(mod) && ncol(X_linear) > 0){
					mod = private$fit_cox_with_formula(
						dat_full[, c("y", "dead", "w"), drop = FALSE],
						"survival::Surv(y, dead) ~ w"
					)
				}
				if (!is.null(mod)){
					private$cached_mod = mod
					private$cached_values$likelihood_test_context = list(
						X = as.matrix(cbind(w = private$w, X_linear)),
						y = private$y,
						dead = private$dead,
						strata = NULL,
						stratified = FALSE,
						j_treat = 1L,
						full_neg_loglik = as.numeric(-stats::logLik(mod))
					)
				}
				return(private$format_mod_output(mod))
			}
			# Highly collinear covariates can make the multivariable Cox fit fail or
			# produce a non-finite treatment variance.  Preserve a usable treatment
			# estimate by refitting the treatment effect alone before declaring the
			# inference non-estimable.
			X_treat = matrix(private$w, ncol = 1L, dimnames = list(NULL, "w"))
			if (length(informative_rows) >= 4L) {
				rows_t = informative_rows
				strata_t = as.integer(strata_info$strata_id[rows_t])
				X_t = X_treat[rows_t, , drop = FALSE]
				y_t = as.numeric(private$y[rows_t])
				d_t = as.numeric(private$dead[rows_t])
				fit_t = tryCatch(
					fast_stratified_coxph_regression_cpp(
						X_t, y_t, d_t, strata = strata_t,
						estimate_only = estimate_only,
						optimization_alg = "newton_raphson"
					),
					error = function(e) NULL
				)
				if (is.null(fit_t) || !isTRUE(fit_t$converged)) {
					fit_t = tryCatch(
						.fit_survival_coxph_kernel(X_t, y_t, d_t, strata = strata_t, estimate_only = FALSE),
						error = function(e) NULL
					)
				}
				if (!is.null(fit_t) && isTRUE(fit_t$converged)) {
					private$cached_mod = fit_t
					private$cached_values$likelihood_test_context = list(
						X = X_t, y = y_t, dead = d_t, strata = strata_t,
						stratified = TRUE, j_treat = 1L,
						full_neg_loglik = fit_t$neg_ll %||% fit_t$neg_log_lik
					)
					return(private$format_rcpp_output(fit_t))
				}
			} else {
				fit_t = tryCatch(
					fast_coxph_regression_cpp(
						X_treat, as.numeric(private$y), as.numeric(private$dead),
						estimate_only = estimate_only,
						optimization_alg = "newton_raphson"
					),
					error = function(e) NULL
				)
				if (is.null(fit_t) || !isTRUE(fit_t$converged)) {
					fit_t = tryCatch(
						.fit_survival_coxph_kernel(X_treat, as.numeric(private$y), as.numeric(private$dead), estimate_only = FALSE),
						error = function(e) NULL
					)
				}
				if (!is.null(fit_t) && isTRUE(fit_t$converged)) {
					private$cached_mod = fit_t
					private$cached_values$likelihood_test_context = list(
						X = X_treat, y = as.numeric(private$y), dead = as.numeric(private$dead),
						strata = NULL, stratified = FALSE, j_treat = 1L,
						full_neg_loglik = fit_t$neg_ll %||% fit_t$neg_log_lik
					)
					return(private$format_rcpp_output(fit_t))
				}
			}
			list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_, neg_log_lik = NA_real_)
		}
	)
)
