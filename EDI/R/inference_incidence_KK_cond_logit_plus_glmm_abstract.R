#' Abstract Conditional Logistic Plus GLMM Inference
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs and reservoir subjects.
#'
#' @keywords internal
InferenceAbstractKKCondLogitPlusGLMM = R6::R6Class("InferenceAbstractKKCondLogitPlusGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with an incidence or proportion response.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param max_abs_reasonable_coef Cap for reasonable coefficient estimates.
		#' @param max_abs_reasonable_se Cap for reasonable treatment standard errors.
		#' @param max_abs_log_sigma Cap for reasonable log random effect variance.
		#' @param verbose Logical. Whether to print progress messages.
		#' @param smart_cold_start_default Logical. Whether to use smart starting values for the optimizer.
		#' @param optimization_alg Character. Optimization algorithm (default "lbfgs").
		initialize = function(des_obj, model_formula = NULL, max_abs_reasonable_coef = 50, max_abs_reasonable_se = 10, max_abs_log_sigma = 8, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), c("incidence", "proportion"))
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE, default = "lbfgs")
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$max_abs_reasonable_coef = max_abs_reasonable_coef
			private$max_abs_reasonable_se = max_abs_reasonable_se
			private$max_abs_log_sigma = max_abs_log_sigma
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Returns the standard error of the treatment effect estimate.
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			se = private$cached_values$s_beta_hat_T
			if (is.null(se) || length(se) == 0L) {
				return(NA_real_)
			}
			as.numeric(se)[1L]
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Numeric vector. Row weights for bootstrap.
		#' @param estimate_only Logical. If TRUE, skip variance component calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$beta_hat_T = beta_hat_T
					private$cached_values$s_beta_hat_T = NA_real_
					return(private$cached_values$beta_hat_T)
				}
			}
			if (!check_package_installed("glmmTMB")) return(NA_real_)
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			res_idx = which(group_id == 0L)
			if (length(res_idx) > 0L) {
				group_id[res_idx] = max(group_id) + seq_along(res_idx)
			}
			X_fit = private$create_design_matrix()
			if (is.null(X_fit)) return(NA_real_)
			X_fit = as.data.frame(X_fit[, -1, drop = FALSE])
			if ("treatment" %in% colnames(X_fit)) colnames(X_fit)[colnames(X_fit) == "treatment"] = "w"
			ok = is.finite(row_weights) & row_weights > 0 & is.finite(as.numeric(private$y))
			if (!any(ok)) return(NA_real_)
			dat = data.frame(
				y = as.numeric(private$y[ok]),
				X_fit[ok, , drop = FALSE],
				group_id = factor(group_id[ok]),
				.bootstrap_weight__ = as.numeric(row_weights[ok])
			)
			formula_glmm = stats::as.formula(paste("y ~", paste(c(setdiff(colnames(dat), c("y", "group_id", ".bootstrap_weight__")), "(1 | group_id)"), collapse = " + ")))
			mod = tryCatch(
				suppressWarnings(glmmTMB::glmmTMB(formula_glmm, family = stats::binomial(link = "logit"), data = dat, weights = .bootstrap_weight__, se = FALSE)),
				error = function(e) NULL
			)
			beta = tryCatch(glmmTMB::fixef(mod)$cond[["w"]], error = function(e) NA_real_)
			max_abs_boot_coef = min(private$max_abs_reasonable_coef, private$bootstrap_extreme_estimate_threshold)
			if (!is.finite(beta) || abs(beta) > max_abs_boot_coef) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = if (is.finite(beta)) as.numeric(beta) else NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
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
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		max_abs_reasonable_coef = 50,
		max_abs_reasonable_se = 10,
		max_abs_log_sigma = 8,
		bootstrap_extreme_estimate_threshold = 8,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			private$clear_nonestimable_state()
			private$cached_mod = NULL
			private$cached_values$likelihood_test_context = NULL
			
			d = private$prepare_clogit_plus_glmm_data()
			if (!d$has_discordant && !d$has_concordant) {
				private$cache_nonestimable_estimate("no_data_for_clogit_plus_glmm")
				return(invisible(NULL))
			}
			
			n_params = ncol(d$X_conc) + 1L # betas + log_sigma
			fit = tryCatch(
				fast_clogit_plus_glmm_cpp(
					X_disc = d$X_disc, y_disc = d$y_disc,
					X_conc = d$X_conc, y_conc = d$y_conc,
					group_conc = d$group_conc,
					has_discordant = d$has_discordant,
					has_concordant = d$has_concordant,
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					estimate_only = estimate_only,
					max_abs_log_sigma = private$max_abs_log_sigma,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				private$cache_nonestimable_estimate("joint_likelihood_failed_to_converge")
				return(invisible(NULL))
			}
			
			# Treatment effect index:
			# If has_concordant is true, shared beta starts at par[0] (intercept), so treatment is par[1].
			# If has_concordant is false, shared beta starts at par[0], so treatment is par[0].
			j_T = if (d$has_concordant) 2L else 1L
			
			params = as.numeric(fit$params)
			coef_params = if (length(params) > 1L) params[-length(params)] else params
			if (any(!is.finite(coef_params)) || any(abs(coef_params) > private$max_abs_reasonable_coef)) {
				private$cache_nonestimable_estimate("joint_likelihood_nonestimable")
				return(invisible(NULL))
			}
			log_sigma = if (length(params) > 1L) params[length(params)] else NA_real_
			if (!is.finite(log_sigma) || abs(log_sigma) > private$max_abs_log_sigma) {
				private$cache_nonestimable_estimate("joint_likelihood_nonestimable")
				return(invisible(NULL))
			}
			beta_hat_T = as.numeric(params[j_T])
			
			private$cached_mod = fit
			private$set_fit_warm_start(as.numeric(fit$params), "params", fisher = fit$fisher_information)
			private$cached_values$likelihood_test_context = list(
				d = d,
				j_T = j_T,
				start = as.numeric(fit$params)
			)
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_j 
			se = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
			if (!is.finite(se) || se <= 0 || se > private$max_abs_reasonable_se) {
				private$cache_nonestimable_se("joint_likelihood_standard_error_unavailable")
				return(invisible(NULL))
			}
			private$cached_values$s_beta_hat_T = se
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			d = ctx$d
			j_treat = as.integer(ctx$j_T)
			list(
				X = d$X_conc, 
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					fast_clogit_plus_glmm_cpp(
						X_disc = d$X_disc, y_disc = d$y_disc,
						X_conc = d$X_conc, y_conc = d$y_conc,
						group_conc = d$group_conc,
						has_discordant = d$has_discordant,
						has_concordant = d$has_concordant,
						warm_start_params = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start,
						warm_start_fisher_info = private$get_fit_warm_start_fisher(length(ctx$start)),
						estimate_only = FALSE,
						max_abs_log_sigma = private$max_abs_log_sigma,
						fixed_idx = j_treat, fixed_values = delta,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){ as.numeric(fit$params) },
				score = function(fit){
					as.numeric(get_clogit_plus_glmm_score_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				observed_information = function(fit){
					as.matrix(get_clogit_plus_glmm_hessian_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				fisher_information = function(fit){
					as.matrix(get_clogit_plus_glmm_hessian_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				information = function(fit){
					as.matrix(get_clogit_plus_glmm_hessian_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		},
		prepare_clogit_plus_glmm_data = function(){
			private$compute_basic_match_data()
			KKstats = private$cached_values$KKstats
			
			# Discordant pairs
			yTs = as.numeric(KKstats$yTs_matched)
			yCs = as.numeric(KKstats$yCs_matched)
			y_m = yTs - yCs # 1 if (1,0), -1 if (0,1)
			
			i_m_disc = which(abs(y_m) == 1)
			X_disc_cov = KKstats$X_matched_diffs_full[i_m_disc, , drop = FALSE]
			# Prepend 1 for the treatment effect in clogit part. 
			# Use a guard to avoid recycling warnings if i_m_disc is empty.
			X_disc = if (length(i_m_disc) > 0L) {
				cbind(treatment = 1, X_disc_cov)
			} else {
				matrix(0, 0, ncol(X_disc_cov) + 1L)
			}
			# y for clogit should be 1 if (1,0) and 0 if (0,1)
			y_disc = (y_m[i_m_disc] + 1) / 2
			
			# Concordant pairs and reservoir for GLMM part
			m_vec = private$m
			# For proportions, differences might not be exactly 0, 1, or -1. 
			# We treat all non-discordant pairs as "concordant" (meaning they go to the GLMM).
			i_m_conc = setdiff(seq_along(y_m), i_m_disc)
			
			i_conc_pair = which(m_vec %in% i_m_conc)
			i_res = which(is.na(m_vec) | m_vec == 0L)
			
			i_glmm = if (private$combine_reservoir_into_glmm()) c(i_conc_pair, i_res) else i_conc_pair
			
			X_glmm_cov = private$get_X()[i_glmm, , drop = FALSE]
			y_glmm = private$y[i_glmm]
			w_glmm = private$w[i_glmm]
			group_glmm = m_vec[i_glmm]
			i_glmm_res = which(is.na(group_glmm) | group_glmm == 0L)
			if (length(i_glmm_res) > 0L) {
				group_glmm[i_glmm_res] = max(c(0L, m_vec), na.rm = TRUE) + seq_along(i_glmm_res)
			}
			
			X_glmm = if (length(i_glmm) > 0L) {
				cbind(`(Intercept)` = 1, treatment = w_glmm, X_glmm_cov)
			} else {
				matrix(0, 0, ncol(X_glmm_cov) + 2L)
			}
			
			list(
				X_disc = as.matrix(X_disc),
				y_disc = as.numeric(y_disc),
				X_conc = as.matrix(X_glmm),
				y_conc = as.numeric(y_glmm),
				group_conc = as.integer(group_glmm),
				has_discordant = length(i_m_disc) > 0L,
				has_concordant = length(i_glmm) > 0L
			)
		},
		combine_reservoir_into_glmm = function() stop("must implement combine_reservoir_into_glmm()"),
		log_sum_exp = function(x){
			m = max(x)
			if (!is.finite(m)) return(m)
			m + log(sum(exp(x - m)))
		},
		log1pexp = function(x){
			ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
		}
	)))
)
