#' Marginal Standardization / G-Computation for Binary Responses
#'
#' Internal base class for incidence-outcome g-computation estimators. A
#' logistic regression working model is fit, then potential-outcome risks under
#' all-treated and all-control assignments are standardized over the empirical
#' covariate distribution. Inference uses a sandwich-robust covariance for the
#' regression coefficients and the delta method for the marginal estimand.
#'
#' @details
#' The implementation is optimized for resampling-based inference (bootstrap and
#' randomization). It utilizes a fast C++ IRLS solver for the underlying logistic
#' regression. During resampling draws, it bypasses the expensive calculation of
#' the sandwich covariance matrix and delta-method standard errors, as only the
#' point estimate is required for the distribution.
#'
#' @keywords internal
#' @noRd
InferenceIncidGCompAbstract = R6::R6Class("InferenceIncidGCompAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
				
		#' @description Initialize the g-computation inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param prob_clip_eps Probability clamp applied to fitted values before computing
		#'   the sandwich covariance. Predicted probabilities (mu) are clipped to
		#'   \code{[prob_clip_eps, 1 - prob_clip_eps]} so that the IWLS weight
		#'   \eqn{mu(1-mu)} is bounded away from zero. Must be in \code{[0, 0.5)}.
		#'   Default \code{.Machine$double.eps} (essentially no clamping).
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, prob_clip_eps = .Machine$double.eps, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
				assertFormula(model_formula, null.ok = TRUE)
				assertNumber(prob_clip_eps, lower = 0, upper = 0.5)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$prob_clip_eps = prob_clip_eps
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		},
		#' @description Computes the g-computation treatment-effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$get_effect_estimate()
		},
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			estimand = private$get_estimand_type()
			se = if (identical(estimand, "RD")) {
				private$cached_values$se_rd
			} else {
				private$cached_values$rr * private$cached_values$se_log_rr
			}
			if (is.null(se) || length(se) == 0L) {
				return(NA_real_)
			}
			as.numeric(se)[1L]
		},
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			effects = private$weighted_gcomp_effects_from_row_weights(row_weights)
			if (is.null(effects) || !private$effects_are_usable(effects, estimate_only = TRUE)) {
				private$set_failed_fit_cache(inference_ready = FALSE)
				private$cached_values$beta_hat_T = NA_real_
				return(NA_real_)
			}
			private$cached_values = gcomp_cache_standardized_effects(
				private$cached_values,
				effects,
				inference_ready = FALSE
			)
			private$cached_values$beta_hat_T = if (identical(private$get_estimand_type(), "RD")) effects$rd else effects$rr
			private$cached_values$beta_hat_T
		},
		#' @description Computes a 1 - \code{alpha} confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared(estimate_only = FALSE)
			private$compute_effect_confidence_interval(alpha)
		},
		#' @description Computes a two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect. Defaults to 0 for RD and 1 for RR.
		compute_asymp_two_sided_pval = function(delta = NULL){
			private$shared(estimate_only = FALSE)
			private$compute_effect_pvalue(delta)
		},
		#' @description Computes a Wald two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect. Defaults to 0 for RD and 1 for RR.
		compute_wald_two_sided_pval = function(delta = NULL){
			private$shared(estimate_only = FALSE)
			private$compute_effect_pvalue(delta)
		},
		#' @description Computes a Wald confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_wald_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared(estimate_only = FALSE)
			private$compute_effect_confidence_interval(alpha)
		},
		#' @description Computes a bootstrap confidence interval for the treatment effect.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of bootstrap samples.
		#' @param type Bootstrap CI type. See \code{InferenceNonParamBootstrap$compute_bootstrap_confidence_interval}.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param min_number_usable_samples Minimum number of finite bootstrap samples required.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L){
			type_resolved = tolower(type %||% "percentile")
			if (identical(private$get_estimand_type(), "RR") && identical(type_resolved, "basic")) {
				return(private$compute_rr_bootstrap_basic_confidence_interval(alpha = alpha, B = B, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples))
			}
			super$compute_bootstrap_confidence_interval(alpha = alpha, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples)
		},
		#' @description Computes a bootstrap two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect. Defaults to 0 for RD and 1 for RR.
		#' @param B Number of bootstrap samples.
		#' @param type Bootstrap p-value type. See \code{InferenceNonParamBootstrap$compute_bootstrap_two_sided_pval}.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param min_number_usable_samples Minimum number of finite bootstrap samples required.
		compute_bootstrap_two_sided_pval = function(delta = NULL, B = 501, type = "symmetric", na.rm = FALSE, show_progress = TRUE, min_number_usable_samples = 5L){
			if (is.null(delta)){
				delta = private$default_null_value()
			}
			super$compute_bootstrap_two_sided_pval(delta = delta, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples)
		},
		#' @description Computes a Bayesian-bootstrap two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect. Defaults to 0 for RD and 1 for RR.
		#' @param B Number of Bayesian-bootstrap samples.
		#' @param type Bayesian-bootstrap p-value type. See \code{InferenceBayesianBootstrap$compute_bayesian_bootstrap_two_sided_pval}.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param min_number_usable_samples Minimum number of finite bootstrap samples required.
		compute_bayesian_bootstrap_two_sided_pval = function(delta = NULL, B = 501, type = NULL, na.rm = FALSE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
			if (is.null(delta)){
				delta = private$default_null_value()
			}
			super$compute_bayesian_bootstrap_two_sided_pval(delta = delta, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples, weighting_unit_type = weighting_unit_type)
		},
		#' @description Computes a Bayesian-bootstrap confidence interval for the treatment effect.
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of Bayesian-bootstrap samples.
		#' @param type Bayesian-bootstrap CI type. See \code{InferenceBayesianBootstrap$compute_bayesian_bootstrap_confidence_interval}.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param min_number_usable_samples Minimum number of finite bootstrap samples required.
		compute_bayesian_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
			type_resolved = tolower(type %||% "percentile")
			if (identical(private$get_estimand_type(), "RR") && type_resolved %in% c("basic", "wald")) {
				return(private$compute_rr_bayesian_bootstrap_log_confidence_interval(alpha = alpha, B = B, type = type_resolved, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples, weighting_unit_type = weighting_unit_type))
			}
			super$compute_bayesian_bootstrap_confidence_interval(alpha = alpha, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples, weighting_unit_type = weighting_unit_type)
		},
		#' @description Computes a jackknife-Wald two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect. Defaults to 0 for RD and 1 for RR.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_wald_two_sided_pval = function(delta = NULL, unit = "auto"){
			if (is.null(delta)){
				delta = private$default_null_value()
			}
			if (identical(private$get_estimand_type(), "RR")) {
				return(private$compute_rr_jackknife_wald_two_sided_pval(delta = delta, unit = unit))
			}
			super$compute_jackknife_wald_two_sided_pval(delta = delta, unit = unit)
		},
		#' @description Computes a jackknife-Wald confidence interval for the treatment effect.
		#' @param alpha Significance level. Default \code{0.05}.
		#' @param unit Deletion unit. Default \code{"auto"}.
		compute_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			if (identical(private$get_estimand_type(), "RR")) {
				return(private$compute_rr_jackknife_wald_confidence_interval(alpha = alpha, unit = unit))
			}
			super$compute_jackknife_wald_confidence_interval(alpha = alpha, unit = unit)
		}
	),
	private = list(
		best_X_colnames = NULL,
		gcomp_boot_beta = NULL,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			# Use the same design matrix structure as the original fit
			X_cols = private$best_X_colnames
			X_data = private$get_X()
			
			if (length(X_cols) == 0L){
				# Univariate case
				X = cbind("(Intercept)" = 1, treatment = private$w)
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind("(Intercept)" = 1, treatment = private$w, X_cov)
			}
			fit = tryCatch(
				fast_logistic_regression_cpp(
					X = X, 
					y = as.numeric(private$y),
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X))
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !private$coefficients_are_usable(as.numeric(fit$b))){
				return(NA_real_)
			}
			private$set_fit_warm_start(fit$b, "beta", fisher = fit$fisher_information)
			
			# Standardized effect
			coef_hat = as.numeric(fit$b)
			X1 = X
			X0 = X
			X1[, 2L] = 1
			X0[, 2L] = 0
			risk1 = mean(stats::plogis(as.numeric(X1 %*% coef_hat)))
			risk0 = mean(stats::plogis(as.numeric(X0 %*% coef_hat)))
			estimand = private$get_estimand_type()
			if (identical(estimand, "RD")) return(risk1 - risk0)
			if (risk1 > 0 && risk0 > 0) return(risk1 / risk0)
			NA_real_
		},
		prob_clip_eps = .Machine$double.eps,
		max_abs_reasonable_coef = 25,
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		},
		get_estimand_type = function() stop(class(self)[1], " must implement get_estimand_type()."),
		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},
		default_null_value = function(){
			if (identical(private$get_estimand_type(), "RR")) 1 else 0
		},
		compute_rr_bootstrap_basic_confidence_interval = function(alpha = 0.05, B = 501, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L){
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			est = as.numeric(self$compute_estimate(estimate_only = FALSE))[1L]
			if (!is.finite(est) || est <= 0) {
				return(private$missing_bootstrap_ci(alpha, "bootstrap_original_estimate_unavailable", stage = "estimate"))
			}
			boot_distr = as.numeric(self$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = show_progress))
			if (isTRUE(na.rm)) {
				boot_distr = boot_distr[is.finite(boot_distr) & boot_distr > 0]
			} else if (any(!is.finite(boot_distr) | boot_distr <= 0)) {
				return(ci)
			}
			if (length(boot_distr) < as.integer(min_number_usable_samples)) {
				return(private$missing_bootstrap_ci(alpha, "bootstrap_too_few_finite_estimates", stage = "estimate"))
			}
			q = stats::quantile(log(boot_distr), probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
			ci[] = exp(2 * log(est) - q)
			ci
		},
		compute_rr_bayesian_bootstrap_log_confidence_interval = function(alpha = 0.05, B = 501, type = "basic", na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			est = as.numeric(self$compute_estimate(estimate_only = FALSE))[1L]
			if (!is.finite(est) || est <= 0) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_original_estimate_unavailable", stage = "estimate"))
			}
			boot_distr = as.numeric(self$approximate_bayesian_bootstrap_distribution_beta_hat_T(B = B, show_progress = show_progress, weighting_unit_type = weighting_unit_type))
			if (isTRUE(na.rm)) {
				boot_distr = boot_distr[is.finite(boot_distr) & boot_distr > 0]
			} else if (any(!is.finite(boot_distr) | boot_distr <= 0)) {
				return(ci)
			}
			if (length(boot_distr) < as.integer(min_number_usable_samples)) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_too_few_finite_estimates", stage = "estimate"))
			}
			log_boot = log(boot_distr)
			if (identical(type, "wald")) {
				se_log = stats::sd(log_boot)
				if (!is.finite(se_log) || se_log <= 0) return(ci)
				z = stats::qnorm(1 - alpha / 2)
				ci[] = exp(log(est) + c(-1, 1) * z * se_log)
			} else {
				q = stats::quantile(log_boot, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
				ci[] = exp(2 * log(est) - q)
			}
			ci
		},
		compute_rr_jackknife_log_se = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			if (private$mark_jackknife_nonestimable_if_block_unsupported(unit = unit)) return(NA_real_)
			private$assert_jackknife_supported(unit = unit)
			jack = as.numeric(private$approximate_jackknife_distribution_beta_hat_T_private(unit = unit))
			jack = jack[is.finite(jack) & jack > 0]
			n_units = length(jack)
			if (n_units <= 1L) {
				private$cache_nonestimable_se("jackknife_too_few_positive_risk_ratio_estimates")
				return(NA_real_)
			}
			log_jack = log(jack)
			var_j = ((n_units - 1) / n_units) * sum((log_jack - mean(log_jack))^2)
			if (!is.finite(var_j) || var_j <= 0) {
				private$cache_nonestimable_se("jackknife_log_risk_ratio_standard_error_unavailable")
				return(NA_real_)
			}
			sqrt(var_j)
		},
		compute_rr_jackknife_wald_two_sided_pval = function(delta = 1, unit = "auto"){
			if (!is.finite(delta) || delta <= 0) {
				private$cache_nonestimable_se("jackknife_log_risk_ratio_null_unavailable")
				return(NA_real_)
			}
			est = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
			if (!is.finite(est) || est <= 0) {
				private$cache_nonestimable_estimate("jackknife_original_risk_ratio_unavailable")
				return(NA_real_)
			}
			se_log = private$compute_rr_jackknife_log_se(unit = unit)
			if (!is.finite(se_log) || se_log <= 0) return(NA_real_)
			2 * stats::pnorm(-abs((log(est) - log(delta)) / se_log))
		},
		compute_rr_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			est = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
			if (!is.finite(est) || est <= 0) {
				private$cache_nonestimable_estimate("jackknife_original_risk_ratio_unavailable")
				return(ci)
			}
			se_log = private$compute_rr_jackknife_log_se(unit = unit)
			if (!is.finite(se_log) || se_log <= 0) return(ci)
			z = stats::qnorm(1 - alpha / 2)
			ci[] = exp(log(est) + c(-1, 1) * z * se_log)
			ci
		},
		set_failed_fit_cache = function(inference_ready = TRUE){
			private$cached_values = gcomp_cache_failed_standardized_effects(
				private$cached_values,
				inference_ready = inference_ready
			)
		},
		effects_are_usable = function(effects, estimate_only = FALSE){
			estimand = private$get_estimand_type()
			if (estimate_only) {
				if (identical(estimand, "RD")){
					return(is.finite(effects$rd))
				} else {
					return(is.finite(effects$rr) && effects$rr > 0)
				}
			}
			if (identical(estimand, "RD")){
				is.finite(effects$rd) && is.finite(effects$se_rd) && effects$se_rd > 0
			} else {
				is.finite(effects$rr) && effects$rr > 0 &&
					is.finite(effects$log_rr) &&
					is.finite(effects$se_log_rr) && effects$se_log_rr > 0
			}
		},
		weighted_gcomp_fit = function(X_full, row_weights){
			X_curr = X_full
			repeat {
				reduced = private$reduce_design_matrix_preserving_treatment(X_curr)
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
					return(NULL)
				}
				ok = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
				if (sum(ok) <= ncol(X_fit)) return(NULL)
				p_fit = ncol(X_fit)
				boot_ws = if (!is.null(private$gcomp_boot_beta) && length(private$gcomp_boot_beta) == p_fit) {
					private$gcomp_boot_beta
				} else {
					private$get_fit_warm_start_for_length("beta", p_fit)
				}
				mod = tryCatch(
					fast_logistic_regression_weighted_cpp(
						X = X_fit[ok, , drop = FALSE],
						y = as.numeric(private$y[ok]),
						weights = as.numeric(row_weights[ok]),
						warm_start_beta = boot_ws,
						warm_start_fisher_info = private$get_fit_warm_start_fisher(p_fit)
					),
					error = function(e) NULL
				)
				if (is.null(mod)) {
					private$gcomp_boot_beta = NULL
					return(NULL)
				}
				coef_hat = as.numeric(mod$b)
				if (private$coefficients_are_usable(coef_hat)) {
					private$set_fit_warm_start(coef_hat, "beta", fisher = mod$fisher_information)
					private$gcomp_boot_beta = coef_hat
					names(coef_hat) = colnames(X_fit)
					return(list(X = X_fit, j_treat = j_treat, coefficients = coef_hat, estimate_only = TRUE))
				}
				if (ncol(X_curr) <= 2L) return(NULL)
				drop_col = gcomp_select_covariate_to_drop(X_curr, coef_hat)
				if (!is.finite(drop_col)) return(NULL)
				X_curr = X_curr[, -drop_col, drop = FALSE]
			}
		},
		weighted_gcomp_effects_from_row_weights = function(row_weights){
			X_full = gcomp_normalize_treatment_design_matrix(
				private$build_design_matrix(),
				covariate_names = private$get_covariate_names
			)
			fit = private$weighted_gcomp_fit(X_full, row_weights)
			if (is.null(fit) && private$harden && ncol(X_full) > 2L) {
				fit = private$weighted_gcomp_fit(X_full[, 1:2, drop = FALSE], row_weights)
			}
			if (is.null(fit)) return(NULL)
			effects = private$compute_standardized_effects_r(fit)
			if (!private$effects_are_usable(effects, estimate_only = TRUE)) return(NULL)
			effects
		},
		coefficients_are_usable = function(coef_hat){
			length(coef_hat) > 0L &&
				all(is.finite(coef_hat)) &&
				max(abs(coef_hat), na.rm = TRUE) <= private$max_abs_reasonable_coef
		},
		fit_logistic_with_sandwich = function(X_full, estimate_only = FALSE){
			X_curr = X_full
			repeat {
				reduced = private$reduce_design_matrix_preserving_treatment(X_curr)
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
					return(NULL)
				}
				mod = tryCatch(
					fast_logistic_regression_cpp(
						X = X_fit, 
						y = as.numeric(private$y),
						warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
					),
					error = function(e) NULL
				)
				if (is.null(mod)){
					return(NULL)
				}
				coef_hat = as.numeric(mod$b)
				converged = private$coefficients_are_usable(coef_hat)
				if (!converged){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = gcomp_select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}
				private$set_fit_warm_start(coef_hat, "beta", fisher = mod$fisher_information)
				
				if (estimate_only){
					return(list(
						X = X_fit,
						j_treat = j_treat,
						coefficients = coef_hat,
						estimate_only = TRUE
					))
				}
				mu_hat = inv_logit(X_fit %*% coef_hat)
				mu_hat = pmin(pmax(as.numeric(mu_hat), private$prob_clip_eps), 1 - private$prob_clip_eps)
				W = mu_hat * (1 - mu_hat)
				if (any(!is.finite(W)) || any(W <= 0)){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = gcomp_select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}
				post_fit = tryCatch(
					gcomp_logistic_post_fit_cpp(
						X_fit = X_fit,
						y = as.numeric(private$y),
						coef_hat = coef_hat,
						mu_hat = mu_hat,
						j_treat = j_treat
					),
					error = function(e) NULL
				)
				if (is.null(post_fit)){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = gcomp_select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}
				coef_names = colnames(X_fit)
				names(coef_hat) = coef_names
				vcov_robust = post_fit$vcov
				colnames(vcov_robust) = rownames(vcov_robust) = coef_names
				return(list(
					X = X_fit,
					j_treat = j_treat,
					coefficients = coef_hat,
					vcov = vcov_robust,
					mu_hat = mu_hat,
					post_fit = post_fit,
					estimate_only = FALSE
				))
			}
		},
		compute_standardized_effects_r = function(fit){
			X_fit = fit$X
			coef_hat = fit$coefficients
			vcov_robust = fit$vcov
			j_treat = fit$j_treat
			estimate_only = isTRUE(fit$estimate_only)
			potential_outcomes = gcomp_logistic_potential_outcomes(X_fit, coef_hat, j_treat)
			X1 = potential_outcomes$X1
			X0 = potential_outcomes$X0
			risk1_i = potential_outcomes$risk1_i
			risk0_i = potential_outcomes$risk0_i
			risk1 = potential_outcomes$risk1
			risk0 = potential_outcomes$risk0
			rd = risk1 - risk0
			log_rr = if (risk1 > 0 && risk0 > 0) log(risk1) - log(risk0) else NA_real_
			rr = if (is.finite(log_rr)) exp(log_rr) else NA_real_
			if (estimate_only){
				return(list(
					risk1 = risk1,
					risk0 = risk0,
					rd = rd,
					se_rd = NA_real_,
					log_rr = log_rr,
					rr = rr,
					se_log_rr = NA_real_,
					full_coefficients = coef_hat,
					full_vcov = NULL,
					summary_table = NULL
				))
			}
			grad1 = as.numeric(crossprod(X1, risk1_i * (1 - risk1_i))) / nrow(X1)
			grad0 = as.numeric(crossprod(X0, risk0_i * (1 - risk0_i))) / nrow(X0)
			grad_rd = grad1 - grad0
			var_rd = as.numeric(t(grad_rd) %*% vcov_robust %*% grad_rd)
			grad_log_rr = if (risk1 > 0 && risk0 > 0) grad1 / risk1 - grad0 / risk0 else rep(NA_real_, length(grad1))
			var_log_rr = if (all(is.finite(grad_log_rr))) as.numeric(t(grad_log_rr) %*% vcov_robust %*% grad_log_rr) else NA_real_
			std_err = sqrt(pmax(diag(vcov_robust), 0))
			z_vals = coef_hat / std_err
			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)
			list(
				risk1 = risk1,
				risk0 = risk0,
				rd = rd,
				se_rd = if (is.finite(var_rd) && var_rd >= 0) sqrt(var_rd) else NA_real_,
				log_rr = log_rr,
				rr = rr,
				se_log_rr = if (is.finite(var_log_rr) && var_log_rr >= 0) sqrt(var_log_rr) else NA_real_,
				full_coefficients = coef_hat,
				full_vcov = vcov_robust,
				summary_table = summary_table
			)
		},
		compute_standardized_effects = function(fit){
			if (isTRUE(fit$estimate_only)){
				return(private$compute_standardized_effects_r(fit))
			}
			X_fit = fit$X
			coef_hat = fit$coefficients
			j_treat = fit$j_treat
			fast = fit$post_fit
			if (is.null(fast)){
				return(private$compute_standardized_effects_r(fit))
			}
			vcov_robust = fast$vcov
			colnames(vcov_robust) = rownames(vcov_robust) = names(coef_hat)
			std_err = fast$std_err
			names(std_err) = names(coef_hat)
			z_vals = fast$z_vals
			names(z_vals) = names(coef_hat)
			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)
			list(
				risk1 = fast$risk1,
				risk0 = fast$risk0,
				rd = fast$rd,
				se_rd = fast$se_rd,
				log_rr = fast$log_rr,
				rr = fast$rr,
				se_log_rr = fast$se_log_rr,
				full_coefficients = coef_hat,
				full_vcov = vcov_robust,
				summary_table = summary_table
			)
		},
		get_effect_estimate = function(){
			estimand = private$get_estimand_type()
			if (identical(estimand, "RD")) return(private$cached_values$rd)
			private$cached_values$rr
		},
		compute_effect_confidence_interval = function(alpha){
			z = stats::qnorm(1 - alpha / 2)
			estimand = private$get_estimand_type()
			if (identical(estimand, "RD")){
				est = private$cached_values$rd
				se = private$cached_values$se_rd
				if (!is.finite(est) || !is.finite(se) || se <= 0){
					return(c(NA_real_, NA_real_))
				}
				ci = est + c(-1, 1) * z * se
			} else {
				log_rr = private$cached_values$log_rr
				se_log_rr = private$cached_values$se_log_rr
				if (!is.finite(log_rr) || !is.finite(se_log_rr) || se_log_rr <= 0){
					return(c(NA_real_, NA_real_))
				}
				ci_log = log_rr + c(-1, 1) * z * se_log_rr
				if (!all(is.finite(exp(ci_log)))){
					stop("G-computation RR: could not compute a finite delta-method confidence interval.")
				}
				ci = exp(ci_log)
			}
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		compute_effect_pvalue = function(delta){
			estimand = private$get_estimand_type()
			if (is.null(delta)){
				delta = private$default_null_value()
			}
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
			}
			if (identical(estimand, "RD")){
				est = private$cached_values$rd
				se = private$cached_values$se_rd
				if (!is.finite(est) || !is.finite(se) || se <= 0){
					return(NA_real_)
				}
				z_stat = (est - delta) / se
			} else {
				log_rr = private$cached_values$log_rr
				se_log_rr = private$cached_values$se_log_rr
				if (should_run_asserts()) {
					if (delta <= 0){
						stop("For RR inference, delta must be strictly positive.")
					}
				}
				if (!is.finite(log_rr) || !is.finite(se_log_rr) || se_log_rr <= 0){
					return(NA_real_)
				}
				z_stat = (log_rr - log(delta)) / se_log_rr
			}
			2 * stats::pnorm(-abs(z_stat))
		},
		shared = function(estimate_only = FALSE){
			if (gcomp_standardized_effect_cache_is_ready(private$cached_values, estimate_only = estimate_only)) return(invisible(NULL))
			if (estimate_only && isFALSE(private$harden)) {
				X = private$build_design_matrix()
				if (!is.null(dim(X)) && is.null(colnames(X))) {
					colnames(X) = c("(Intercept)", "treatment", if (ncol(X) > 2L) private$get_covariate_names() else NULL)
				}
				b = glm.fit(x = X, y = as.numeric(private$y), family = binomial())$coefficients
				if (length(b) >= 2L && all(is.finite(b))) {
					eta_base = drop(X %*% b) - X[, 2L] * b[2L]
					risk1 = mean(plogis(eta_base + b[2L]))
					risk0 = mean(plogis(eta_base))
					private$best_X_colnames = setdiff(colnames(X), c("(Intercept)", "treatment"))
					private$cached_values$risk1 = risk1
					private$cached_values$risk0 = risk0
					private$cached_values$rd = as.numeric(risk1 - risk0)
					private$cached_values$rr = as.numeric(risk1 / risk0)
					private$cached_values$log_rr = log(private$cached_values$rr)
					private$cached_values$beta_hat_T = if (identical(private$get_estimand_type(), "RD")) private$cached_values$rd else private$cached_values$rr
					private$cached_values$gcomp_standardized_effects_inference_ready = FALSE
					return(invisible(NULL))
				}
			}
			X_full = gcomp_normalize_treatment_design_matrix(
				private$build_design_matrix(),
				covariate_names = private$get_covariate_names
			)
			fit = private$fit_logistic_with_sandwich(X_full, estimate_only = estimate_only)
			if (!is.null(fit)) {
				private$best_X_colnames = setdiff(colnames(fit$X), c("(Intercept)", "treatment"))
			}
			effects = if (!is.null(fit)) private$compute_standardized_effects(fit) else NULL
			if (private$harden && (is.null(fit) || is.null(effects) || !private$effects_are_usable(effects, estimate_only)) && ncol(X_full) > 2L){
				fit = private$fit_logistic_with_sandwich(X_full[, 1:2, drop = FALSE], estimate_only = estimate_only)
				effects = if (!is.null(fit)) private$compute_standardized_effects(fit) else NULL
			}
			if (is.null(fit) || is.null(effects) || !private$effects_are_usable(effects, estimate_only)){
				private$set_failed_fit_cache(inference_ready = !estimate_only)
				return(invisible(NULL))
			}
			private$cached_values = gcomp_cache_standardized_effects(
				private$cached_values,
				effects,
				inference_ready = !estimate_only
			)
			# Populate beta_hat_T for the base class methods
			private$cached_values$beta_hat_T = if (identical(private$get_estimand_type(), "RD")) effects$rd else effects$rr
		}
	)
)
