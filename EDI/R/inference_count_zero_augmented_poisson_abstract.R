#' Zero-Augmented Poisson Inference for Count Responses
#'
#' Internal base class for non-KK zero-inflated and hurdle Poisson regression
#' models fit using the \pkg{glmmTMB} fitter. The reported treatment effect is the
#' treatment coefficient from the conditional count component, on the log-rate
#' scale.
#'
#' @keywords internal
#' @noRd
InferenceCountZeroAugmentedPoissonAbstract = R6::R6Class("InferenceCountZeroAugmentedPoissonAbstract",
	lock_objects = FALSE,
	inherit = InferenceCountLikelihood,
	public = list(
				
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param model_formula_zero Formula for the zero/hurdle auxiliary component.
		#'   If \code{NULL} (default), the auxiliary component uses the same formula as
		#'   the conditional count component.
		#' @param use_rcpp Whether to use Rcpp speedup.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, model_formula_zero = NULL, use_rcpp = TRUE, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFormula(model_formula, null.ok = TRUE)
				assertFormula(model_formula_zero, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (is.null(model_formula_zero)) {
				model_formula_zero = private$model_formula
			}
			if (should_run_asserts()) {
				assertFormula(model_formula_zero, null.ok = FALSE)
			}
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
			
			
			private$use_rcpp = use_rcpp
			private$model_formula_zero = model_formula_zero
		},
		#' @description Compute asymp confidence interval
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) {
				return(private$count_likelihood_missing_ci(alpha))
			}
			private$shared()
			se = private$get_standard_error()
			if (is.finite(se) && se > 0) {
				private$cached_values$s_beta_hat_T = se
			}
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute asymp two sided pval for treatment effect
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$shared()
			se = private$get_standard_error()
			if (is.finite(se) && se > 0) {
				private$cached_values$s_beta_hat_T = se
			}
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		compute_wald_two_sided_pval = function(delta = 0){
			private$compute_zero_augmented_robust_wald_pval(delta)
		},
		compute_wald_confidence_interval = function(alpha = 0.05){
			private$compute_zero_augmented_robust_wald_ci(alpha)
		},
		compute_score_two_sided_pval = function(delta = 0){
			private$compute_zero_augmented_robust_wald_pval(delta)
		},
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			private$compute_zero_augmented_robust_wald_pval(delta)
		},
		compute_gradient_two_sided_pval = function(delta = 0){
			private$compute_zero_augmented_robust_wald_pval(delta)
		},
		compute_score_confidence_interval = function(alpha = 0.05){
			private$compute_zero_augmented_robust_wald_ci(alpha)
		},
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			private$compute_zero_augmented_robust_wald_ci(alpha)
		},
		compute_gradient_confidence_interval = function(alpha = 0.05){
			private$compute_zero_augmented_robust_wald_ci(alpha)
		},
		compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L){
			if (private$zero_augmented_model_lrt_bootstrap_disabled()) {
				private$cache_nonestimable_se("zero_augmented_poisson_parametric_lrt_bootstrap_disabled_due_raw_lrt_miscalibration")
				return(NA_real_)
			}
			super$compute_lik_ratio_bootstrap_two_sided_pval(
				delta = delta,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate
			)
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			if (!check_package_installed("glmmTMB")) {
				stop(class(self)[1], " weighted bootstrap estimation requires package 'glmmTMB'.")
			}
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			X_fit = private$build_component_matrix(private$model_formula, private$best_X_colnames, treatment_name = "w")
			if (is.null(X_fit)) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_design_unusable")
				return(NA_real_)
			}
			Xzi_fit = private$build_component_matrix(private$model_formula_zero, private$best_Xzi_colnames, treatment_name = "w")
			if (is.null(Xzi_fit)) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_aux_design_unusable")
				return(NA_real_)
			}
			dat = private$build_component_frame(X_fit, Xzi_fit)
			mod = private$fit_zero_augmented_model(dat, X_fit, Xzi_fit, weights = row_weights)
			if (is.null(mod)) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_weighted_fit_unavailable")
				return(NA_real_)
			}
			cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
			if (is.null(cond_coef) || !("w" %in% names(cond_coef)) || !is.finite(cond_coef["w"])) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_weighted_treatment_missing")
				return(NA_real_)
			}
			private$clear_nonestimable_state()
			private$cached_mod = mod
			private$cached_values$likelihood_test_context = NULL
			private$cached_values$beta_hat_T = as.numeric(cond_coef["w"])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$full_coefficients = cond_coef
			private$cached_values$beta_hat_T
		},
		#' @description Zero-augmented mixture fits are too unstable under delete-one
		#'   refits for jackknife inference; report explicit non-estimability.
		compute_jackknife_estimate = function(unit = "auto"){
			private$cache_nonestimable_se("zero_augmented_poisson_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_corrected_estimate = function(unit = "auto"){
			self$compute_jackknife_estimate(unit = unit)
		},
		compute_jackknife_bias_estimate = function(unit = "auto"){
			private$cache_nonestimable_se("zero_augmented_poisson_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_std_error = function(unit = "auto"){
			private$cache_nonestimable_se("zero_augmented_poisson_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_standard_error = function(unit = "auto"){
			self$compute_jackknife_std_error(unit = unit)
		},
		compute_jackknife_wald_two_sided_pval = function(delta = 0, unit = "auto"){
			private$cache_nonestimable_se("zero_augmented_poisson_jackknife_not_supported")
			NA_real_
		},
		compute_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			private$cache_nonestimable_se("zero_augmented_poisson_jackknife_not_supported")
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L){
			if (!is.null(type) && identical(tolower(type), "bca")) {
				private$cache_nonestimable_se("zero_augmented_poisson_jackknife_not_supported")
				return(NA_real_)
			}
			super$compute_bootstrap_two_sided_pval(delta = delta, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples)
		},
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L){
			if (!is.null(type) && identical(tolower(type), "bca")) {
				return(private$missing_bootstrap_ci(alpha, "zero_augmented_poisson_jackknife_not_supported", stage = "se"))
			}
			super$compute_bootstrap_confidence_interval(alpha = alpha, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples)
		},
		compute_bayesian_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
			if (!is.null(type) && identical(tolower(type), "bca")) {
				private$cache_nonestimable_se("zero_augmented_poisson_jackknife_not_supported")
				return(NA_real_)
			}
			super$compute_bayesian_bootstrap_two_sided_pval(delta = delta, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples, weighting_unit_type = weighting_unit_type)
		},
		compute_bayesian_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
			if (!is.null(type) && identical(tolower(type), "bca")) {
				return(private$missing_bootstrap_ci(alpha, "zero_augmented_poisson_jackknife_not_supported", stage = "se"))
			}
			super$compute_bayesian_bootstrap_confidence_interval(alpha = alpha, B = B, type = type, na.rm = na.rm, show_progress = show_progress, min_number_usable_samples = min_number_usable_samples, weighting_unit_type = weighting_unit_type)
		}
	),
		private = list(
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		# Override the generic design-backed worker to freeze best_X_colnames /
		# best_Xzi_colnames in the worker state.  load_bootstrap_sample_into_design_backed_worker
		# clears best_X_colnames (line 855) which would force a full design-reduction on every
		# BRT draw — saving and restoring it here skips that and goes straight to the C++ fit.
		create_bootstrap_worker_state = function(){
			ws = private$create_design_backed_bootstrap_worker_state()
			ws$base_best_X_colnames   = private$best_X_colnames
			ws$base_best_Xzi_colnames = private$best_Xzi_colnames
			ws
		},
		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
			worker_state$worker_priv$best_X_colnames   = worker_state$base_best_X_colnames
			worker_state$worker_priv$best_Xzi_colnames = worker_state$base_best_Xzi_colnames
		},
		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},
			cached_mod = NULL,
			za_X_cov_all = NULL,
			za_Xzi_cov_all = NULL,
			record_zero_augmented_fit_summary = function(fit, X_full, Xzi_full, X_fit, Xzi_fit, is_hurdle = FALSE, fallback_used = FALSE, fallback_reason = NULL){
				cond_full_names = colnames(X_full)
				aux_full_names = colnames(Xzi_full)
				cond_fit_names = colnames(X_fit)
				aux_fit_names = colnames(Xzi_fit)
				params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
				if (length(params) < length(cond_fit_names) + length(aux_fit_names)) return(invisible(NULL))
				cond_fit = params[seq_along(cond_fit_names)]
				aux_fit = params[length(cond_fit_names) + seq_along(aux_fit_names)]
				
				cond_full = rep(NA_real_, length(cond_full_names))
				names(cond_full) = cond_full_names
				aux_full = rep(NA_real_, length(aux_full_names))
				names(aux_full) = aux_full_names
				cond_full[match(cond_fit_names, cond_full_names)] = cond_fit
				aux_full[match(aux_fit_names, aux_full_names)] = aux_fit
				
				coef_all = c(cond_full, aux_full)
				row_names = c(
					paste0("conditional:", names(cond_full)),
					paste0(if (isTRUE(is_hurdle)) "hurdle:" else "zero:", names(aux_full))
				)
				se_all = rep(NA_real_, length(coef_all))
				vcov_fit = tryCatch({
					vc = fit$vcov %||% solve(as.matrix(fit$fisher_information))
					as.matrix(vc)
				}, error = function(e) NULL)
				if (!is.null(vcov_fit) && nrow(vcov_fit) >= length(cond_fit_names) + length(aux_fit_names)) {
					se_fit = sqrt(pmax(0, diag(vcov_fit)[seq_len(length(cond_fit_names) + length(aux_fit_names))]))
					se_cond = rep(NA_real_, length(cond_full_names))
					se_aux = rep(NA_real_, length(aux_full_names))
					se_cond[match(cond_fit_names, cond_full_names)] = se_fit[seq_along(cond_fit_names)]
					se_aux[match(aux_fit_names, aux_full_names)] = se_fit[length(cond_fit_names) + seq_along(aux_fit_names)]
					se_all = c(se_cond, se_aux)
				}
				
				summary_table = matrix(NA_real_, nrow = length(coef_all), ncol = 4L)
				rownames(summary_table) = row_names
				colnames(summary_table) = c("Value", "Std. Error", "z value", "Pr(>|z|)")
				summary_table[, 1L] = coef_all
				summary_table[, 2L] = se_all
				ok = is.finite(coef_all) & is.finite(se_all) & se_all > 0
				summary_table[ok, 3L] = coef_all[ok] / se_all[ok]
				summary_table[ok, 4L] = 2 * stats::pnorm(-abs(summary_table[ok, 3L]))
				
				private$cached_values$full_coefficients = cond_full
				private$cached_values$zero_coefficients = aux_full
				private$cached_values$summary_table = summary_table
				private$cached_values$model_fit_fallback = NULL
				if (isTRUE(fallback_used)) {
					private$cached_values$model_fit_fallback = list(
						used = TRUE,
						reason = fallback_reason %||% "full_zero_augmented_rcpp_fit_failed_to_converge",
						requested_model = "full covariate-adjusted zero-augmented Rcpp model",
						fitted_model = "treatment-only zero-augmented Rcpp model",
						omitted_conditional = setdiff(cond_full_names, cond_fit_names),
						omitted_auxiliary = setdiff(aux_full_names, aux_fit_names),
						family = private$za_description()
					)
				}
				invisible(NULL)
			},
			invalidate_likelihood_fit = function(reason){
			private$cached_mod = NULL
			private$cached_values$likelihood_test_context = NULL
			private$cache_nonestimable_estimate(reason)
			invisible(NULL)
		},
		get_standard_error = function(){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$shared(estimate_only = FALSE)
			se_cached = private$cached_values$s_beta_hat_T
			if (is.finite(se_cached) && se_cached > 0) return(se_cached)
			se = private$compute_standard_error_from_information_matrix()
			if (is.finite(se)) return(se)
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df %||% NA_real_
		},
		compute_zero_augmented_robust_wald_pval = function(delta = 0){
			if (should_run_asserts()) assertNumeric(delta, len = 1)
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$shared(estimate_only = FALSE)
			se = private$cached_values$s_beta_hat_T
			if (!is.finite(se) || se <= 0) {
				private$cache_nonestimable_se("zero_augmented_poisson_robust_standard_error_unavailable")
				return(NA_real_)
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		compute_zero_augmented_robust_wald_ci = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$mark_count_likelihood_block_asymp_nonestimable()) {
				return(private$count_likelihood_missing_ci(alpha))
			}
			private$shared(estimate_only = FALSE)
			se = private$cached_values$s_beta_hat_T
			if (!is.finite(se) || se <= 0) {
				private$cache_nonestimable_se("zero_augmented_poisson_robust_standard_error_unavailable")
				return(private$count_likelihood_missing_ci(alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		zero_augmented_model_lrt_bootstrap_disabled = function(){
			identical(private$za_description(), "Zero-Inflated Poisson") ||
				identical(private$za_description(), "Hurdle Poisson")
		},
		safe_zero_augmented_vcov_se = function(fit, j_treat = 2L){
			v = tryCatch({
				vc = fit$vcov
				if (is.null(vc) || length(dim(vc)) != 2L ||
						nrow(vc) < j_treat || ncol(vc) < j_treat) {
					NA_real_
				} else {
					as.numeric(vc[j_treat, j_treat])
				}
			}, error = function(e) NA_real_)
			if (is.finite(v) && v > 0) sqrt(v) else NA_real_
		},
		best_X_colnames = NULL,
		best_Xzi_colnames = NULL,
		use_rcpp = TRUE,
		model_formula_zero = NULL,
		get_complexity_tier = function() "light",
		build_component_matrix = function(model_formula, selected_colnames = NULL, treatment_name = "treatment"){
			# Determine if this is the zero-inflation formula or conditional formula
			is_zero_formula = identical(model_formula, private$model_formula_zero)
			
			if (is_zero_formula) {
				if (is.null(private$za_Xzi_cov_all)) {
					if (identical(model_formula, ~ .)) {
						X_cov_all = private$get_X()
					} else {
						X_imp = private$des_obj$get_X_imp()
						X_cov_all = if (is.null(X_imp)) matrix(NA_real_, nrow = private$n, ncol = 0) else create_model_matrix_from_features(model_formula, X_imp)
					}
					private$za_Xzi_cov_all = X_cov_all
				} else {
					X_cov_all = private$za_Xzi_cov_all
				}
			} else {
				if (is.null(private$za_X_cov_all)) {
					if (identical(model_formula, ~ .)) {
						X_cov_all = private$get_X()
					} else {
						X_imp = private$des_obj$get_X_imp()
						X_cov_all = if (is.null(X_imp)) matrix(NA_real_, nrow = private$n, ncol = 0) else create_model_matrix_from_features(model_formula, X_imp)
					}
					private$za_X_cov_all = X_cov_all
				} else {
					X_cov_all = private$za_X_cov_all
				}
			}

			if (is.null(selected_colnames)) {
				X_cov = X_cov_all
				if (is.null(X_cov) || ncol(as.matrix(X_cov)) == 0L) {
					X_fit = cbind(1, private$w)
					colnames(X_fit) = c("(Intercept)", treatment_name)
					return(X_fit)
				}
				X_cov = as.matrix(X_cov)
				if (isTRUE(private$harden)) {
					X_cov = drop_highly_correlated_cols(X_cov, threshold = 0.999)$M
				}
				X_fit = cbind(1, private$w, X_cov)
				colnames(X_fit)[1:2] = c("(Intercept)", treatment_name)
				if (isTRUE(private$harden)) {
					res = drop_linearly_dependent_cols(X_fit)
					X_fit = res$M
					colnames(X_fit) = c("(Intercept)", treatment_name, colnames(X_cov))[res$js]
				}
				X_fit
			} else {
				if (is.null(X_cov_all) || length(selected_colnames) == 0L) {
					X_cov = matrix(NA_real_, nrow = private$n, ncol = 0)
				} else {
					X_cov = as.matrix(X_cov_all[, intersect(selected_colnames, colnames(X_cov_all)), drop = FALSE])
				}
				if (is.null(X_cov) || ncol(as.matrix(X_cov)) == 0L) {
					X_fit = cbind(1, private$w)
					colnames(X_fit) = c("(Intercept)", treatment_name)
					return(X_fit)
				}
				X_fit = cbind(1, private$w, as.matrix(X_cov))
				colnames(X_fit)[1:2] = c("(Intercept)", treatment_name)
				X_fit
			}
		},
		build_component_frame = function(X_cond, Xzi){
			dat = data.frame(y = private$y, w = private$w)
			if (ncol(X_cond) > 2L) {
				Xc = as.data.frame(X_cond[, -c(1, 2), drop = FALSE])
				names(Xc) = make.names(colnames(X_cond)[-c(1, 2)], unique = TRUE)
				dat = cbind(dat, Xc)
			}
			if (ncol(Xzi) > 2L) {
				Xz = as.data.frame(Xzi[, -c(1, 2), drop = FALSE])
				names(Xz) = make.names(colnames(Xzi)[-c(1, 2)], unique = TRUE)
				for (nm in names(Xz)) {
					if (!nm %in% names(dat)) dat[[nm]] = Xz[[nm]]
				}
			}
			dat
		},
		build_formula_from_matrix = function(X_fit, response = "y"){
			vars = colnames(X_fit)[-1]
			rhs = if (!length(vars)) "1" else paste(vars, collapse = " + ")
			if (is.null(response)) return(stats::as.formula(paste("~", rhs)))
			stats::as.formula(paste(response, "~", rhs))
		},
		zero_augmented_sandwich_se = function(fit, X_fit, Xzi_fit, j_treat = 2L, is_hurdle = FALSE){
			params = as.numeric(fit$params %||% NA_real_)
			bread = tryCatch(as.matrix(fit$vcov), error = function(e) NULL)
			if (is.null(bread) || !length(params) || any(!is.finite(params))) return(NA_real_)
			total_p = ncol(X_fit) + ncol(Xzi_fit)
			j_treat = as.integer(j_treat)[1L]
			if (length(params) != total_p || nrow(bread) != total_p || ncol(bread) != total_p ||
					j_treat < 1L || j_treat > ncol(X_fit) || any(!is.finite(bread))) {
				return(NA_real_)
			}
			X_fit = as.matrix(X_fit)
			Xzi_fit = as.matrix(Xzi_fit)
			y = as.numeric(private$y)
			lambda = exp(pmin(as.numeric(X_fit %*% params[seq_len(ncol(X_fit))]), 700))
			pi = stats::plogis(as.numeric(Xzi_fit %*% params[ncol(X_fit) + seq_len(ncol(Xzi_fit))]))
			score = matrix(0, nrow = nrow(X_fit), ncol = total_p)
			if (isTRUE(is_hurdle)) {
				is_zero = y == 0
				if (any(is_zero)) {
					score[is_zero, ncol(X_fit) + seq_len(ncol(Xzi_fit))] =
						(1 - pi[is_zero]) * Xzi_fit[is_zero, , drop = FALSE]
				}
				if (any(!is_zero)) {
					lambda_pos = lambda[!is_zero]
					denom = pmax(1 - exp(-lambda_pos), 1e-15)
					score[!is_zero, seq_len(ncol(X_fit))] =
						(y[!is_zero] - lambda_pos / denom) * X_fit[!is_zero, , drop = FALSE]
					score[!is_zero, ncol(X_fit) + seq_len(ncol(Xzi_fit))] =
						-pi[!is_zero] * Xzi_fit[!is_zero, , drop = FALSE]
				}
			} else {
				is_zero = y == 0
				if (any(is_zero)) {
					e0 = exp(-lambda[is_zero])
					A = pmax(pi[is_zero] + (1 - pi[is_zero]) * e0, 1e-15)
					score[is_zero, seq_len(ncol(X_fit))] =
						-((1 - pi[is_zero]) * e0 * lambda[is_zero] / A) * X_fit[is_zero, , drop = FALSE]
					score[is_zero, ncol(X_fit) + seq_len(ncol(Xzi_fit))] =
						(pi[is_zero] * (1 - pi[is_zero]) * (1 - e0) / A) * Xzi_fit[is_zero, , drop = FALSE]
				}
				if (any(!is_zero)) {
					score[!is_zero, seq_len(ncol(X_fit))] =
						(y[!is_zero] - lambda[!is_zero]) * X_fit[!is_zero, , drop = FALSE]
					score[!is_zero, ncol(X_fit) + seq_len(ncol(Xzi_fit))] =
						-pi[!is_zero] * Xzi_fit[!is_zero, , drop = FALSE]
				}
			}
			meat = crossprod(score)
				vcov_robust = tryCatch(bread %*% meat %*% bread, error = function(e) NULL)
				if (is.null(vcov_robust) || nrow(vcov_robust) < j_treat) return(NA_real_)
				se = sqrt(as.numeric(vcov_robust[j_treat, j_treat]))
				if (is.finite(se) && se > 0) se else NA_real_
			},
			hurdle_poisson_lambda_mle = function(mean_positive){
				if (!is.finite(mean_positive) || mean_positive < 1) return(NA_real_)
				if (mean_positive <= 1 + sqrt(.Machine$double.eps)) return(.Machine$double.eps)
				fn = function(lambda) lambda / (-expm1(-lambda)) - mean_positive
				tryCatch(
					stats::uniroot(fn, lower = .Machine$double.eps, upper = max(2, 2 * mean_positive + 1), tol = 1e-10)$root,
					error = function(e) NA_real_
				)
			},
			hurdle_poisson_neg_loglik = function(params, X_fit, Xzi_fit){
				y = as.numeric(private$y)
				eta_c = as.numeric(X_fit %*% params[seq_len(ncol(X_fit))])
				eta_z = as.numeric(Xzi_fit %*% params[ncol(X_fit) + seq_len(ncol(Xzi_fit))])
				lambda = exp(pmin(eta_c, 700))
				pi0 = stats::plogis(eta_z)
				is_zero = y == 0
				nll = numeric(length(y))
				nll[is_zero] = -log(pmax(pi0[is_zero], 1e-15))
				if (any(!is_zero)) {
					log1m_pi = stats::plogis(eta_z[!is_zero], lower.tail = FALSE, log.p = TRUE)
					log1m_exp_minus_lambda = log(-expm1(-lambda[!is_zero]))
					nll[!is_zero] = -(log1m_pi + y[!is_zero] * eta_c[!is_zero] - lambda[!is_zero] - log1m_exp_minus_lambda)
				}
				sum(nll)
			},
			fit_treatment_only_hurdle_poisson_closed_form = function(X_fit, Xzi_fit, estimate_only = FALSE){
				if (ncol(X_fit) != 2L || ncol(Xzi_fit) != 2L) return(NULL)
				w = as.numeric(X_fit[, 2L])
				y = as.numeric(private$y)
				if (length(unique(w[is.finite(w)])) != 2L || any(!(w %in% c(0, 1))) || any(!is.finite(y))) return(NULL)
				n0 = sum(w == 0)
				n1 = sum(w == 1)
				if (n0 == 0L || n1 == 0L) return(NULL)
				y0_pos = y[w == 0 & y > 0]
				y1_pos = y[w == 1 & y > 0]
				if (length(y0_pos) == 0L || length(y1_pos) == 0L) return(NULL)
				
				lambda0 = private$hurdle_poisson_lambda_mle(mean(y0_pos))
				lambda1 = private$hurdle_poisson_lambda_mle(mean(y1_pos))
				if (!is.finite(lambda0) || !is.finite(lambda1) || lambda0 <= 0 || lambda1 <= 0) return(NULL)
				
				p0 = (sum(w == 0 & y == 0) + 0.5) / (n0 + 1)
				p1 = (sum(w == 1 & y == 0) + 0.5) / (n1 + 1)
				params = c(
					log(lambda0),
					log(lambda1) - log(lambda0),
					stats::qlogis(p0),
					stats::qlogis(p1) - stats::qlogis(p0)
				)
				if (!all(is.finite(params))) return(NULL)
				names(params) = c(colnames(X_fit), colnames(Xzi_fit))
				neg_loglik = private$hurdle_poisson_neg_loglik(params, X_fit, Xzi_fit)
				observed_information = tryCatch(
					as.matrix(get_zero_augmented_poisson_hessian_cpp(X_fit, y, Xzi_fit, params, is_hurdle = TRUE)),
					error = function(e) NULL
				)
				vcov = NULL
				if (!estimate_only && !is.null(observed_information) &&
						nrow(observed_information) == length(params) &&
						all(is.finite(observed_information))) {
					vcov = tryCatch(solve(observed_information), error = function(e) NULL)
					if (is.null(vcov) || !all(is.finite(vcov))) vcov = NULL
				}
				if (is.null(vcov)) {
					vcov = matrix(NA_real_, nrow = length(params), ncol = length(params))
				}
				colnames(vcov) = rownames(vcov) = names(params)
				
				list(
					coefficients = list(
						cond = stats::setNames(params[seq_len(ncol(X_fit))], colnames(X_fit)),
						zi = stats::setNames(params[ncol(X_fit) + seq_len(ncol(Xzi_fit))], colnames(Xzi_fit))
					),
					params = params,
					vcov = vcov,
					converged = TRUE,
					neg_ll = neg_loglik,
					neg_loglik = neg_loglik,
					observed_information = observed_information,
					fisher_information = observed_information,
					information = observed_information,
					information_type = "observed",
					hessian = if (!is.null(observed_information)) -observed_information else NULL,
					closed_form = TRUE
				)
			},
			compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_fit = private$build_component_matrix(private$model_formula, private$best_X_colnames, treatment_name = "treatment")
			Xzi_fit = private$build_component_matrix(private$model_formula_zero, private$best_Xzi_colnames, treatment_name = "treatment")
			if (private$use_rcpp && identical(private$za_description(), "Zero-Inflated Negative Binomial")) {
				n_params = ncol(X_fit) + ncol(Xzi_fit) + 1L
				ws_args = private$get_backend_warm_start_args(n_params)
				vc_start = ncol(X_fit) + 1L
				n_vc = ncol(Xzi_fit) + 1L
				has_vc = !is.null(private$cached_vc_params) && length(private$cached_vc_params) == n_vc && all(is.finite(private$cached_vc_params))
				fit = tryCatch(
					fast_zinb_cpp(
						X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg,
						fixed_idx    = if (has_vc) as.integer(vc_start:(vc_start + n_vc - 1L)) else NULL,
						fixed_values = if (has_vc) as.numeric(private$cached_vc_params) else NULL
					),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				private$set_fit_warm_start(as.numeric(fit$params), "params")
				return(as.numeric(fit$params[2]))
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				n_params = ncol(X_fit) + ncol(Xzi_fit)
				ws_args = private$get_backend_warm_start_args(n_params)
					fit = tryCatch(
						fast_zero_augmented_poisson_cpp(
							X_fit, as.numeric(private$y), Xzi_fit,
						is_hurdle = is_hurdle,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if ((is.null(fit) || !isTRUE(fit$converged)) && is_hurdle) {
						X_fit_fallback = X_fit[, seq_len(min(2L, ncol(X_fit))), drop = FALSE]
						Xzi_fit_fallback = Xzi_fit[, seq_len(min(2L, ncol(Xzi_fit))), drop = FALSE]
						fit_fallback = private$fit_treatment_only_hurdle_poisson_closed_form(
							X_fit_fallback, Xzi_fit_fallback, estimate_only = estimate_only
						)
						if (!is.null(fit_fallback) && isTRUE(fit_fallback$converged)) {
							fit = fit_fallback
						}
					}
					if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
					private$set_fit_warm_start(as.numeric(fit$params), "params")
					return(as.numeric(fit$params[2]))
			} else {
				dat = private$build_component_frame(X_fit, Xzi_fit)
				mod = private$fit_zero_augmented_model(dat, X_fit, Xzi_fit)
				if (is.null(mod)) return(NA_real_)
				
				cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (is.null(cond_coef) || !("w" %in% names(cond_coef))) return(NA_real_)
				return(as.numeric(cond_coef["w"]))
			}
		},
		za_family = function() stop(class(self)[1], " must implement za_family()."),
		za_description = function() stop(class(self)[1], " must implement za_description()."),
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp) && private$za_description() %in% c(
				"Zero-Inflated Negative Binomial",
				"Zero-Inflated Poisson",
				"Hurdle Poisson"
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			beta_hat_T = as.numeric(private$cached_values$beta_hat_T %||% NA_real_)[1L]
			if (!is.finite(beta_hat_T)) {
				private$invalidate_likelihood_fit("zero_augmented_poisson_treatment_missing")
				return(NULL)
			}
			X_fit = ctx$X
			Xzi_fit = ctx$Xzi
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			is_hurdle = isTRUE(ctx$is_hurdle)
			is_zinb = identical(private$za_description(), "Zero-Inflated Negative Binomial")
			start_len = if (is_zinb) ncol(X_fit) + ncol(Xzi_fit) + 1L else ncol(X_fit) + ncol(Xzi_fit)
			list(
				X = X_fit,
				Xzi = Xzi_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm_start_params = start %||% private$get_fit_warm_start_for_length("params", start_len)
					warm_fisher = private$get_fit_warm_start_fisher(start_len)
					fit = if (is_zinb) {
						fast_zinb_cpp(
							X = X_fit,
							y = y,
							Xzi = Xzi_fit,
							warm_start_params = warm_start_params,
							warm_start_fisher_info = warm_fisher,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					} else {
						fast_zero_augmented_poisson_cpp(
							X_fit,
							y,
							Xzi_fit,
							is_hurdle = is_hurdle,
							warm_start_params = warm_start_params,
							warm_start_fisher_info = warm_fisher,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					}
					if (!is.null(fit) && is.null(fit$b) && is.null(fit$params) && !is.null(fit$coefficients)) {
						fit$params = as.numeric(c(fit$coefficients$cond, fit$coefficients$zi))
					}
					fit
				},
				extract_start = function(fit){
					as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
				},
				score = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb) {
						as.numeric(fit$score %||% get_zinb_score_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params))
					} else {
						as.numeric(fit$score %||% get_zero_augmented_poisson_score_cpp(X_fit, y, Xzi_fit, params, is_hurdle = is_hurdle))
					}
				},
				observed_information = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					as.matrix(fit$observed_information %||% (if (is_zinb) -get_zinb_hessian_cpp(X_fit, y, Xzi_fit, params) else -get_zero_augmented_poisson_hessian_cpp(X_fit, y, Xzi_fit, params, is_hurdle)))
				},
				information = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					as.matrix(fit$information %||% fit$fisher_information %||% fit$observed_information %||% (if (is_zinb) -get_zinb_hessian_cpp(X_fit, y, Xzi_fit, params) else -get_zero_augmented_poisson_hessian_cpp(X_fit, y, Xzi_fit, params, is_hurdle)))
				},
				neg_loglik = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb) {
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_zinb_neg_loglik_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params))
					} else {
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% fit$mod$neg_ll %||% fit$mod$neg_loglik)
					}
				}
			)
		},
		predictors_df = function(){
			if (ncol(as.matrix(private$X)) > 0){
				full_X = private$create_design_matrix()
				X_model = full_X[, -1, drop = FALSE]
				colnames(X_model)[1] = "w"
				as.data.frame(X_model)
			} else {
				data.frame(w = private$w)
			}
		},
		fit_zero_augmented_model = function(dat, X_fit, Xzi_fit, weights = NULL){
			formula_cond = private$build_formula_from_matrix(X_fit)
			formula_zi = private$build_formula_from_matrix(Xzi_fit, response = NULL)
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)
			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = formula_zi,
						family = private$za_family(),
						data = dat,
						weights = weights,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
			if (!is.null(mod)) return(mod)
			if (ncol(dat) <= 2L) return(NULL)
			dat_fallback = dat[, c("y", "w"), drop = FALSE]
			tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						y ~ w,
						ziformula = ~ w,
						family = private$za_family(),
						data = dat_fallback,
						weights = weights,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
		},
			generate_mod = function(estimate_only = FALSE){
				private$cached_values$likelihood_test_context = NULL
				private$cached_values$model_fit_fallback = NULL
				private$cached_values$summary_table = NULL
				X_full = private$build_component_matrix(private$model_formula, treatment_name = "w")
				if (is.null(X_full)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_design_unusable")
					return(NULL)
				}
				if (is.null(private$best_X_colnames)) {
					res_reduced = private$reduce_design_matrix_preserving_treatment(X_full)
					X_fit = res_reduced$X
					if (is.null(X_fit)){
						private$cache_nonestimable_estimate("zero_augmented_poisson_design_unusable")
					return(NULL)
				}
				colnames(X_fit) = colnames(X_full)[res_reduced$keep]
				private$best_X_colnames = setdiff(colnames(X_fit), c("(Intercept)", "w"))
				} else {
					X_fit = private$build_component_matrix(private$model_formula, private$best_X_colnames, treatment_name = "w")
				}
				Xzi_full = private$build_component_matrix(private$model_formula_zero, treatment_name = "w")
				if (is.null(Xzi_full)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_aux_design_unusable")
					return(NULL)
				}
				if (is.null(private$best_Xzi_colnames)) {
					res_zi = private$reduce_design_matrix_preserving_treatment(Xzi_full)
					Xzi_fit = res_zi$X
					if (is.null(Xzi_fit)){
						private$cache_nonestimable_estimate("zero_augmented_poisson_aux_design_unusable")
					return(NULL)
				}
				colnames(Xzi_fit) = colnames(Xzi_full)[res_zi$keep]
				private$best_Xzi_colnames = setdiff(colnames(Xzi_fit), c("(Intercept)", "w"))
				} else {
					Xzi_fit = private$build_component_matrix(private$model_formula_zero, private$best_Xzi_colnames, treatment_name = "w")
				}
				fallback_used = FALSE
				fallback_reason = NULL
				
				out = list()
			if (private$use_rcpp && identical(private$za_description(), "Zero-Inflated Negative Binomial")) {
				n_params = ncol(X_fit) + ncol(Xzi_fit) + 1L
				ws_args = private$get_backend_warm_start_args(n_params)
				fit = tryCatch(
					fast_zinb_cpp(
						X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					),
					error = function(e) NULL
				)
					if ((is.null(fit) || !isTRUE(fit$converged)) && (ncol(X_fit) > 2L || ncol(Xzi_fit) > 2L)) {
						X_fit_fallback = X_fit[, seq_len(min(2L, ncol(X_fit))), drop = FALSE]
						Xzi_fit_fallback = Xzi_fit[, seq_len(min(2L, ncol(Xzi_fit))), drop = FALSE]
					if (ncol(X_fit_fallback) == 2L && ncol(Xzi_fit_fallback) == 2L) {
						n_params_fallback = ncol(X_fit_fallback) + ncol(Xzi_fit_fallback) + 1L
						ws_args_fallback = private$get_backend_warm_start_args(n_params_fallback)
						fit_fallback = tryCatch(
							fast_zinb_cpp(
								X = X_fit_fallback, y = as.numeric(private$y), Xzi = Xzi_fit_fallback,
								warm_start_params = ws_args_fallback$start_params,
								warm_start_fisher_info = ws_args_fallback$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = estimate_only, optimization_alg = private$optimization_alg
							),
							error = function(e) NULL
						)
						if (!is.null(fit_fallback) && isTRUE(fit_fallback$converged)) {
							X_fit = X_fit_fallback
							Xzi_fit = Xzi_fit_fallback
							private$best_X_colnames = character(0)
							private$best_Xzi_colnames = character(0)
							fallback_used = TRUE
							fallback_reason = "full_zinb_rcpp_fit_failed_to_converge"
							fit = fit_fallback
						}
					}
				}
				if (is.null(fit) || !isTRUE(fit$converged)) {
					private$cache_nonestimable_estimate("zinb_fit_unavailable")
					return(NULL)
				}
				beta_hat_T = as.numeric(fit$params[2])
				if (!is.finite(beta_hat_T)) {
					private$invalidate_likelihood_fit("zinb_treatment_missing")
					return(NULL)
				}
				
				private$clear_nonestimable_state()
				private$cached_mod = fit
				private$set_fit_warm_start(as.numeric(fit$params), "params")
				vc_vals = as.numeric(fit$params[(ncol(X_fit) + 1L):length(fit$params)])
				if (all(is.finite(vc_vals))) private$cached_vc_params = vc_vals

				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					Xzi = Xzi_fit,
					j_treat = 2L,
					is_hurdle = FALSE
				)
				private$record_zero_augmented_fit_summary(
					fit = fit,
					X_full = X_full,
					Xzi_full = Xzi_full,
					X_fit = X_fit,
					Xzi_fit = Xzi_fit,
					is_hurdle = FALSE,
					fallback_used = fallback_used,
					fallback_reason = fallback_reason
				)
				out$beta_hat_T = beta_hat_T
				if (!estimate_only) {
					se = private$zero_augmented_sandwich_se(fit, X_fit, Xzi_fit, j_treat = 2L, is_hurdle = FALSE)
					if (!is.finite(se) || se <= 0) {
						se = private$safe_zero_augmented_vcov_se(fit, j_treat = 2L)
					}
					out$ssq_b_j = if (is.finite(se) && se > 0) se^2 else NA_real_
				}
				out$params = as.numeric(fit$params)
				out$neg_loglik = fit$neg_loglik %||% fit$neg_ll
				out$fisher_information = fit$fisher_information
				out$mod = fit
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				n_params = ncol(X_fit) + ncol(Xzi_fit)
				ws_args = private$get_backend_warm_start_args(n_params)
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(
						X_fit, as.numeric(private$y), Xzi_fit,
						is_hurdle = is_hurdle,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					),
					error = function(e) NULL
				)
				if ((is.null(fit) || !isTRUE(fit$converged)) && (ncol(X_fit) > 2L || ncol(Xzi_fit) > 2L)) {
					X_fit_fallback = X_fit[, seq_len(min(2L, ncol(X_fit))), drop = FALSE]
					Xzi_fit_fallback = Xzi_fit[, seq_len(min(2L, ncol(Xzi_fit))), drop = FALSE]
					if (ncol(X_fit_fallback) == 2L && ncol(Xzi_fit_fallback) == 2L) {
						n_params_fallback = ncol(X_fit_fallback) + ncol(Xzi_fit_fallback)
						ws_args_fallback = private$get_backend_warm_start_args(n_params_fallback)
						fit_fallback = tryCatch(
							fast_zero_augmented_poisson_cpp(
								X_fit_fallback, as.numeric(private$y), Xzi_fit_fallback,
								is_hurdle = is_hurdle,
								warm_start_params = ws_args_fallback$start_params,
								warm_start_fisher_info = ws_args_fallback$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = estimate_only, optimization_alg = private$optimization_alg
							),
							error = function(e) NULL
						)
						if (!is.null(fit_fallback) && isTRUE(fit_fallback$converged)) {
							X_fit = X_fit_fallback
							Xzi_fit = Xzi_fit_fallback
							private$best_X_colnames = character(0)
							private$best_Xzi_colnames = character(0)
							fallback_used = TRUE
							fallback_reason = if (is_hurdle) "full_hurdle_poisson_rcpp_fit_failed_to_converge" else "full_zip_rcpp_fit_failed_to_converge"
							fit = fit_fallback
							}
						}
					}
					if ((is.null(fit) || !isTRUE(fit$converged)) && is_hurdle) {
						X_fit_fallback = X_fit[, seq_len(min(2L, ncol(X_fit))), drop = FALSE]
						Xzi_fit_fallback = Xzi_fit[, seq_len(min(2L, ncol(Xzi_fit))), drop = FALSE]
						fit_fallback = private$fit_treatment_only_hurdle_poisson_closed_form(
							X_fit_fallback, Xzi_fit_fallback, estimate_only = estimate_only
						)
						if (!is.null(fit_fallback) && isTRUE(fit_fallback$converged)) {
							X_fit = X_fit_fallback
							Xzi_fit = Xzi_fit_fallback
							private$best_X_colnames = character(0)
							private$best_Xzi_colnames = character(0)
							fallback_used = TRUE
							fallback_reason = "hurdle_poisson_treatment_only_closed_form_after_rcpp_nonconvergence"
							fit = fit_fallback
						}
					}
					if (is.null(fit) || !isTRUE(fit$converged)) {
						private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
						return(NULL)
				}
				beta_hat_T = as.numeric(fit$params[2])
				if (!is.finite(beta_hat_T)) {
					private$invalidate_likelihood_fit("zero_augmented_poisson_treatment_missing")
					return(NULL)
				}
				
				private$clear_nonestimable_state()
				private$cached_mod = fit
				full_params = as.numeric(fit$params)
				private$set_fit_warm_start(full_params, "params")
				
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					Xzi = Xzi_fit,
					j_treat = 2L,
					is_hurdle = is_hurdle
				)
				private$record_zero_augmented_fit_summary(
					fit = fit,
					X_full = X_full,
					Xzi_full = Xzi_full,
					X_fit = X_fit,
					Xzi_fit = Xzi_fit,
					is_hurdle = is_hurdle,
					fallback_used = fallback_used,
					fallback_reason = fallback_reason
				)
				out$beta_hat_T = beta_hat_T
				if (!estimate_only) {
					se = private$zero_augmented_sandwich_se(fit, X_fit, Xzi_fit, j_treat = 2L, is_hurdle = is_hurdle)
					if (!is.finite(se) || se <= 0) {
						se = private$safe_zero_augmented_vcov_se(fit, j_treat = 2L)
					}
					out$ssq_b_j = if (is.finite(se) && se > 0) se^2 else NA_real_
				}
				out$params = full_params
				out$fisher_information = fit$fisher_information
				out$mod = fit
			} else {
				dat = private$build_component_frame(X_fit, Xzi_fit)
				mod = private$fit_zero_augmented_model(dat, X_fit, Xzi_fit)
				if (is.null(mod)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
					return(NULL)
				}
				
				private$clear_nonestimable_state()
				private$cached_values$likelihood_test_context = NULL
				cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (is.null(cond_coef) || !("w" %in% names(cond_coef)) || !is.finite(cond_coef["w"])){
					private$cache_nonestimable_estimate("zero_augmented_poisson_treatment_missing")
					return(NULL)
				}
				out$beta_hat_T = as.numeric(cond_coef["w"])
				if (!estimate_only) {
					coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
					se = if (!is.null(coef_table) && ("w" %in% rownames(coef_table))) as.numeric(coef_table["w", "Std. Error"]) else NA_real_
					out$ssq_b_j = if (is.finite(se) && se > 0) se^2 else NA_real_
				}
				out$mod = mod
			}
			out
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		supports_lik_ratio_param_bootstrap = function(){
			isTRUE(private$use_rcpp) && private$za_description() %in% c(
				"Zero-Inflated Negative Binomial",
				"Zero-Inflated Poisson",
				"Hurdle Poisson"
			)
		},
		supports_lik_ratio_param_bootstrap_confidence_interval = function(){
			FALSE
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			is_zinb   = identical(private$za_description(), "Zero-Inflated Negative Binomial")
			is_hurdle = identical(private$za_description(), "Hurdle Poisson")
			X   = spec$X
			Xzi = spec$Xzi
			j   = spec$j
			n   = nrow(X)
			if (is_zinb) {
				p_null = as.numeric(null_fit$params)
				b_cond = p_null[seq_len(ncol(X))]
				b_zi   = p_null[ncol(X) + seq_len(ncol(Xzi))]
			} else {
				b_cond = as.numeric(null_fit$coefficients$cond)
				b_zi   = as.numeric(null_fit$coefficients$zi)
			}
			lambda = exp(pmin(as.numeric(X %*% b_cond), 20))
			pi     = plogis(as.numeric(Xzi %*% b_zi))
			if (is_zinb){
				full_params = as.numeric(null_fit$params)
				theta = exp(min(full_params[length(full_params)], 15))
				if (!is.finite(theta) || theta <= 0) return(NULL)
				u = rbinom(n, 1L, pi)
				counts = rnbinom(n, size = theta, mu = lambda)
				y_sim = as.integer(ifelse(u == 1L, 0L, counts))
			} else if (is_hurdle){
				p0 = ppois(0, lambda)
				u  = rbinom(n, 1L, pi)
				y_pos = as.integer(qpois(p0 + (1 - p0) * runif(n), lambda))
				y_pos = pmax(y_pos, 1L)
				y_sim = as.integer(ifelse(u == 1L, 0L, y_pos))
			} else {
				u = rbinom(n, 1L, pi)
				y_sim = as.integer(ifelse(u == 1L, 0L, rpois(n, lambda)))
			}
			n_params = if (is_zinb) ncol(X) + ncol(Xzi) + 1L else ncol(X) + ncol(Xzi)
			full_res = tryCatch(
				if (is_zinb) {
					fast_zinb_cpp(
						X = X, y = as.numeric(y_sim), Xzi = Xzi,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg
					)
				} else {
					fast_zero_augmented_poisson_cpp(
						X, as.numeric(y_sim), Xzi,
						is_hurdle = is_hurdle,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg
					)
				},
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged)) return(NULL)
			beta_j = if (is_zinb) as.numeric(full_res$params[j]) else as.numeric(full_res$coefficients$cond[j])
			if (!is.finite(beta_j)) return(NULL)
			if (is.null(full_res$params) && !is.null(full_res$coefficients)){
				full_res$params = as.numeric(c(full_res$coefficients$cond, full_res$coefficients$zi))
			}
			list(
				worker_data = list(y = as.numeric(y_sim)),
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					ws = start %||% private$get_fit_warm_start_for_length("params", n_params)
					fit = tryCatch(
						if (is_zinb) {
							fast_zinb_cpp(
								X = X, y = as.numeric(y_sim), Xzi = Xzi,
								warm_start_params = ws,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = FALSE,
								optimization_alg = private$optimization_alg,
								fixed_idx = j, fixed_values = d
							)
						} else {
							fast_zero_augmented_poisson_cpp(
								X, as.numeric(y_sim), Xzi,
								is_hurdle = is_hurdle,
								warm_start_params = ws,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = FALSE,
								optimization_alg = private$optimization_alg,
								fixed_idx = j, fixed_values = d
							)
						},
						error = function(e) NULL
					)
					if (!is.null(fit) && is.null(fit$params) && !is.null(fit$coefficients)){
						fit$params = as.numeric(c(fit$coefficients$cond, fit$coefficients$zi))
					}
					fit
				},
				neg_loglik = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb){
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_zinb_neg_loglik_cpp(X = X, y = as.numeric(y_sim), Xzi = Xzi, params))
					} else {
						as.numeric(fit$neg_loglik %||% fit$neg_ll)
					}
				}
			)
		}
	)
)
