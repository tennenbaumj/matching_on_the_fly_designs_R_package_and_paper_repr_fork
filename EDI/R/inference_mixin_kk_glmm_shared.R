#' GLMM behaviour bundle for KK-design inference classes
#'
#' Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots).
#' Splice into a daughter class via
#' \code{public  = c(InferenceMixinKKGLMMShared$public,  list(...))} and
#' \code{private = c(InferenceMixinKKGLMMShared$private, list(...))}.
#' The daughter must inherit from \code{InferenceAsympLik} and call
#' \code{private$init_kk_glmm_shared(des_obj)} from its \code{initialize}
#' after \code{super$initialize(...)}.  Each daughter is responsible for
#' providing its own \code{shared()} dispatcher that calls either its
#' Rcpp path or \code{private$shared_glmm_tmb()} as the glmmTMB fallback.
#'
#' @keywords internal
#' @noRd
InferenceMixinKKGLMMShared = list(
	public = list(
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				self$compute_estimate(estimate_only = estimate_only)
				beta_hat_T = as.numeric(private$cached_values$beta_hat_T)[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$df = Inf
					private$cached_values$summary_table = NULL
					return(private$cached_values$beta_hat_T)
				}
			}
			result = private$compute_weighted_glmm_bootstrap_estimate(row_weights, estimate_only = estimate_only)
			if (is.list(result)) {
				private$cached_values$beta_hat_T = result$beta
				private$cached_values$s_beta_hat_T = result$se
			} else {
				private$cached_values$beta_hat_T = as.numeric(result)[1L]
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = Inf
			private$cached_values$summary_table = NULL
			private$cached_values$beta_hat_T
		},
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		m = NULL,
		optimization_alg = "lbfgs",
		# Subclasses that provide their own Rcpp fitter set this to TRUE before
		# calling super$initialize() to suppress the glmmTMB package check.
		skip_glmm_pkg_check = FALSE,
		max_abs_reasonable_coef = 1e4,
		kk_glmm_engine = TRUE,
		init_kk_glmm_shared = function(des_obj){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), private$glmm_response_type())
			}
			if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or DesignFixedBinaryMatch.")
			}
			if (inherits(des_obj, "DesignFixedBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			if (identical(private$glmm_response_type(), "proportion")) {
				private$y = .sanitize_proportion_response(private$y, interior = FALSE)
			}
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !isTRUE(private$skip_glmm_pkg_check)) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
				}
			}
		},
		# Default (multivariate): all covariates + treatment.
		# Univariate subclasses override this to return data.frame(w = private$w).
		glmm_predictors_df = function(){
			df = as.data.frame(private$create_design_matrix()[, -1, drop = FALSE])
			# create_design_matrix uses "treatment"; glmmTMB path expects "w"
			if ("treatment" %in% colnames(df))
				colnames(df)[colnames(df) == "treatment"] = "w"
			df
		},
		glmm_predictors_df_candidates = function(){
			predictors_df = private$glmm_predictors_df()
			if (!private$harden || is.null(predictors_df) || ncol(predictors_df) <= 1L){
				return(list(predictors_df))
			}
			X_full = as.matrix(predictors_df)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = match("w", colnames(X_full)),
				fit_fun = function(X_fit){
					private$fit_glmm_on_data(as.data.frame(X_fit), se = TRUE)
				},
				fit_ok = function(mod, X_fit, keep){
					private$.is_usable_glmm_fit(mod, se = TRUE)
				}
			)
			if (is.null(attempt$fit)) return(list(data.frame(w = predictors_df$w)))
			candidates = list(as.data.frame(attempt$X))
			if (ncol(attempt$X) > 1L){
				if (!("w" %in% unlist(lapply(candidates, colnames), use.names = FALSE))){
					candidates[[length(candidates) + 1L]] = data.frame(w = predictors_df$w)
				}
			}
			candidates
		},
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			se = private$compute_standard_error_from_information_matrix()
			if (is.finite(se)) return(se)
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		},
		# glmmTMB fallback implementation (used by daughter shared() dispatchers).
		# Daughters call private$shared_glmm_tmb() rather than super$shared().
		shared_glmm_tmb = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			private$clear_nonestimable_state()
			mod = private$fit_glmm(se = !estimate_only)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_glmm_fit_failed")
				return(invisible(NULL))
			}
			# glmmTMB fixed effects for the conditional model
			beta = glmmTMB::fixef(mod)$cond
			if ("w" %in% names(beta)){
				private$cached_values$beta_hat_T = as.numeric(beta["w"])
			} else {
				private$cached_values$beta_hat_T = NA_real_
			}
			if (estimate_only) return(invisible(NULL))
			coef_table = summary(mod)$coefficients$cond
			if ("w" %in% rownames(coef_table) && "Std. Error" %in% colnames(coef_table)){
				private$cached_values$s_beta_hat_T = as.numeric(coef_table["w", "Std. Error"])
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = Inf
			private$cached_values$summary_table = coef_table
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				return(invisible(NULL))
		},
		fit_glmm_on_data = function(predictors_df, se = TRUE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			if (!check_package_installed("glmmTMB")){
				return(NULL)
			}
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)
			y_dat = private$y
			if (identical(private$glmm_response_type(), "ordinal")) {
				y_dat = ordered(y_dat, levels = sort(unique(y_dat)))
			}
			dat = data.frame(y = y_dat, predictors_df, group_id = factor(group_id))
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			glmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))
			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					glmmTMB::glmmTMB(
						glmm_formula,
						family  = private$glmm_family(),
						data    = dat,
						control = glmm_control,
						se      = se
					)
				)))
				mod
			}, error = function(e) {
				message(paste("GLMM FIT ERROR:", e$message))
				NULL
			})
		},
		fit_weighted_glmm_on_data = function(predictors_df, row_weights, se = FALSE){
			if (!check_package_installed("glmmTMB")){
				return(NULL)
			}
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L) {
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			}
			ok = is.finite(row_weights) & row_weights > 0
			if (identical(private$glmm_response_type(), "ordinal")) {
				ok = ok & is.finite(private$y)
			} else {
				ok = ok & is.finite(as.numeric(private$y))
			}
			if (!any(ok)) return(NULL)
			y_dat = private$y
			if (identical(private$glmm_response_type(), "ordinal")) {
				y_dat = ordered(y_dat, levels = sort(unique(y_dat[ok])))
			}
			dat = data.frame(
				y = y_dat[ok],
				predictors_df[ok, , drop = FALSE],
				group_id = factor(group_id[ok]),
				.bootstrap_weight__ = as.numeric(row_weights[ok])
			)
			fixed_terms = setdiff(colnames(dat), c("y", "group_id", ".bootstrap_weight__"))
			glmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)
			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					glmmTMB::glmmTMB(
						glmm_formula,
						family = private$glmm_family(),
						data = dat,
						weights = .bootstrap_weight__,
						control = glmm_control,
						se = se
					)
				)))
				mod
			}, error = function(e) NULL)
		},
		fit_glmm = function(se = TRUE){
			for (predictors_df in private$glmm_predictors_df_candidates()){
				mod = private$fit_glmm_on_data(predictors_df, se = se)
				if (private$.is_usable_glmm_fit(mod, se)) return(mod)
			}
			NULL
		},
		compute_weighted_glmm_bootstrap_estimate = function(row_weights, estimate_only = TRUE){
			for (predictors_df in private$glmm_predictors_df_candidates()) {
				mod = private$fit_weighted_glmm_on_data(predictors_df, row_weights = row_weights, se = !estimate_only)
				if (!private$.is_usable_glmm_fit(mod, se = FALSE)) next
				beta = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (!is.null(beta) && "w" %in% names(beta) && is.finite(beta["w"])) {
					beta_val = as.numeric(beta["w"])
					if (estimate_only) return(beta_val)
					se = tryCatch({
						ct = summary(mod)$coefficients$cond
						if (!is.null(ct) && "w" %in% rownames(ct)) {
							se_val = suppressWarnings(as.numeric(ct["w", "Std. Error"]))
							if (is.finite(se_val) && se_val > 0) se_val else NA_real_
						} else NA_real_
					}, error = function(e) NA_real_)
					return(list(beta = beta_val, se = se))
				}
			}
			if (estimate_only) NA_real_ else list(beta = NA_real_, se = NA_real_)
		},
		.is_usable_glmm_fit = function(mod, se){
			if (is.null(mod)) return(FALSE)
			beta = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
			if (is.null(beta) || !("w" %in% names(beta)) || !is.finite(beta["w"])) return(FALSE)
			if (!se) return(TRUE)
			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table)) || !("Std. Error" %in% colnames(coef_table))) return(FALSE)
			se_w = suppressWarnings(as.numeric(coef_table["w", "Std. Error"]))
			is.finite(se_w) && se_w > 0
		}
	)
)
