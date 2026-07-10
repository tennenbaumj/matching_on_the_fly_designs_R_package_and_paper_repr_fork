#' GEE behaviour bundle for KK-design inference classes
#'
#' Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots).
#' Splice into a daughter class via
#' \code{public  = c(InferenceMixinKKGEEShared$public,  list(...))} and
#' \code{private = c(InferenceMixinKKGEEShared$private, list(...))}.
#' The daughter must inherit from \code{InferenceAsymp} and call
#' \code{private$init_kk_gee_shared(des_obj, use_rcpp, model_formula)}
#' from its \code{initialize} after \code{super$initialize(...)}.
#'
#' Under \code{harden = TRUE}, multivariate GEE fits preserve the treatment column
#' and retry reduced covariate sets after QR-based rank reduction and
#' correlation-based pruning. Extreme finite coefficients / standard errors are
#' rejected and treated as non-estimable.
#'
#' @keywords internal
InferenceMixinKKGEEShared = list(
	public = list(
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			beta_hat_T = private$fit_weighted_gee_with_fallback(row_weights)
			private$cached_values$beta_hat_T = as.numeric(beta_hat_T)[1L]
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = Inf
			private$cached_values$summary_table = NULL
			private$cached_values$nonestimable = !is.finite(private$cached_values$beta_hat_T)
			private$cached_values$nonestimable_reason = if (is.finite(private$cached_values$beta_hat_T)) NULL else "weighted_gee_estimate_unavailable"
			private$cached_values$nonestimable_stage = if (is.finite(private$cached_values$beta_hat_T)) NULL else "estimate"
			private$cached_values$beta_hat_T
		},
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (private$use_kk_gee_jackknife_wald_calibration()) {
				return(private$compute_kk_gee_jackknife_wald_confidence_interval(alpha = alpha))
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
			if (private$use_kk_gee_jackknife_wald_calibration()) {
				return(private$compute_kk_gee_jackknife_wald_two_sided_pval(delta = delta))
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				if (should_run_asserts()) {
					stop("TO-DO")
				}
				NA_real_
			}
		}
	),
	private = list(
		m = NULL,
		use_rcpp = TRUE,
		max_abs_reasonable_coef = 1e4,
		kk_gee_engine = TRUE,
		get_complexity_tier = function() "medium",
		use_kk_gee_jackknife_wald_calibration = function(){
			identical(private$gee_response_type(), "count")
		},
		compute_kk_gee_jackknife_wald_two_sided_pval = function(delta = 0){
			p = tryCatch(
				self$compute_jackknife_wald_two_sided_pval(delta = delta),
				error = function(e) NA_real_
			)
			if (!is.finite(p) || p < 0 || p > 1) {
				private$cache_nonestimable_se("kk_count_gee_jackknife_wald_calibration_unavailable")
				return(NA_real_)
			}
			p
		},
		compute_kk_gee_jackknife_wald_confidence_interval = function(alpha = 0.05){
			ci = tryCatch(
				self$compute_jackknife_wald_confidence_interval(alpha = alpha),
				error = function(e) c(NA_real_, NA_real_)
			)
			if (length(ci) < 2L || any(!is.finite(ci[1:2])) || ci[1] > ci[2]) {
				private$cache_nonestimable_se("kk_count_gee_jackknife_wald_calibration_unavailable")
				ci = c(NA_real_, NA_real_)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}
			ci
		},
		compute_wald_two_sided_pval_impl = function(delta){
			if (private$use_kk_gee_jackknife_wald_calibration()) {
				return(private$compute_kk_gee_jackknife_wald_two_sided_pval(delta = delta))
			}
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		compute_wald_confidence_interval_impl = function(alpha){
			if (private$use_kk_gee_jackknife_wald_calibration()) {
				return(private$compute_kk_gee_jackknife_wald_confidence_interval(alpha = alpha))
			}
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		gee_warm_start_args = function(expected_length, expected_fisher_dim = expected_length){
			start_beta = private$get_fit_warm_start_for_length("beta", expected_length)
			list(
				start_beta = start_beta,
				warm_start_beta = start_beta,
				start_params = private$get_fit_warm_start_for_length("params", expected_length),
				warm_start_weights = private$get_fit_warm_start_weights(private$n),
				warm_start_fisher_info = private$get_fit_warm_start_fisher(expected_fisher_dim)
			)
		},
		init_kk_gee_shared = function(des_obj, use_rcpp = TRUE, model_formula = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), private$gee_response_type())
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or DesignFixedBinaryMatch.")
				}
			}
			if (inherits(des_obj, "DesignFixedBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			if (identical(private$gee_response_type(), "proportion")) {
				private$y = .sanitize_proportion_response(private$y, interior = FALSE)
			}
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$use_rcpp = use_rcpp
			if (should_run_asserts() && !use_rcpp) {
				if (!check_package_installed("geepack")){
					stop("Package 'geepack' is required for ", class(self)[1], " when use_rcpp=FALSE. Please install it.")
				}
			}
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			mod = private$fit_gee_with_fallback(std_err = FALSE, estimate_only = estimate_only)
			private$extract_gee_treatment_estimate(mod)
		},
		# Default (multivariate): intercept dropped, treatment column named "w".
		# Univariate subclasses override this to return data.frame(w = private$w).
		gee_predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},
		gee_predictors_df_candidates = function(){
			predictors_df = private$gee_predictors_df()
			if (!private$harden || is.null(predictors_df) || ncol(predictors_df) <= 1L){
				return(list(predictors_df))
			}
			normalize_candidate = function(X_fit){
				if (is.null(X_fit)) return(as.data.frame(predictors_df, check.names = FALSE))
				X_df = as.data.frame(X_fit, check.names = FALSE)
				if (ncol(X_df) == 0L || !("w" %in% colnames(X_df))) {
					return(as.data.frame(predictors_df, check.names = FALSE))
				}
				X_df
			}
			X_full = as.matrix(predictors_df)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = match("w", colnames(X_full)),
				fit_fun = function(X_fit) X_fit,
				fit_ok = function(mod, X_fit, keep) TRUE
			)
			candidates = list(as.data.frame(predictors_df, check.names = FALSE))
			keys = paste(colnames(candidates[[1L]]), collapse = "|")
			first_qr = normalize_candidate(attempt$X_fit)
			first_key = paste(colnames(first_qr), collapse = "|")
			if (!(first_key %in% keys)) {
				candidates[[length(candidates) + 1L]] = first_qr
				keys = c(keys, first_key)
			}
			other_idx = setdiff(seq_len(ncol(X_full)), 1L)
			if (length(other_idx) > 0L){
				thresholds = c(0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
				for (thresh in thresholds){
					X_cov = drop_highly_correlated_cols(X_full[, other_idx, drop = FALSE], threshold = thresh)$M
					X_try = cbind(w = X_full[, 1], X_cov)
					attempt_try = private$fit_with_hardened_qr_column_dropping(
						X_full = X_try,
						required_cols = match("w", colnames(X_try)),
						fit_fun = function(X_fit) X_fit,
						fit_ok = function(mod, X_fit, keep) TRUE
					)
					X_try_df = normalize_candidate(attempt_try$X_fit)
					key = paste(colnames(X_try_df), collapse = "|")
					if (!(key %in% keys)){
						candidates[[length(candidates) + 1L]] = X_try_df
						keys = c(keys, key)
					}
				}
			}
			if (!("w" %in% unlist(lapply(candidates, colnames), use.names = FALSE))){
				candidates[[length(candidates) + 1L]] = data.frame(w = predictors_df$w)
			}
			candidates
		},
		gee_treatment_index = function(beta){
			if (is.null(beta) || !length(beta)) return(NA_integer_)
			beta_names = names(beta)
			if (!is.null(beta_names) && ("w" %in% beta_names)) return(match("w", beta_names))
			if (length(beta) >= 2L) return(2L)
			NA_integer_
		},
		gee_coefficients_are_usable = function(beta){
			length(beta) > 0L &&
				all(is.finite(beta)) &&
				max(abs(beta), na.rm = TRUE) <= private$max_abs_reasonable_coef
		},
		extract_gee_treatment_estimate = function(mod){
			if (is.null(mod)) return(NA_real_)
			beta = tryCatch(
				if (inherits(mod, "geeglm")) {
					stats::coef(mod)
				} else if (!is.null(mod$beta)) {
					mod$beta
				} else {
					stats::coef(mod)
				},
				error = function(e) NULL
			)
			if (is.null(beta) || !private$gee_coefficients_are_usable(beta)) return(NA_real_)
			j_treat = private$gee_treatment_index(beta)
			if (!is.finite(j_treat) || is.na(j_treat) || j_treat < 1L || j_treat > length(beta)) return(NA_real_)
			as.numeric(beta[[j_treat]])
		},
		extract_gee_treatment_se = function(mod, j_treat = NA_integer_, coef_table = NULL){
			if (is.null(mod)) return(NA_real_)
			beta = tryCatch(
				if (inherits(mod, "geeglm")) {
					stats::coef(mod)
				} else if (!is.null(mod$beta)) {
					mod$beta
				} else {
					stats::coef(mod)
				},
				error = function(e) NULL
			)
			if (is.na(j_treat) || !is.finite(j_treat)) j_treat = private$gee_treatment_index(beta)
			if (!is.finite(j_treat) || is.na(j_treat) || j_treat < 1L) return(NA_real_)
			if (!inherits(mod, "geeglm")) {
				vc = mod$vcov
				if (!is.null(vc) && is.matrix(vc) && j_treat <= nrow(vc) && j_treat <= ncol(vc)) {
					v = suppressWarnings(as.numeric(vc[j_treat, j_treat]))
					if (is.finite(v) && v > 0) return(sqrt(v))
				}
				return(NA_real_)
			}
			if (!is.null(coef_table)) {
				se_col = intersect(c("Std.err", "Std.error", "Robust S.E."), colnames(coef_table))
				if (length(se_col) > 0L) {
					row_idx = if (!is.null(rownames(coef_table)) && ("w" %in% rownames(coef_table))) match("w", rownames(coef_table)) else j_treat
					if (is.finite(row_idx) && !is.na(row_idx) && row_idx >= 1L && row_idx <= nrow(coef_table)) {
						se = suppressWarnings(as.numeric(coef_table[row_idx, se_col[1]]))
						if (is.finite(se) && se > 0) return(se)
					}
				}
			}
			vc = tryCatch(stats::vcov(mod), error = function(e) NULL)
			if (!is.null(vc) && is.matrix(vc) && j_treat <= nrow(vc) && j_treat <= ncol(vc)) {
				v = suppressWarnings(as.numeric(vc[j_treat, j_treat]))
				if (is.finite(v) && v > 0) return(sqrt(v))
			}
			NA_real_
		},
		# Dispatches to private$shared_gee_dispatch() which each daughter defines.
		shared = function(estimate_only = FALSE) private$shared_gee_dispatch(estimate_only),
		# Default geepack/Rcpp implementation. Daughters that use this path define:
		#   shared_gee_dispatch = function(e) private$shared_gee_default(e)
		# Ordinal daughters define their own shared_gee_dispatch with ordLORgee logic.
		shared_gee_default = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			private$clear_nonestimable_state()
			mod = private$fit_gee_with_fallback(std_err = !estimate_only, estimate_only = estimate_only)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_gee_fit_failed")
				return(invisible(NULL))
			}
			beta = tryCatch(if (inherits(mod, "geeglm")) stats::coef(mod) else mod$beta, error = function(e) NULL)
			j_treat = private$gee_treatment_index(beta)
			private$cached_values$beta_hat_T = if (is.finite(j_treat) && !is.na(j_treat) && j_treat >= 1L && j_treat <= length(beta)) as.numeric(beta[[j_treat]]) else NA_real_

			if (!is.finite(private$cached_values$beta_hat_T)){
				private$cache_nonestimable_estimate("kk_gee_estimate_unavailable")
				return(invisible(NULL))
			}
			private$cached_values$df   = Inf
			if (estimate_only) {
				if (!inherits(mod, "geeglm")) {
					private$set_fit_warm_start(mod$beta, "beta", fisher = mod$fisher_information)
				} else {
					private$set_fit_warm_start(stats::coef(mod), "beta")
				}
				return(invisible(NULL))
			}
			if (inherits(mod, "geeglm")) {
				coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
				private$cached_values$s_beta_hat_T = private$extract_gee_treatment_se(mod, j_treat = j_treat, coef_table = coef_table)
				private$set_fit_warm_start(stats::coef(mod), "beta")
			} else {
				# Rcpp result
				vc = mod$vcov
				private$cached_values$s_beta_hat_T = if (j_treat <= nrow(vc)) sqrt(as.numeric(vc[j_treat, j_treat])) else NA_real_
				private$cached_values$quasi_loglik = mod$quasi_loglik
				private$cached_values$score        = mod$score

				private$set_fit_warm_start(mod$beta, "beta", fisher = mod$fisher_information)
			}
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				private$cache_nonestimable_se("kk_gee_standard_error_unavailable")
				return(invisible(NULL))
			}
			private$clear_nonestimable_state()
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				return(invisible(NULL))
		},
		gee_has_reservoir = function(){
			m_vec = private$m
			if (is.null(m_vec)) return(FALSE)
			m_vec[is.na(m_vec)] = 0L
			any(m_vec == 0L)
		},
		fit_gee_rcpp = function(fit_data, estimate_only = FALSE){
			family_str = private$gee_family_str()
			X_rcpp = cbind(`(Intercept)` = 1, as.matrix(fit_data$dat))
			tryCatch({
				gee_pairs_singletons_cpp(
					X = X_rcpp,
					y = as.numeric(fit_data$y_sorted),
					group_id = as.integer(fit_data$id_sorted),
					family_str = family_str,
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_rcpp)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_rcpp))
				)
			}, error = function(e) NULL)
		},
		build_gee_fit_data = function(include_reservoir = TRUE, predictors_df = NULL, row_weights = NULL){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			keep = if (isTRUE(include_reservoir)) rep(TRUE, length(m_vec)) else m_vec > 0L
			if (!any(keep)) return(NULL)
			m_keep = m_vec[keep]
			group_id = m_keep
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L) {
				max_group = if (all(group_id == 0L)) 0L else max(group_id)
				group_id[reservoir_idx] = max_group + seq_along(reservoir_idx)
			}
			pred_df = predictors_df
			if (is.null(pred_df)) pred_df = private$gee_predictors_df()
			if (!is.data.frame(pred_df)) pred_df = as.data.frame(pred_df)

			y_keep = private$y[keep]
			wt_keep = if (is.null(row_weights)) NULL else as.numeric(row_weights[keep])
			dat = data.frame(y = y_keep, pred_df[keep, , drop = FALSE], group_id = group_id)
			if (!is.null(wt_keep)) {
				dat$.bootstrap_weight = wt_keep
			}
			dat = dat[order(dat$group_id), , drop = FALSE]
			id_sorted = dat$group_id
			y_sorted = dat$y
			weights_sorted = if (".bootstrap_weight" %in% colnames(dat)) as.numeric(dat$.bootstrap_weight) else NULL
			dat$group_id = NULL
			dat$y = NULL
			if (".bootstrap_weight" %in% colnames(dat)) {
				dat$.bootstrap_weight = NULL
			}
			list(dat = dat, id_sorted = id_sorted, y_sorted = y_sorted, weights_sorted = weights_sorted)
		},
		fit_gee_on_data = function(fit_data, std_err = TRUE, estimate_only = FALSE){
			if (private$use_rcpp) {
				X_rcpp = cbind(`(Intercept)` = 1, as.matrix(fit_data$dat))
				ws_args = private$gee_warm_start_args(ncol(X_rcpp))
				res = tryCatch({
					gee_pairs_singletons_cpp(
						X = X_rcpp,
						y = as.numeric(fit_data$y_sorted),
						group_id = as.integer(fit_data$id_sorted),
						family_str = private$gee_family_str(),
						warm_start_beta = ws_args$warm_start_beta,
						warm_start_fisher_info = ws_args$warm_start_fisher_info
					)
				}, error = function(e) NULL)
				if (!is.null(res) && isTRUE(res$converged)) return(res)
			}
			# Fallback to geepack
			std_err_arg = if (is.character(std_err)) std_err[1] else "san.se"
			tryCatch({
				dat_geepack = data.frame(y = fit_data$y_sorted, fit_data$dat)
				X_geepack = cbind(`(Intercept)` = 1, as.matrix(fit_data$dat))
				ws_args = private$gee_warm_start_args(ncol(X_geepack))

				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					geepack::geeglm(
						y ~ .,
						family = private$gee_family(),
						data   = dat_geepack,
						id     = fit_data$id_sorted,
						corstr = "exchangeable",
						std.err = std_err_arg,
						start  = ws_args$warm_start_beta
					)
				)))
				mod
			}, error = function(e) NULL)
		},
		fit_gee = function(std_err = TRUE, include_reservoir = TRUE, predictors_df = NULL, estimate_only = FALSE){
			fit_data = private$build_gee_fit_data(include_reservoir = include_reservoir, predictors_df = predictors_df)
			if (is.null(fit_data)) return(NULL)
			private$fit_gee_on_data(fit_data, std_err = std_err, estimate_only = estimate_only)
		},
		fit_weighted_gee_on_data = function(fit_data){
			if (is.null(fit_data) || is.null(fit_data$weights_sorted)) return(NULL)
			weights_sorted = as.numeric(fit_data$weights_sorted)
			keep = is.finite(weights_sorted) & weights_sorted > 0
			if (!any(keep)) return(NULL)
			weights_kept = weights_sorted[keep]
			id_sorted = fit_data$id_sorted[keep]
			if (length(unique(id_sorted)) < 1L) return(NULL)
			dat_kept = fit_data$dat[keep, , drop = FALSE]
			y_kept = fit_data$y_sorted[keep]
			if (private$use_rcpp) {
				X_rcpp = cbind(`(Intercept)` = 1, as.matrix(dat_kept))
				ws_args = private$gee_warm_start_args(ncol(X_rcpp))
				return(tryCatch({
					gee_pairs_singletons_weighted_cpp(
						X = X_rcpp,
						y = as.numeric(y_kept),
						group_id = as.integer(id_sorted),
						family_str = private$gee_family_str(),
						weights = as.numeric(weights_kept),
						warm_start_beta = ws_args$warm_start_beta,
						warm_start_fisher_info = ws_args$warm_start_fisher_info
					)
				}, error = function(e) NULL))
			}
			if (!requireNamespace("geepack", quietly = TRUE)) return(NULL)
			dat_geepack = data.frame(y = y_kept, dat_kept, check.names = FALSE)
			X_geepack = cbind(`(Intercept)` = 1, as.matrix(dat_geepack[, setdiff(colnames(dat_geepack), "y"), drop = FALSE]))
			ws_args = private$gee_warm_start_args(ncol(X_geepack))
			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					geepack::geeglm(
						y ~ .,
						family = private$gee_family(),
						data   = dat_geepack,
						weights = weights_kept,
						id     = id_sorted,
						corstr = "exchangeable",
						std.err = "none",
						start  = ws_args$warm_start_beta
					)
				)))
				mod
			}, error = function(e) NULL)
		},
		fit_weighted_gee_with_fallback = function(row_weights){
			gee_fit_ok = function(mod){
				if (is.null(mod)) return(FALSE)
				beta = tryCatch(
					if (inherits(mod, "geeglm")) {
						stats::coef(mod)
					} else if (!is.null(mod$beta)) {
						mod$beta
					} else {
						stats::coef(mod)
					},
					error = function(e) NULL
				)
				if (is.null(beta) || !private$gee_coefficients_are_usable(beta)) return(FALSE)
				beta_hat = private$extract_gee_treatment_estimate(mod)
				is.finite(beta_hat)
			}
			for (predictors_df in private$gee_predictors_df_candidates()) {
				fit_data = private$build_gee_fit_data(
					include_reservoir = TRUE,
					predictors_df = predictors_df,
					row_weights = row_weights
				)
				mod = private$fit_weighted_gee_on_data(fit_data)
				if (gee_fit_ok(mod)) {
					beta = if (inherits(mod, "geeglm")) stats::coef(mod) else mod$beta
					fisher = if (inherits(mod, "geeglm")) NULL else mod$fisher_information
					private$set_fit_warm_start(beta, "beta", fisher = fisher)
					return(private$extract_gee_treatment_estimate(mod))
				}
				if (private$gee_has_reservoir()) {
					fit_data_fb = private$build_gee_fit_data(
						include_reservoir = FALSE,
						predictors_df = predictors_df,
						row_weights = row_weights
					)
					mod_fb = private$fit_weighted_gee_on_data(fit_data_fb)
					if (gee_fit_ok(mod_fb)) {
						beta_fb = if (inherits(mod_fb, "geeglm")) stats::coef(mod_fb) else mod_fb$beta
						fisher_fb = if (inherits(mod_fb, "geeglm")) NULL else mod_fb$fisher_information
						private$set_fit_warm_start(beta_fb, "beta", fisher = fisher_fb)
						return(private$extract_gee_treatment_estimate(mod_fb))
					}
				}
			}
			NA_real_
		},
		fit_gee_with_fallback = function(std_err = TRUE, estimate_only = FALSE){
			gee_fit_ok = function(mod){
				if (is.null(mod)) return(FALSE)
				beta = tryCatch(
					if (inherits(mod, "geeglm")) {
						stats::coef(mod)
					} else if (!is.null(mod$beta)) {
						mod$beta
					} else {
						stats::coef(mod)
					},
					error = function(e) NULL
				)
				if (is.null(beta) || !private$gee_coefficients_are_usable(beta)) return(FALSE)
				beta_hat = private$extract_gee_treatment_estimate(mod)
				if (!is.finite(beta_hat)) return(FALSE)
				if (estimate_only) return(TRUE)
				j_treat = private$gee_treatment_index(beta)
				coef_table = if (inherits(mod, "geeglm")) tryCatch(summary(mod)$coefficients, error = function(e) NULL) else NULL
				se_hat = private$extract_gee_treatment_se(mod, j_treat = j_treat, coef_table = coef_table)
				is.finite(se_hat) && se_hat > 0 && se_hat <= private$max_abs_reasonable_coef
			}
			for (predictors_df in private$gee_predictors_df_candidates()) {
				mod = private$fit_gee(std_err = std_err, include_reservoir = TRUE, predictors_df = predictors_df, estimate_only = estimate_only)
				if (gee_fit_ok(mod)) return(mod)
				if (private$gee_has_reservoir()) {
					mod_fb = private$fit_gee(std_err = std_err, include_reservoir = FALSE, predictors_df = predictors_df, estimate_only = estimate_only)
					if (gee_fit_ok(mod_fb)) return(mod_fb)
				}
			}
			NULL
		},
		gee_family_str = function() {
			fam = private$gee_family()
			if (fam$family == "gaussian") return("gaussian")
			if (fam$family == "binomial") return("binomial")
			if (fam$family == "poisson") return("poisson")
			"gaussian"
		}
	)
)
