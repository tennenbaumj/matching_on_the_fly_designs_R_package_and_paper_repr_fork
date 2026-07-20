#' Marginal (Cluster-Robust) Weibull Inference for KK Matched-Pair Survival Designs
#'
#' Fits a single pooled Weibull Accelerated Failure Time (AFT) model across all
#' subjects (treatment plus, optionally, all recorded covariates), ignoring the
#' matched-pair structure in the mean model. Standard errors are computed via a
#' cluster-robust (sandwich) covariance estimator: matched pairs from a KK
#' matching-on-the-fly or binary-match design form size-2 clusters, and
#' unmatched reservoir subjects each form their own singleton cluster.
#'
#' @details
#' This is the "marginal" competitor to \code{\link{InferenceSurvivalKKWeibullFrailtyOneLik}}:
#' rather than modeling the within-pair correlation explicitly via a frailty
#' term, it fits an ordinary (working-independence) Weibull AFT model and
#' corrects the treatment-effect standard error post hoc for the within-pair
#' dependence. The model is fit via the package's fast C++ Weibull AFT backend
#' (\code{fast_weibull_regression_cpp}) and the cluster-robust sandwich is
#' assembled from per-subject dfbeta contributions collapsed within clusters,
#' which is numerically equivalent to
#' \code{survival::survreg(..., cluster = ..., robust = TRUE)} (retained as a
#' fallback if the C++ fit fails to converge).
#'
#' @examples
#' \donttest{
#' des = DesignSeqOneByOneKK14$new(n = 20, response_type = 'survival')
#' for (i in 1:20) {
#'   des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' des$add_all_subject_responses(rexp(20))
#' inf = InferenceSurvivalKKWeibullMarginal$new(des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKWeibullMarginal = R6::R6Class("InferenceSurvivalKKWeibullMarginal",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the marginal (cluster-robust) Weibull inference object.
		#' @param des_obj A completed KK or BinaryMatch design with survival response.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj = des_obj, verbose = verbose, model_formula = model_formula)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the pooled treatment effect estimate (log-time-ratio scale).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the treatment estimate under Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, compute only the weighted point estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) return(beta_hat_T)
			}
			X_cov = private$get_X()
			X_fit = if (!is.null(X_cov) && ncol(as.matrix(X_cov)) > 0L) {
				cbind(treatment = private$w, X_cov)
			} else {
				matrix(private$w, ncol = 1L, dimnames = list(NULL, "treatment"))
			}
			fit = weighted_weibull_bootstrap_surrogate_fit(
				private$y, private$dead, X_fit, row_weights, cluster = private$get_cluster_ids(),
				warm_start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + 2L)
			)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes the asymptotic (cluster-robust) confidence interval.
		#' @param alpha The significance level. Default 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the asymptotic (cluster-robust) two-sided p-value.
		#' @param delta Null treatment effect value. Default 0.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the treatment effect estimate.
		#' @param B Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		},
		#' @description Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		cached_mod = NULL,
		best_X_colnames = NULL,
		cached_vc_params = NULL,
		max_abs_reasonable_coef = 1e4,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		# Matched-pair members share a cluster id; reservoir (unmatched) subjects each
		# get a unique singleton id. Mirrors InferenceAbstractKKMarginalIncid's private
		# get_cluster_ids() (inference_incidence_KK_marginal_abstract.R), duplicated here
		# since that class is incidence-only (asserts no censoring).
		get_cluster_ids = function(){
			des_priv = private$des_obj_priv_int
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec_int = as.integer(m_vec)
			m_vec_int[is.na(m_vec_int)] = 0L
			des_m = des_priv$m
			if (is.null(des_m)) des_m = rep(NA_integer_, private$n)
			des_m_int = as.integer(des_m)
			des_m_int[is.na(des_m_int)] = 0L
			if (!is.null(des_priv$cluster_id) && identical(m_vec_int, des_m_int)){
				return(des_priv$cluster_id)
			}
			if (!is.null(private$cached_values$cluster_id) &&
				identical(m_vec_int, private$cached_values$cluster_id_m_vec)){
				return(private$cached_values$cluster_id)
			}
			cluster_id = des_priv$compute_matching_cluster_ids(m_vec_int)
			if (identical(m_vec_int, des_m_int)){
				des_priv$cluster_id = cluster_id
				des_priv$cluster_id_m_vec = m_vec_int
			} else {
				private$cached_values$cluster_id = cluster_id
				private$cached_values$cluster_id_m_vec = m_vec_int
			}
			cluster_id
		},
		# Fits the pooled Weibull AFT model via fast_weibull_regression_cpp (X_fit must
		# carry an explicit "(Intercept)" column plus "treatment"). When robust=TRUE the
		# cluster-robust sandwich is built from per-subject dfbeta rows (score %*% vcov)
		# summed within clusters and crossprod'ed -- the same estimator survreg computes
		# from resid(fit, "dfbeta") under cluster= / robust=TRUE. Returns NULL on any
		# failure so callers can fall back to the survreg path.
		fit_weibull_marginal_cpp = function(X_fit, robust = TRUE, cluster_ids = NULL){
			y = private$y
			dead = private$dead
			ok = is.finite(y) & is.finite(dead)
			if (!any(ok)) return(NULL)
			X_ok = as.matrix(X_fit[ok, , drop = FALSE])
			j_treat = match("treatment", colnames(X_ok))
			if (is.na(j_treat)) return(NULL)
			p = ncol(X_ok)
			n_params = p + 1L
			res = tryCatch(
				fast_weibull_regression_cpp(
					y = as.numeric(y[ok]), dead = as.numeric(dead[ok]), X = X_ok,
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					estimate_only = !robust
				),
				error = function(e) NULL
			)
			if (is.null(res) || !isTRUE(res$converged)) return(NULL)
			if (!robust) {
				b = as.numeric(res$b)
				beta_T = b[j_treat]
				if (!is.finite(beta_T)) return(NULL)
				private$set_fit_warm_start(c(b, as.numeric(res$log_sigma)), "params")
				return(list(beta_T = beta_T, se_T = NA_real_, fit_obj = res))
			}
			params = as.numeric(res$params)
			vc = res$vcov
			if (length(params) != n_params || !all(is.finite(params)) ||
				is.null(vc) || !all(is.finite(vc))) return(NULL)
			beta_T = params[j_treat]
			se_T = NA_real_
			if (!is.null(cluster_ids)) {
				# Per-obs score of the Weibull AFT loglik at the MLE, params = [beta, log_sigma]:
				# w_i = (log y_i - x_i'beta)/sigma, s_beta_i = x_i (e^{w_i} - dead_i)/sigma,
				# s_logsigma_i = w_i (e^{w_i} - dead_i) - dead_i.
				sigma = exp(params[n_params])
				d_vec = as.numeric(dead[ok])
				w_vec = (log(as.numeric(y[ok])) - as.numeric(X_ok %*% params[seq_len(p)])) / sigma
				resid_vec = exp(w_vec) - d_vec
				U = cbind(X_ok * (resid_vec / sigma), w_vec * resid_vec - d_vec)
				if (all(is.finite(U))) {
					D_c = rowsum(U %*% vc, group = as.integer(cluster_ids[ok]))
					v_T = crossprod(D_c)[j_treat, j_treat]
					if (is.finite(v_T) && v_T > 0) se_T = sqrt(v_T)
				}
			}
			private$set_fit_warm_start(params, "params", fisher = res$information)
			list(beta_T = beta_T, se_T = se_T, fit_obj = res)
		},
		# Fits a pooled Weibull AFT model on X_fit (first column must be "treatment").
		# When robust=TRUE, computes a cluster-robust sandwich SE clustering on cluster_ids;
		# when FALSE, skips the robust/cluster machinery entirely for a cheap point-estimate-only fit.
		# Fallback path when the C++ backend fails to converge.
		fit_weibull_marginal_survreg = function(X_fit, robust = TRUE, cluster_ids = NULL){
			y = private$y
			dead = private$dead
			ok = is.finite(y) & is.finite(dead)
			if (!any(ok)) return(NULL)
			X_ok = X_fit[ok, , drop = FALSE]
			dat = as.data.frame(X_ok, check.names = FALSE)
			dat$.time__ = as.numeric(y[ok])
			dat$.dead__ = as.numeric(dead[ok])
			rhs = paste(colnames(X_ok), collapse = " + ")
			fmla = stats::as.formula(paste0("survival::Surv(.time__, .dead__) ~ ", rhs))
			fit = if (isTRUE(robust) && !is.null(cluster_ids)) {
				dat$.cluster__ = cluster_ids[ok]
				tryCatch(
					suppressWarnings(survival::survreg(fmla, data = dat, dist = "weibull", cluster = .cluster__, robust = TRUE)),
					error = function(e) NULL
				)
			} else {
				tryCatch(
					suppressWarnings(survival::survreg(fmla, data = dat, dist = "weibull")),
					error = function(e) NULL
				)
			}
			if (is.null(fit)) return(NULL)
			coefs = tryCatch(stats::coef(fit), error = function(e) NULL)
			if (is.null(coefs) || !("treatment" %in% names(coefs))) return(NULL)
			beta_T = as.numeric(coefs[["treatment"]])
			se_T = NA_real_
			if (isTRUE(robust) && !is.null(cluster_ids)) {
				vc = tryCatch(stats::vcov(fit), error = function(e) NULL)
				if (!is.null(vc) && "treatment" %in% rownames(vc)) {
					v = vc["treatment", "treatment"]
					if (is.finite(v) && v > 0) se_T = sqrt(v)
				}
			}
			list(beta_T = beta_T, se_T = se_T, fit_obj = fit)
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

			if (sum(private$dead) == 0L){
				private$cache_nonestimable_estimate("kk_weibull_marginal_no_events")
				return(invisible(NULL))
			}

			X_cov = private$get_X()
			X_full = if (is.null(X_cov) || ncol(as.matrix(X_cov)) == 0L) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, as.matrix(X_cov))
			}
			cluster_ids = if (!estimate_only) private$get_cluster_ids() else NULL

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit){
					fit = private$fit_weibull_marginal_cpp(X_fit, robust = !estimate_only, cluster_ids = cluster_ids)
					if (!is.null(fit)) return(fit)
					# survreg builds its own intercept from the formula, so drop the explicit column
					private$fit_weibull_marginal_survreg(X_fit[, -1L, drop = FALSE], robust = !estimate_only, cluster_ids = cluster_ids)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || !is.finite(mod$beta_T) || abs(mod$beta_T) > private$max_abs_reasonable_coef) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$se_T) && mod$se_T > 0
				}
			)

			fit = attempt$fit
			if (is.null(fit)) {
				private$cache_nonestimable_estimate("kk_weibull_marginal_fit_failed")
				return(invisible(NULL))
			}

			private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
			private$cached_mod = fit$fit_obj
			log_s = tryCatch(as.numeric(fit$fit_obj$log_sigma), error = function(e) NULL)
			if (isTRUE(is.finite(log_s))) private$cached_vc_params = log_s
			private$cached_values$beta_hat_T = fit$beta_T

			if (!estimate_only) {
				private$cached_values$s_beta_hat_T = fit$se_T
				private$cached_values$df = Inf
			}
			invisible(NULL)
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			se = private$cached_values$s_beta_hat_T
			if (is.null(se) || length(se) == 0L) return(NA_real_)
			se
		},
		get_degrees_of_freedom = function() Inf,
		# Cheap point-estimate-only refit under a re-randomized private$w, reusing the
		# covariate columns selected by the original shared() fit. No cluster/robust SE
		# is computed here since randomization inference only needs the point estimate.
		# Runs entirely through the C++ backend with warm starts chained across
		# permutations (set_fit_warm_start allows this for the "rand" operation).
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_X_colnames)) {
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)) {
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_data = private$get_X()
			X_cov = if (length(private$best_X_colnames) == 0L) {
				NULL
			} else {
				X_data[, intersect(private$best_X_colnames, colnames(X_data)), drop = FALSE]
			}
			X_fit = if (is.null(X_cov) || ncol(X_cov) == 0L) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			# Fixed-VC fast path
			if (!is.null(private$cached_vc_params) && is.finite(private$cached_vc_params[1L])) {
				ok = is.finite(private$y) & is.finite(private$dead)
				X_ok = as.matrix(X_fit[ok, , drop = FALSE])
				n_params = ncol(X_ok) + 1L
				res_fast = tryCatch(
					fast_weibull_regression_cpp(
						y    = as.numeric(private$y[ok]),
						dead = as.numeric(private$dead[ok]),
						X    = X_ok,
						warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
						estimate_only = TRUE,
						fixed_idx    = as.integer(n_params),
						fixed_values = private$cached_vc_params[1L]
					),
					error = function(e) NULL
				)
				if (!is.null(res_fast) && isTRUE(res_fast$converged)) {
					b = as.numeric(res_fast$b)
					j_treat = match("treatment", colnames(X_fit))
					if (!is.na(j_treat) && length(b) >= j_treat && is.finite(b[j_treat])) {
						private$set_fit_warm_start(c(b, private$cached_vc_params[1L]), "params")
						return(as.numeric(b[j_treat]))
					}
				}
			}
			fit = private$fit_weibull_marginal_cpp(X_fit, robust = FALSE, cluster_ids = NULL)
			if (is.null(fit)) {
				fit = private$fit_weibull_marginal_survreg(X_fit[, -1L, drop = FALSE], robust = FALSE, cluster_ids = NULL)
			}
			if (is.null(fit) || !is.finite(fit$beta_T)) return(NA_real_)
			fit$beta_T
		},
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_data = private$get_X()
			Xc = if (!is.null(private$best_X_colnames)) {
				cols = intersect(private$best_X_colnames, colnames(X_data))
				if (length(cols) > 0L) as.matrix(X_data[, cols, drop = FALSE]) else matrix(numeric(0), nrow = private$n, ncol = 0L)
			} else if (!is.null(X_data) && ncol(as.matrix(X_data)) > 0L) {
				as.matrix(X_data)
			} else {
				matrix(numeric(0), nrow = private$n, ncol = 0L)
			}
			compute_weibull_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		}
	)))
)
