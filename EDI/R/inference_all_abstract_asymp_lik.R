#' Likelihood-Backed Asymptotic Inference
#'
#' @name InferenceAsympLik
#' @description Intermediate base class for asymptotic inference families that expose
#' likelihood / partial-likelihood / working-likelihood test paths in addition
#' to Wald inference. The term "likelihood" is used broadly: subclasses may be
#' backed by a true full likelihood, a partial likelihood (e.g. Cox PH), a
#' quasi-likelihood (e.g. GEE, quasi-Poisson), or a composite/combined
#' likelihood. Classes requiring a full generative likelihood — i.e. those
#' supporting parametric-bootstrap LR calibration — inherit instead from
#' \code{InferenceParamBootstrap}.
#'
#' @keywords internal
InferenceAsympLik = R6::R6Class("InferenceAsympLik",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMSummaryTable,
	public = list(
		#' @description Computes an asymptotic confidence interval using the configured test.
		#'
		#' @param alpha Significance level 1 - \code{alpha}. Default 0.05.
		#'
		#' @return A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			switch(
				private$testing_type,
				wald = private$compute_wald_confidence_interval_impl(alpha),
				score = private$compute_score_confidence_interval_impl(alpha),
				gradient = private$compute_gradient_confidence_interval_impl(alpha),
				lik_ratio = private$compute_lik_ratio_confidence_interval_impl(alpha)
			)
		},
		#' @description Computes an asymptotic two-sided p-value using the configured test.
		#'
		#' @param delta Null treatment effect to test against. Default 0.
		#'
		#' @return The asymptotic p-value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			switch(
				private$testing_type,
				wald = private$compute_wald_two_sided_pval_impl(delta),
				score = private$compute_score_two_sided_pval_impl(delta),
				gradient = private$compute_gradient_two_sided_pval_impl(delta),
				lik_ratio = private$compute_lik_ratio_two_sided_pval_impl(delta)
			)
		},
		#' @description Sets the asymptotic testing method used by p-values and CIs.
		#'
		#' @param testing_type One of \code{"wald"}, \code{"score"}, \code{"gradient"}, or \code{"lik_ratio"}.
		#'
		#' @return The inference object, invisibly.
		set_testing_type = function(testing_type = c("wald", "score", "gradient", "lik_ratio")){
			testing_type = private$normalize_testing_type(testing_type)
			supported = private$get_supported_testing_types_impl()
			if (!testing_type %in% supported) {
				stop(
					class(self)[1], " does not support testing_type = \"", testing_type,
					"\". Supported values are: ", paste(supported, collapse = ", "),
					call. = FALSE
				)
			}
			private$testing_type = testing_type
			invisible(self)
		},
		#' @description Sets the information matrix preference used by score-test dispatch.
		#'
		#' @param information_preference One of \code{"auto"}, \code{"fisher"}, or \code{"observed"}.
		#'
		#' @return The inference object, invisibly.
		set_information_preference = function(information_preference = c("auto", "fisher", "observed")){
			information_preference = private$normalize_information_preference(information_preference)
			supported = private$get_supported_information_preferences_impl()
			if (!information_preference %in% supported) {
				stop(
					class(self)[1], " does not support information_preference = \"", information_preference,
					"\". Supported values are: ", paste(supported, collapse = ", "),
					call. = FALSE
				)
			}
			private$information_preference = information_preference
			private$information_source_used = NULL
			private$clear_likelihood_test_eval_cache()
			invisible(self)
		},
		#' @description Gets the asymptotic testing method used by p-values and CIs.
		get_testing_type = function(){
			private$testing_type
		},
		#' @description Gets the score-test information matrix preference.
		get_information_preference = function(){
			private$information_preference
		},
		#' @description Gets the actual information source used by the most recent information-backed computation.
		get_information_source_used = function(){
			private$information_source_used
		},
		#' @description Gets the asymptotic testing methods supported by this inference object.
		get_supported_testing_types = function(){
			private$get_supported_testing_types_impl()
		},
		#' @description Gets the score-test information matrix preferences supported by this inference object.
		get_supported_information_preferences = function(){
			private$get_supported_information_preferences_impl()
		},
		#' @description Computes the score two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_score_two_sided_pval = function(delta = 0){
			private$compute_score_two_sided_pval_impl(delta)
		},
		#' @description Computes the score confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_score_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_score_confidence_interval_impl(alpha)
		},

		#' @description Computes the likelihood-ratio two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			private$compute_lik_ratio_two_sided_pval_impl(delta)
		},

		#' @description Computes the likelihood-ratio confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_lik_ratio_confidence_interval_impl(alpha)
		},

		#' @description Computes the gradient two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_gradient_two_sided_pval = function(delta = 0){
			private$compute_gradient_two_sided_pval_impl(delta)
		},

		#' @description Computes the gradient confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_gradient_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_gradient_confidence_interval_impl(alpha)
		}
	),
	private = list(
		likelihood_ci_max_abs = 10,
		testing_type = "wald",
		information_preference = "auto",
		information_source_used = NULL,

		get_information_matrix = function(spec = NULL, fit = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec)) return(NULL)
			if (is.null(fit)) {
				fit = spec$full_fit %||% private$cached_mod
			}
			if (is.null(fit)) return(NULL)

			extract_fisher = function(){
				tryCatch({
					if (!is.null(spec$fisher_information)) {
						spec$fisher_information(fit)
					} else if (!is.null(fit$fisher_information)) {
						fit$fisher_information
					} else if (identical(fit$information_type %||% "", "fisher") && !is.null(fit$information)) {
						fit$information
					} else {
						NULL
					}
				}, error = function(e) NULL)
			}

			extract_observed = function(){
				tryCatch({
					if (!is.null(spec$observed_information)) {
						spec$observed_information(fit)
					} else if (!is.null(fit$observed_information)) {
						fit$observed_information
					} else if (identical(fit$information_type %||% "", "observed") && !is.null(fit$information)) {
						fit$information
					} else {
						NULL
					}
				}, error = function(e) NULL)
			}

			extract_legacy = function(){
				tryCatch({
					if (!is.null(spec$information)) {
						spec$information(fit)
					} else {
						fit$information
					}
				}, error = function(e) NULL)
			}

			preference = private$information_preference
			if (identical(preference, "auto")) {
				preference = private$get_default_information_source()
			}
			if (identical(preference, "fisher")) {
				information = extract_fisher()
				if (is.null(information)) {
					stop(class(self)[1], " does not expose Fisher information for information-backed inference.", call. = FALSE)
				}
				private$information_source_used = "fisher"
				return(information)
			}

			if (identical(preference, "observed")) {
				information = extract_observed()
				if (is.null(information)) {
					stop(class(self)[1], " does not expose observed information for information-backed inference.", call. = FALSE)
				}
				private$information_source_used = "observed"
				return(information)
			}

			information = extract_fisher()
			if (!is.null(information)) {
				private$information_source_used = "fisher"
				return(information)
			}
			information = extract_observed()
			if (!is.null(information)) {
				private$information_source_used = "observed"
				return(information)
			}
			information = extract_legacy()
			if (!is.null(information)) {
				private$information_source_used = "legacy"
			}
			information
		},

		compute_variance_from_information_matrix = function(information, j){
			information = as.matrix(information)
			if (!is.matrix(information) || nrow(information) != ncol(information) || length(j) != 1L ||
				!is.finite(j) || j < 1L || j > nrow(information)) {
				return(NA_real_)
			}
			if (nrow(information) == 1L) {
				val = as.numeric(information[1L, 1L])
				return(if (is.finite(val) && val > 0) 1 / val else NA_real_)
			}
			# Use specialized C++ helper for single diagonal entry to avoid full matrix inversion in R
			res = tryCatch(eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(information, as.integer(j)), error = function(e) NA_real_)
			if (is.finite(res) && res > 0) return(res)
			
			# Fallback to full inversion if C++ helper fails or returns invalid result
			vcov = tryCatch(solve(information), error = function(e) NULL)
			if (is.null(vcov) || any(!is.finite(vcov))) {
				vcov = tryCatch(qr.solve(information, diag(nrow(information))), error = function(e) NULL)
			}
			if (is.null(vcov) || any(!is.finite(vcov))) return(NA_real_)
			vcov = (vcov + t(vcov)) / 2
			as.numeric(vcov[j, j])
		},

		compute_standard_error_from_information_matrix = function(spec = NULL, fit = NULL, j = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec)) return(NA_real_)
			if (is.null(j)) {
				j = as.integer(spec$j)
			}
			information = tryCatch(private$get_information_matrix(spec = spec, fit = fit), error = function(e) NULL)
			if (is.null(information)) return(NA_real_)
			variance = private$compute_variance_from_information_matrix(information, j)
			if (is.finite(variance) && variance >= 0) sqrt(variance) else NA_real_
		},

		compute_score_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_score_confidence_interval_impl = function(alpha){
			private$invert_test_pval_confidence_interval(alpha, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_gradient_confidence_interval_impl = function(alpha){
			private$invert_gradient_ci_uniroot(alpha)
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		compute_lik_ratio_confidence_interval_impl = function(alpha){
			private$invert_lik_ratio_ci_newton(alpha)
		},


		get_likelihood_test_spec = function(){
			NULL
		},

		get_score_test_information_matrix = function(spec, fit){
			private$get_information_matrix(spec = spec, fit = fit)
		},

		make_warm_fit_null_wrapper = function(spec, cache_key){
			last_start = NULL
			last_delta = NULL
			fit_null_formals = tryCatch(names(formals(spec$fit_null)), error = function(e) character())
			accepts_start = "start" %in% fit_null_formals
			function(delta){
				warm_enabled = isTRUE(private$null_fit_warm_start_enabled)
				cache_state = if (warm_enabled) private$get_likelihood_null_warm_state(cache_key) else NULL
				start = if (warm_enabled) last_start else NULL
				if (warm_enabled && is.null(start) && !is.null(cache_state)) start = cache_state$start
				if (!is.null(start) && length(start) == 0L) start = NULL
				fit = tryCatch(
					if (accepts_start) spec$fit_null(delta, start = start) else spec$fit_null(delta),
					error = function(e) NULL
				)
				extract_start = spec$extract_start %||% function(fit_obj) {
					fit_obj$params %||% fit_obj$b %||% {
						co = fit_obj$coefficients
						if (is.numeric(co)) as.numeric(co) else NULL
					}
				}
				last_start_val = if (warm_enabled && accepts_start && !is.null(fit)) {
					tryCatch(extract_start(fit), error = function(e) NULL)
				} else NULL
				if (!is.null(last_start_val) && length(last_start_val) == 0L) last_start_val = NULL
				last_start <<- last_start_val
				last_delta <<- delta
				if (warm_enabled && accepts_start) {
					private$set_likelihood_null_warm_state(cache_key, delta = delta, start = last_start)
				}
				fit
			}
		},

		get_memoized_likelihood_test_eval = function(delta, testing_type, spec = NULL, warm_cache_key = NULL, include_null_fit = TRUE, include_score = FALSE, include_information = FALSE, include_full_negloglik = FALSE, include_null_negloglik = FALSE){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec)) {
				stop(class(self)[1], " does not expose a likelihood-test specification.", call. = FALSE)
			}

			if (is.null(warm_cache_key)) warm_cache_key = paste0("likelihood_test:", testing_type)
			delta = private$normalize_likelihood_test_delta(delta)
			entry = private$get_likelihood_test_eval_entry(testing_type, delta)
			if (is.null(entry)) {
				entry = list(
					delta = delta,
					testing_type = testing_type
				)
			}

			j = as.integer(spec$j)
			if (length(j) != 1L || !is.finite(j) || j < 1L) {
				entry$invalid = TRUE
				private$set_likelihood_test_eval_entry(testing_type, delta, entry)
				return(entry)
			}
			entry$j = j

			if (isTRUE(include_null_fit) && is.null(entry$null_fit)) {
				fit_null = private$make_warm_fit_null_wrapper(spec, cache_key = warm_cache_key)
				null_fit = tryCatch(fit_null(delta), error = function(e) NULL)
				null_params = if (!is.null(null_fit)) {
					null_fit$params %||% null_fit$b %||% {
						co = null_fit$coefficients
						if (is.numeric(co)) as.numeric(co) else NULL
					}
				} else NULL
				if (is.null(null_fit) || is.null(null_params) || length(null_params) < j || !is.finite(null_params[j])) {
					entry$invalid = TRUE
					private$set_likelihood_test_eval_entry(testing_type, delta, entry)
					return(entry)
				}
				attr(null_fit, "edi_likelihood_test_delta") = delta
				attr(null_fit, "edi_likelihood_test_testing_type") = testing_type
				entry$null_fit = null_fit
				entry$null_params = as.numeric(null_params)
				entry$invalid = FALSE
			}
			if (isTRUE(include_null_fit) && !is.null(warm_cache_key) && !is.null(entry$null_params)) {
				cache_state = private$get_likelihood_null_warm_state(warm_cache_key)
				if (is.null(cache_state) || is.null(cache_state$start)) {
					private$set_likelihood_null_warm_state(warm_cache_key, delta = delta, start = as.numeric(entry$null_params))
				}
			}

			if (isTRUE(include_score) && !is.null(entry$null_fit) && is.null(entry$score)) {
				entry$score = tryCatch(spec$score(entry$null_fit), error = function(e) NULL)
			}

			if (isTRUE(include_information) && !is.null(entry$null_fit) && is.null(entry$information)) {
				entry$information = tryCatch(private$get_score_test_information_matrix(spec, entry$null_fit), error = function(e) NULL)
			}

			if (isTRUE(include_full_negloglik) && is.null(entry$full_negloglik)) {
				entry$full_negloglik = tryCatch(spec$neg_loglik(spec$full_fit), error = function(e) NA_real_)
			}

			if (isTRUE(include_null_negloglik) && !is.null(entry$null_fit) && is.null(entry$null_negloglik)) {
				entry$null_negloglik = tryCatch(spec$neg_loglik(entry$null_fit), error = function(e) NA_real_)
			}

			private$set_likelihood_test_eval_entry(testing_type, delta, entry)
			entry
		},

		get_memoized_likelihood_test_pval = function(delta, testing_type, spec = NULL, warm_cache_key = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec)) {
				stop(class(self)[1], " does not expose a likelihood-test specification.", call. = FALSE)
			}

			need_score = testing_type %in% c("score", "gradient")
			need_information = identical(testing_type, "score")
			need_full_negloglik = identical(testing_type, "lik_ratio")
			need_null_negloglik = identical(testing_type, "lik_ratio")
			entry = private$get_memoized_likelihood_test_eval(
				delta = delta,
				testing_type = testing_type,
				spec = spec,
				warm_cache_key = warm_cache_key,
				include_score = need_score,
				include_information = need_information,
				include_full_negloglik = need_full_negloglik,
				include_null_negloglik = need_null_negloglik
			)
			if (isTRUE(entry$invalid) || !is.finite(entry$j)) return(NA_real_)
			if (!is.null(entry$p_value) && is.finite(entry$p_value)) return(entry$p_value)

			j = entry$j
			p_value = NA_real_
			if (testing_type == "score") {
				if (is.null(entry$score) || is.null(entry$information)) return(NA_real_)
				res = score_test_from_score_information_cpp(as.numeric(entry$score), as.matrix(entry$information), j)
				p_value = as.numeric(res$p_value %||% res)
			} else if (testing_type == "gradient") {
				est = self$compute_estimate()
				if (is.null(entry$score) || !is.finite(est)) return(NA_real_)
				res = gradient_test_from_restricted_score_cpp(as.numeric(entry$score), est, delta, j)
				p_value = as.numeric(res$p_value %||% res)
			} else if (testing_type == "lik_ratio") {
				f_nll = entry$full_negloglik
				n_nll = entry$null_negloglik
				if (is.null(f_nll) || !length(f_nll) || is.null(n_nll) || !length(n_nll)) return(NA_real_)
				if (!is.finite(f_nll) || !is.finite(n_nll)) return(NA_real_)
				res = likelihood_ratio_test_from_negloglik_cpp(f_nll, n_nll, df = 1L)
				p_value = as.numeric(res$p_value %||% res)
			} else {
				stop("Unsupported testing_type: ", testing_type, call. = FALSE)
			}

			entry$p_value = p_value
			private$set_likelihood_test_eval_entry(testing_type, delta, entry)
			p_value
		},

		compute_likelihood_test_two_sided_pval = function(delta, testing_type){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) {
				if (!isTRUE(self$is_nonestimable())) {
					private$cache_nonestimable_estimate("likelihood_test_spec_unavailable")
				}
				return(NA_real_)
			}
			p_value = private$get_memoized_likelihood_test_pval(
				delta = delta,
				testing_type = testing_type,
				spec = spec,
				warm_cache_key = paste0("likelihood_test:", testing_type)
			)
			if (!is.finite(p_value) && !isTRUE(self$is_nonestimable("estimate"))) {
				private$cache_nonestimable_se(paste0(testing_type, "_test_unavailable"))
			}
			p_value
		},

		invert_test_pval_confidence_interval = function(alpha, testing_type = private$testing_type){
			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			se = private$get_standard_error()
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = private$compute_wald_confidence_interval_impl(alpha)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else NA_real_
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else NA_real_

			testing_type = private$normalize_testing_type(testing_type)
			spec = private$get_likelihood_test_spec()
			pval_fn = if (!is.null(spec) && testing_type %in% c("score", "gradient", "lik_ratio")) {
				function(delta) {
					private$get_memoized_likelihood_test_pval(
						delta = delta,
						testing_type = testing_type,
						spec = spec,
						warm_cache_key = if (identical(testing_type, "score")) {
							"likelihood_test:score"
						} else {
							paste0(testing_type, "_ci")
						}
					)
				}
			} else {
				function(delta) self$compute_asymp_two_sided_pval(delta)
			}

			# A likelihood-backed CI cannot be inverted when the test is unavailable
			# at the fitted estimate.  In particular, some Cox score paths expose
			# an explicit non-estimable standard error and otherwise pass NA values
			# into the C++ root search, which may return an empty endpoint vector.
			p_est = tryCatch(pval_fn(est), error = function(e) NA_real_)
			if (!is.finite(p_est)) {
				private$cache_nonestimable_se(paste0(testing_type, "_test_unavailable"))
				return(c(NA_real_, NA_real_))
			}

			ci_vals = pval_invert_ci_cpp(
				pval_fn    = pval_fn,
				est        = est,
				alpha      = alpha,
				step       = step,
				lower_seed = lower_seed,
				upper_seed = upper_seed
			)

			ci = if (length(ci_vals) == 2L) sort(as.numeric(ci_vals[1:2]), na.last = TRUE) else c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			if (length(ci) < 2L || !all(is.finite(ci[1:2])) || ci[1L] > est || ci[2L] < est || any(abs(ci[1:2]) > private$likelihood_ci_max_abs)) {
				fallback = sort(as.numeric(wald_ci[1:2]), na.last = TRUE)
				if (length(fallback) >= 2L && all(is.finite(fallback)) && fallback[1L] <= est && fallback[2L] >= est && all(abs(fallback[1:2]) <= private$likelihood_ci_max_abs)) {
					ci = fallback
				} else {
					private$cache_nonestimable_se("test_inversion_confidence_interval_unavailable")
					ci = c(NA_real_, NA_real_)
				}
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			}
			ci
		},

		invert_gradient_ci_uniroot = function(alpha){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) return(private$invert_test_pval_confidence_interval(alpha))

			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			j = as.integer(spec$j)
			if (length(j) != 1L || !is.finite(j) || j < 1L) return(c(NA_real_, NA_real_))

			se = private$get_standard_error()
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = private$compute_wald_confidence_interval_impl(alpha)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else est - step
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else est + step

			pval_fn = function(delta){
				private$get_memoized_likelihood_test_pval(
					delta = delta,
					testing_type = "gradient",
					spec = spec,
					warm_cache_key = "gradient_ci"
				)
			}

			p_est = tryCatch(pval_fn(est), error = function(e) NA_real_)
			if (!is.finite(p_est)) {
				private$cache_nonestimable_se("gradient_test_unavailable")
				return(c(NA_real_, NA_real_))
			}
			if (p_est < alpha) {
				ci = c(est, est)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}

			ci_vals = tryCatch(
				pval_invert_ci_cpp(
					pval_fn    = pval_fn,
					est        = est,
					alpha      = alpha,
					step       = step,
					lower_seed = lower_seed,
					upper_seed = upper_seed
				),
				error = function(e) {
					private$cache_nonestimable_se("gradient_ci_inversion_failed")
					numeric(0)
				}
			)
			ci = if (length(ci_vals) == 2L) sort(as.numeric(ci_vals[1:2]), na.last = TRUE) else c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			if (length(ci) < 2L || !all(is.finite(ci[1:2])) || ci[1L] > est || ci[2L] < est || any(abs(ci[1:2]) > private$likelihood_ci_max_abs)) {
				fallback = sort(as.numeric(wald_ci[1:2]), na.last = TRUE)
				if (length(fallback) >= 2L && all(is.finite(fallback)) && fallback[1L] <= est && fallback[2L] >= est && all(abs(fallback[1:2]) <= private$likelihood_ci_max_abs)) {
					ci = fallback
				} else {
					private$cache_nonestimable_se("gradient_confidence_interval_unavailable")
					ci = c(NA_real_, NA_real_)
				}
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			}
			ci
		},

		invert_lik_ratio_ci_newton = function(alpha){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) return(private$invert_test_pval_confidence_interval(alpha))

			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			full_eval = tryCatch(private$get_memoized_likelihood_test_eval(
				delta = est,
				testing_type = "lik_ratio",
				spec = spec,
				warm_cache_key = "lik_ratio_ci",
				include_null_fit = FALSE,
				include_full_negloglik = TRUE
			), error = function(e) list(full_negloglik = NA_real_))
			full_negloglik = full_eval$full_negloglik %||% NA_real_
			if (!is.finite(full_negloglik)) return(private$invert_test_pval_confidence_interval(alpha))

			se = private$get_standard_error()
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = private$compute_wald_confidence_interval_impl(alpha)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else est - step
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else est + step

			j = as.integer(spec$j)
				fit_null_fn = function(delta){
					eval = tryCatch(
						private$get_memoized_likelihood_test_eval(
							delta = delta,
							testing_type = "lik_ratio",
							spec = spec,
							warm_cache_key = "lik_ratio_ci",
							include_score = FALSE,
							include_full_negloglik = FALSE,
							include_null_negloglik = FALSE
						),
						error = function(e) NULL
					)
					if (is.null(eval) || isTRUE(eval$invalid) || is.null(eval$null_fit)) return(NULL)
					eval$null_fit
				}
				neg_loglik_fn = function(fit){
					delta_fit = attr(fit, "edi_likelihood_test_delta", exact = TRUE)
					eval = if (is.finite(delta_fit %||% NA_real_)) {
						tryCatch(
							private$get_memoized_likelihood_test_eval(
								delta = delta_fit,
								testing_type = "lik_ratio",
								spec = spec,
								warm_cache_key = "lik_ratio_ci",
								include_score = FALSE,
								include_full_negloglik = FALSE,
								include_null_negloglik = TRUE
							),
							error = function(e) NULL
						)
					} else NULL
					if (!is.null(eval) && is.finite(eval$null_negloglik %||% NA_real_)) return(eval$null_negloglik)
					tryCatch(spec$neg_loglik(fit), error = function(e) NA_real_)
				}
				score_fn = function(fit){
					delta_fit = attr(fit, "edi_likelihood_test_delta", exact = TRUE)
					eval = if (is.finite(delta_fit %||% NA_real_)) {
						tryCatch(
							private$get_memoized_likelihood_test_eval(
								delta = delta_fit,
								testing_type = "lik_ratio",
								spec = spec,
								warm_cache_key = "lik_ratio_ci",
								include_score = TRUE,
								include_full_negloglik = FALSE,
								include_null_negloglik = FALSE
							),
							error = function(e) NULL
						)
					} else NULL
					if (!is.null(eval) && !is.null(eval$score)) return(eval$score)
					tryCatch(spec$score(fit), error = function(e) NULL)
				}

			ci_vals = tryCatch(lrt_ci_nr_cpp(
				fit_null_fn    = fit_null_fn,
				neg_loglik_fn  = neg_loglik_fn,
				score_fn       = score_fn,
				est            = est,
				full_negloglik = full_negloglik,
				alpha          = alpha,
				step           = step,
				lower_seed     = lower_seed,
				upper_seed     = upper_seed,
				j              = j
			), error = function(e) numeric(0))

			ci = if (length(ci_vals) == 2L) sort(as.numeric(ci_vals[1:2]), na.last = TRUE) else c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			if (length(ci) < 2L || !all(is.finite(ci[1:2])) || ci[1L] > est || ci[2L] < est || any(abs(ci[1:2]) > private$likelihood_ci_max_abs)) {
				fallback = sort(as.numeric(wald_ci[1:2]), na.last = TRUE)
				if (length(fallback) >= 2L && all(is.finite(fallback)) && fallback[1L] <= est && fallback[2L] >= est && all(abs(fallback[1:2]) <= private$likelihood_ci_max_abs)) {
					ci = fallback
				} else {
					private$cache_nonestimable_se("lik_ratio_confidence_interval_unavailable")
					ci = c(NA_real_, NA_real_)
				}
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			}
			ci
		},

		supports_likelihood_tests = function(){
			TRUE
		},

		supports_information_preference = function(){
			isTRUE(private$supports_likelihood_tests())
		},

		supports_observed_information = function(){
			isTRUE(private$supports_information_preference())
		},

		supports_fisher_information = function(){
			FALSE
		},

		get_supported_testing_types_impl = function(){
			if (isTRUE(private$supports_likelihood_tests())) {
				c("wald", "score", "gradient", "lik_ratio")
			} else {
				"wald"
			}
		},

		get_supported_information_preferences_impl = function(){
			if (isTRUE(private$supports_information_preference())) {
				c("auto", "observed")
			} else {
				"auto"
			}
		},

		get_default_information_source = function(){
			if (isTRUE(private$supports_fisher_information())) return("fisher")
			if (isTRUE(private$supports_observed_information())) return("observed")
			"legacy"
		},

		normalize_testing_type = function(testing_type){
			if (length(testing_type) != 1L) testing_type = testing_type[1L]
			testing_type = tolower(as.character(testing_type))
			switch(
				testing_type,
				wald = "wald",
				score = "score",
				gradient = "gradient",
				lr = "lik_ratio",
				lrt = "lik_ratio",
				lik_ratio = "lik_ratio",
				likelihood_ratio = "lik_ratio",
				stop("testing_type must be one of: wald, score, gradient, lik_ratio", call. = FALSE)
			)
		},
		normalize_information_preference = function(information_preference){
			if (length(information_preference) != 1L) information_preference = information_preference[1L]
			information_preference = tolower(as.character(information_preference))
			switch(
				information_preference,
				auto = "auto",
				fisher = "fisher",
				observed = "observed",
				obs = "observed",
				stop("information_preference must be one of: auto, fisher, observed", call. = FALSE)
			)
		}
	)
)
