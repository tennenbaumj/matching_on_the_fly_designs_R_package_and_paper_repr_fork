#' Mixin for Memoized Likelihood-Test Evaluation
#'
#' A Pattern-1 mixin (plain list with code{$public} and code{$private} slots)
#' providing warm constrained fits, per-null cached likelihood-test components,
#' and score, gradient, likelihood-ratio, and Bartlett p-values. Consumers must
#' provide code{get_likelihood_test_spec()}, likelihood-test cache helpers from
#' code{Inference}, code{get_score_test_information_matrix()}, and the
#' testing-type capability hooks used by the Bartlett paths.
#'
#' Splice into a class with
#' code{private = c(InferenceMixinLikelihoodTestMemoization$private, list(...))}.
#'
#' @keywords internal
#' @noRd
InferenceMixinLikelihoodTestMemoization = list(
	public = list(),
	private = list(
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

		get_memoized_likelihood_test_pval = function(delta, testing_type, spec = NULL, warm_cache_key = NULL, bartlett_B = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec)) {
				stop(class(self)[1], " does not expose a likelihood-test specification.", call. = FALSE)
			}

			bartlett_types = c("lik_ratio_bartlett_approx", "lik_ratio_bartlett_exact")
			need_score = testing_type %in% c("score", "gradient")
			need_information = identical(testing_type, "score")
			need_full_negloglik = testing_type %in% c("lik_ratio", bartlett_types)
			need_null_negloglik = testing_type %in% c("lik_ratio", bartlett_types)
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

			is_bartlett_approx = identical(testing_type, "lik_ratio_bartlett_approx")
			bartlett_B_resolved = if (is_bartlett_approx) as.integer(bartlett_B %||% 99L) else NA_integer_
			if (!is.null(entry$p_value) && is.finite(entry$p_value)) {
				if (!is_bartlett_approx || identical(entry$bartlett_B, bartlett_B_resolved)) {
					return(entry$p_value)
				}
			}

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
			} else if (testing_type == "lik_ratio_bartlett_approx") {
				if (!isTRUE(private$supports_bartlett_likelihood_ratio_approx())) return(NA_real_)
				f_nll = entry$full_negloglik
				n_nll = entry$null_negloglik
				if (is.null(f_nll) || !length(f_nll) || is.null(n_nll) || !length(n_nll)) return(NA_real_)
				if (!is.finite(f_nll) || !is.finite(n_nll)) return(NA_real_)
				factor = tryCatch(
					private$get_bartlett_factor_approx(spec = spec, delta = delta, full_fit = spec$full_fit, null_fit = entry$null_fit, B = bartlett_B_resolved),
					error = function(e) NULL
				)
				if (is.null(factor) || length(factor) != 1L || !is.finite(factor) || factor <= 0) return(NA_real_)
				res = likelihood_ratio_test_from_negloglik_cpp(f_nll, n_nll, df = 1L)
				statistic = as.numeric(res$statistic %||% NA_real_)
				if (!is.finite(statistic)) return(NA_real_)
				p_value = as.numeric(stats::pchisq(statistic / factor, df = 1, lower.tail = FALSE))
			} else if (testing_type == "lik_ratio_bartlett_exact") {
				if (!isTRUE(private$supports_bartlett_likelihood_ratio_exact())) return(NA_real_)
				f_nll = entry$full_negloglik
				n_nll = entry$null_negloglik
				if (is.null(f_nll) || !length(f_nll) || is.null(n_nll) || !length(n_nll)) return(NA_real_)
				if (!is.finite(f_nll) || !is.finite(n_nll)) return(NA_real_)
				factor = tryCatch(
					private$get_bartlett_factor_exact(spec = spec, delta = delta, full_fit = spec$full_fit, null_fit = entry$null_fit),
					error = function(e) NULL
				)
				if (is.null(factor) || length(factor) != 1L || !is.finite(factor) || factor <= 0) return(NA_real_)
				res = likelihood_ratio_test_from_negloglik_cpp(f_nll, n_nll, df = 1L)
				statistic = as.numeric(res$statistic %||% NA_real_)
				if (!is.finite(statistic)) return(NA_real_)
				p_value = as.numeric(stats::pchisq(statistic / factor, df = 1, lower.tail = FALSE))
			} else {
				stop("Unsupported testing_type: ", testing_type, call. = FALSE)
			}

			entry$p_value = p_value
			if (is_bartlett_approx) entry$bartlett_B = bartlett_B_resolved
			private$set_likelihood_test_eval_entry(testing_type, delta, entry)
			p_value
		},

		compute_likelihood_test_two_sided_pval = function(delta, testing_type, bartlett_B = NULL){
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
				warm_cache_key = paste0("likelihood_test:", testing_type),
				bartlett_B = bartlett_B
			)
			if (!is.finite(p_value) && !isTRUE(self$is_nonestimable("estimate"))) {
				private$cache_nonestimable_se(paste0(testing_type, "_test_unavailable"))
			}
			p_value
		}
	)
)
