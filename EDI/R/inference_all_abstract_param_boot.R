#' Parametric-Bootstrap-Capable Likelihood Inference
#'
#' Intermediate abstract base for the subset of likelihood-backed inference
#' families that are plausible targets for parametric null-bootstrap
#' likelihood-ratio calibration.
#'
#' This class sits between \code{InferenceAsympLik} and
#' \code{InferenceAsympLikStdModCache} in the hierarchy.  Families with highly
#' bespoke partial-likelihood, quadrature, frailty, copula, or custom
#' combined-likelihood geometry remain direct children of
#' \code{InferenceAsympLik} and do not pass through here.
#'
#' The only operational user-facing parametric-bootstrap LR methods on this
#' surface are
#' \code{compute_lik_ratio_bootstrap_two_sided_pval(...)} and
#' \code{compute_lik_ratio_bootstrap_confidence_interval(...)}. Diagnostic
#' accessors such as \code{get_last_param_bootstrap_diagnostics()} are
#' supplementary and not alternative execution entry points.
#'
#' Parametric-bootstrap LR calibration is available only for concrete classes
#' that inherit from \code{InferenceParamBootstrap} and whose private method
#' \code{supports_lik_ratio_param_bootstrap()} returns \code{TRUE}. Families
#' that are intentionally unsupported are kept off this branch entirely.
#'
#' @keywords internal
InferenceParamBootstrap = R6::R6Class("InferenceParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(
		#' @description Returns diagnostics from the most recent parametric-bootstrap LR run.
		#' @return A list of diagnostics, or \code{NULL} if no parametric-bootstrap LR run has been executed.
		get_last_param_bootstrap_diagnostics = function(){
			private$cached_values$last_param_bootstrap_diagnostics
		},
		#' @description Bootstrap-calibrated likelihood-ratio two-sided p-value.
		#'
		#' Fits the null model at \code{delta}, simulates \code{B} datasets from
		#' that fitted null, refits unrestricted and null models on each, and
		#' returns the empirical tail probability of the observed LR statistic.
		#' This is the primary user-facing entry point for bootstrap-calibrated
		#' likelihood-ratio p-values.
		#'
		#' This method is available only for classes whose private
		#' \code{supports_lik_ratio_param_bootstrap()} method returns \code{TRUE}.
		#' For unsupported classes it errors immediately rather than silently
		#' falling back to another procedure.
		#'
		#' The standard user-facing arguments are \code{delta = 0},
		#' \code{B = 199}, and \code{show_progress = FALSE}. The remaining
		#' arguments control replicate-quality thresholds and retry behavior.
		#'
		#' Runtime cost is roughly one unrestricted fit plus \code{B} simulated
		#' unrestricted/null refit pairs, so this is typically much more
		#' expensive than the asymptotic LR p-value.
		#'
		#' @param delta        Null treatment effect. Default 0.
		#' @param B            Number of bootstrap replicates. Default 199.
		#' @param show_progress Logical; show a progress bar. Default \code{FALSE}.
		#' @param min_number_usable_samples Minimum number of usable bootstrap
		#'   replicates required to return a finite p-value. Default \code{5L}.
		#' @param max_attempts_per_replicate Maximum number of simulation/refit
		#'   retries per bootstrap replicate. Default \code{2L}.
		#' @return A scalar p-value, or \code{NA_real_} if the computation fails.
		compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L){
			private$active_resampling_operation = "param_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (!isTRUE(private$supports_lik_ratio_param_bootstrap())){
				stop(
					class(self)[1], " does not support parametric-bootstrap LR calibration. ",
					"Override private$supports_lik_ratio_param_bootstrap() and simulate_under_lik_null().",
					call. = FALSE
				)
			}
			if (should_run_asserts()){
				assertNumeric(delta, len = 1)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertCount(max_attempts_per_replicate, positive = TRUE)
				if (as.integer(B) < as.integer(min_number_usable_samples)) {
					stop("B must be at least min_number_usable_samples for bootstrap LR calibration.", call. = FALSE)
				}
			}
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) {
				if (!isTRUE(self$is_nonestimable("estimate")))
					private$cache_nonestimable_se("lik_ratio_bootstrap_spec_unavailable")
				return(NA_real_)
			}

			eval_obs = private$get_memoized_likelihood_test_eval(
				delta         = delta,
				testing_type  = "lik_ratio",
				spec          = spec,
				include_full_negloglik = TRUE,
				include_null_negloglik = TRUE
			)
			if (isTRUE(eval_obs$invalid)) {
				if (!isTRUE(self$is_nonestimable("estimate")))
					private$cache_nonestimable_se("lik_ratio_bootstrap_observed_lr_invalid")
				return(NA_real_)
			}
			if (!is.finite(eval_obs$full_negloglik) || !is.finite(eval_obs$null_negloglik)) {
				if (!isTRUE(self$is_nonestimable("estimate")))
					private$cache_nonestimable_se("lik_ratio_bootstrap_observed_negloglik_nonfinite")
				return(NA_real_)
			}
			lr_obs = 2 * (eval_obs$null_negloglik - eval_obs$full_negloglik)
			if (!is.finite(lr_obs)) {
				if (!isTRUE(self$is_nonestimable("estimate")))
					private$cache_nonestimable_se("lik_ratio_bootstrap_observed_lr_nonfinite")
				return(NA_real_)
			}
			if (private$param_bootstrap_lr_extreme(lr_obs)) {
				if (!isTRUE(self$is_nonestimable("estimate")))
					private$cache_nonestimable_se("lik_ratio_bootstrap_observed_lr_extreme")
				return(NA_real_)
			}
			null_fit = eval_obs$null_fit

			actual_cores = private$effective_parallel_cores("param_bootstrap", self$num_cores)
			if (!is.null(private$seed)) {
				had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				on.exit(
					if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv),
					add = TRUE
				)
				set.seed(private$seed)
			}
			replicate_seeds = sample.int(.Machine$integer.max, as.integer(B), replace = TRUE)

			use_worker_path = isTRUE(private$use_reusable_param_bootstrap_worker())
			deterministic_mode = isTRUE(private$use_deterministic_param_bootstrap())
			run_one_lr = function(b){
				if (isTRUE(deterministic_mode)) {
					private$with_param_bootstrap_thread_budget(1L, {
						private$compute_param_bootstrap_lr_deterministic(
							spec = spec,
							delta = delta,
							null_fit = null_fit,
							seed = replicate_seeds[[b]],
							max_attempts_per_replicate = max_attempts_per_replicate
						)
					})
				} else {
					private$compute_param_bootstrap_lr_impl(
						spec = spec,
						delta = delta,
						null_fit = null_fit,
						seed = replicate_seeds[[b]],
						max_attempts_per_replicate = max_attempts_per_replicate
					)
				}
			}

			run_worker_chunk = function(idxs){
				worker_state = private$create_param_bootstrap_worker_state(spec, delta, null_fit)
				if (is.null(worker_state)) {
					return(lapply(idxs, run_one_lr))
				}
				lapply(idxs, function(idx){
					private$compute_param_bootstrap_worker_lrt(
						worker_state = worker_state,
						delta = delta,
						seed = replicate_seeds[[idx]],
						max_attempts_per_replicate = max_attempts_per_replicate
					)
				})
			}

			results = if (isTRUE(deterministic_mode) && actual_cores <= 1L) {
				lapply(seq_len(B), run_one_lr)
			} else if (isTRUE(deterministic_mode) && actual_cores > 1L) {
				unlist(private$par_lapply(
					as.list(seq_len(B)),
					function(idx) list(run_one_lr(as.integer(idx)[1L])),
					n_cores = actual_cores,
					budget = 1L,
					show_progress = show_progress
				), recursive = FALSE, use.names = FALSE)
			} else if (use_worker_path && actual_cores <= 1L) {
				run_worker_chunk(seq_len(B))
			} else if (!use_worker_path && actual_cores <= 1L) {
				lapply(seq_len(B), run_one_lr)
			} else {
				chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
				chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
				chunks = split(seq_len(B), chunk_id)
				unlist(private$par_lapply(
					chunks,
					if (use_worker_path) run_worker_chunk else function(idxs) lapply(idxs, run_one_lr),
					n_cores = actual_cores,
					budget = 1L,
					show_progress = show_progress
				), recursive = FALSE, use.names = FALSE)
			}

			lr_boots = vapply(results, function(res) as.numeric(res$lr %||% NA_real_)[1L], numeric(1))
			extreme_lr = is.finite(lr_boots) & private$param_bootstrap_lr_extreme(lr_boots)
			if (any(extreme_lr)) {
				for (idx in which(extreme_lr)) {
					results[[idx]] = private$param_boot_failure_result(
						"extreme_lr",
						attempts = results[[idx]]$attempts %||% NA_integer_,
						details = results[[idx]]$details %||% NULL
					)
				}
				lr_boots[extreme_lr] = NA_real_
			}
			finite_lr = lr_boots[is.finite(lr_boots)]
			n_finite = length(finite_lr)
			private$cached_values$last_param_bootstrap_diagnostics =
				private$summarize_param_bootstrap_diagnostics(
					results = results,
					B = as.integer(B),
					min_number_usable_samples = as.integer(min_number_usable_samples),
					max_attempts_per_replicate = as.integer(max_attempts_per_replicate),
					used_reusable_worker = isTRUE(use_worker_path),
					used_deterministic_mode = isTRUE(deterministic_mode)
				)
			private$cached_values$last_param_bootstrap_summary = list(
				B = as.integer(B),
				n_finite = as.integer(n_finite),
				finite_fraction = n_finite / as.integer(B),
				min_number_usable_samples = as.integer(min_number_usable_samples),
				max_attempts_per_replicate = as.integer(max_attempts_per_replicate),
				used_reusable_worker = isTRUE(use_worker_path),
				used_deterministic_mode = isTRUE(deterministic_mode)
			)
			if (n_finite < as.integer(min_number_usable_samples)) {
				private$cache_nonestimable_se("lik_ratio_bootstrap_too_few_converged_samples")
				return(NA_real_)
			}
			n_exceed = sum(finite_lr >= lr_obs)
			min(1, max(2 / (1 + n_finite), (1 + n_exceed) / (1 + n_finite)))
		},

		#' @description Bootstrap-calibrated likelihood-ratio confidence interval.
		#'
		#' Inverts \code{compute_lik_ratio_bootstrap_two_sided_pval} via a
		#' bracket-and-bisect search seeded with the Wald interval.  Each p-value
		#' evaluation costs \code{B} bootstrap refits, so this method is
		#' substantially more expensive than the p-value alone.
		#' This is the primary user-facing entry point for bootstrap-calibrated
		#' likelihood-ratio confidence intervals.
		#'
		#' This method is available only for classes whose private
		#' \code{supports_lik_ratio_param_bootstrap_confidence_interval()} method
		#' returns \code{TRUE}.
		#'
		#' The standard user-facing arguments are \code{B = 199} and
		#' \code{show_progress = FALSE}. Runtime cost is high because each
		#' confidence-interval bound requires repeated bootstrap p-value
		#' evaluations.
		#'
		#' @param alpha        Significance level. Default 0.05.
		#' @param B            Bootstrap replicates per p-value evaluation. Default 199.
		#' @param show_progress Logical; show a progress bar. Default \code{FALSE}.
		#' @param min_number_usable_samples Minimum number of usable bootstrap
		#'   replicates required within each p-value evaluation. Default
		#'   \code{5L}.
		#' @param max_attempts_per_replicate Maximum number of simulation/refit
		#'   retries per bootstrap replicate. Default \code{2L}.
		#' @param root_tolerance Effect-scale tolerance for the inversion root.
		#'   If \code{NULL}, a tolerance proportional to the Wald standard error is
		#'   used. The bootstrap p-value being inverted is Monte Carlo-discrete, so
		#'   solving to machine precision is not meaningful.
		#' @param max_root_iterations Maximum number of bisection iterations per
		#'   bound during interval inversion. Use \code{0L} to return the first
		#'   finite outer bracket. Default \code{8L}.
		#' @return Named two-element numeric vector with the confidence-interval bounds.
		compute_lik_ratio_bootstrap_confidence_interval = function(alpha = 0.05, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L, root_tolerance = NULL, max_root_iterations = 8L){
			private$active_resampling_operation = "param_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (!isTRUE(private$supports_lik_ratio_param_bootstrap_confidence_interval())){
				stop(
					class(self)[1], " does not support parametric-bootstrap LR confidence intervals.",
					call. = FALSE
				)
			}
			if (should_run_asserts()){
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertCount(max_attempts_per_replicate, positive = TRUE)
				assertNumeric(root_tolerance, lower = .Machine$double.xmin, null.ok = TRUE)
				assertIntegerish(max_root_iterations, lower = 0, len = 1, any.missing = FALSE)
			}
			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			pval_fn = function(d) self$compute_lik_ratio_bootstrap_two_sided_pval(
				d,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate
			)

			se = tryCatch(private$get_standard_error(), error = function(e) NA_real_)
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			root_tolerance = if (is.null(root_tolerance)) max(1e-4, 0.10 * step) else as.numeric(root_tolerance)[1L]
			if (!is.finite(root_tolerance) || root_tolerance <= 0) root_tolerance = 1e-4
			max_root_iterations = max(0L, as.integer(max_root_iterations)[1L])
			pval_mc_tolerance = max(1 / (as.numeric(B) + 1), sqrt(alpha * (1 - alpha) / as.numeric(B)))
			wald_ci = tryCatch(
				private$compute_wald_confidence_interval_impl(alpha),
				error = function(e) c(NA_real_, NA_real_)
			)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else NA_real_
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else NA_real_

			p_est = pval_fn(est)
			if (!is.finite(p_est)) {
				ci = c(NA_real_, NA_real_)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}
			if (p_est < alpha) {
				ci = c(est, est)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}

			finite_seed_offsets = abs(c(lower_seed, upper_seed) - est)
			finite_seed_offsets = finite_seed_offsets[is.finite(finite_seed_offsets)]
			max_radius = max(
				8 * step,
				if (length(finite_seed_offsets)) 2 * max(finite_seed_offsets) else 0
			)

			find_bound = function(direction, seed){
				outer = NA_real_
				f_outer = NA_real_
				if (is.finite(seed) && ((direction < 0 && seed < est) || (direction > 0 && seed > est))) {
					if (abs(seed - est) <= max_radius) {
						f_seed = pval_fn(seed) - alpha
						if (is.finite(f_seed)) {
							if (abs(f_seed) <= 1e-6) return(seed)
							if (f_seed <= 0 && abs(f_seed) <= pval_mc_tolerance) return(seed)
							if (f_seed <= 0) {
								outer = seed
								f_outer = f_seed
							}
						}
					}
				}

				if (!is.finite(outer)) {
					for (i in 0:19) {
						d = est + direction * step * 2^i
						if (abs(d - est) > max_radius) break
						f_d = pval_fn(d) - alpha
						if (!is.finite(f_d)) next
						if (abs(f_d) <= 1e-6) return(d)
						if (f_d <= 0 && abs(f_d) <= pval_mc_tolerance) return(d)
						if (f_d <= 0) {
							outer = d
							f_outer = f_d
							break
						}
					}
				}

				if (!is.finite(outer) || !is.finite(f_outer)) return(NA_real_)

				f_est = p_est - alpha
				if (!is.finite(f_est)) return(NA_real_)
				if (abs(f_est) <= 1e-6) return(est)

				lower = min(est, outer)
				upper = max(est, outer)
				f.lower = if (isTRUE(all.equal(lower, outer))) f_outer else f_est
				f.upper = if (isTRUE(all.equal(upper, outer))) f_outer else f_est
				if (!is.finite(f.lower) || !is.finite(f.upper)) return(NA_real_)
				if (f.lower * f.upper > 0) return(NA_real_)
				if (max_root_iterations == 0L) return(outer)

				best = (lower + upper) / 2
				for (iter in seq_len(max_root_iterations)) {
					mid = (lower + upper) / 2
					best = mid
					if (abs(upper - lower) <= root_tolerance) return(mid)
					f.mid = pval_fn(mid) - alpha
					if (!is.finite(f.mid)) next
					if (abs(f.mid) <= max(1e-6, pval_mc_tolerance)) return(mid)
					if ((f.mid >= 0 && f.lower >= 0) || (f.mid <= 0 && f.lower <= 0)) {
						lower = mid
						f.lower = f.mid
					} else {
						upper = mid
						f.upper = f.mid
					}
				}
				best
			}

			ci = c(
				find_bound(direction = -1, seed = lower_seed),
				find_bound(direction = 1, seed = upper_seed)
			)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		param_bootstrap_extreme_lr_threshold = 1e6,
		param_bootstrap_lr_extreme = function(lr, max_abs = private$param_bootstrap_extreme_lr_threshold){
			lr = as.numeric(lr)
			max_abs = as.numeric(max_abs)[1L]
			if (!is.finite(max_abs) || max_abs <= 0) max_abs = 1e6
			is.finite(lr) & abs(lr) > max_abs
		},
		simulate_param_boot_bernoulli_y = function(mu){
			mu = as.numeric(mu)
			if (!length(mu) || any(!is.finite(mu))) return(NULL)
			mu = pmin(pmax(mu, 0), 1)
			as.numeric(stats::rbinom(length(mu), 1L, mu))
		},
		simulate_param_boot_poisson_y = function(mu){
			mu = as.numeric(mu)
			if (!length(mu) || any(!is.finite(mu)) || any(mu < 0)) return(NULL)
			as.numeric(stats::rpois(length(mu), mu))
		},
		simulate_param_boot_gaussian_y = function(mu, sigma2){
			mu = as.numeric(mu)
			sigma2 = as.numeric(sigma2)[1L]
			if (!length(mu) || any(!is.finite(mu)) || !is.finite(sigma2) || sigma2 <= 0) return(NULL)
			as.numeric(mu + stats::rnorm(length(mu), 0, sqrt(sigma2)))
		},
		simulate_param_boot_ordinal_y = function(X, params_null, y_template, cdf_fn){
			X = as.matrix(X)
			params_null = as.numeric(params_null)
			cat_vals = sort(unique(as.numeric(y_template)))
			K = length(cat_vals)
			if (!is.function(cdf_fn) || K < 2L) return(NULL)
			n_alpha = K - 1L
			if (length(params_null) <= n_alpha) return(NULL)
			thresholds = params_null[seq_len(n_alpha)]
			betas = params_null[(n_alpha + 1L):length(params_null)]
			eta = as.numeric(X %*% betas)
			cum_probs = outer(thresholds, eta, function(a, e) cdf_fn(a - e))
			cat_probs = pmax(rbind(cum_probs, 1) - rbind(0, cum_probs), 0)
			y_sim = cat_vals[apply(cat_probs, 2, function(p){
				s = sum(p)
				if (!is.finite(s) || s <= 0) return(1L)
				sample.int(K, 1L, prob = p / s)
			})]
			y_sim = as.numeric(y_sim)
			if (length(unique(y_sim)) < K) return(NULL)
			y_sim
		},
		simulate_param_boot_weibull_observed = function(X, b_null, log_sigma, y_obs, dead_obs){
			X = as.matrix(X)
			b_null = as.numeric(b_null)
			log_sigma = as.numeric(log_sigma)[1L]
			y_obs = as.numeric(y_obs)
			dead_obs = as.numeric(dead_obs)
			sigma = exp(log_sigma)
			if (!is.finite(sigma) || sigma <= 0) return(NULL)
			mu_log = as.numeric(X %*% b_null)
			T_sim = stats::rweibull(nrow(X), shape = 1 / sigma, scale = exp(mu_log))
			if (!all(is.finite(T_sim)) || any(T_sim <= 0)) return(NULL)
			C_i = ifelse(dead_obs == 0, y_obs, Inf)
			y_sim = pmin(T_sim, C_i)
			dead_sim = as.numeric(T_sim <= C_i)
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			list(y = y_sim, dead = dead_sim)
		},
		supports_reusable_param_bootstrap_worker = function(){
			TRUE
		},
		use_reusable_param_bootstrap_worker = function(){
			isTRUE(private$reusable_bootstrap_worker_enabled) &&
				private$has_private_method("supports_reusable_param_bootstrap_worker") &&
				isTRUE(tryCatch(private$supports_reusable_param_bootstrap_worker(), error = function(e) FALSE))
		},
		use_deterministic_param_bootstrap = function(){
			!is.null(private$seed) && is.finite(private$seed)
		},
		with_param_bootstrap_thread_budget = function(budget, expr){
			budget = max(1L, as.integer(budget)[1L])
			ns = asNamespace("EDI")
			edi_env = ns$edi_env
			prev_override = edi_env$num_cores_override
			prev_threads = getOption(".edi_last_set_threads")
			if (is.null(prev_threads) || length(prev_threads) != 1L || !is.finite(prev_threads)) {
				prev_threads = 1L
			}
			assign("num_cores_override", budget, envir = edi_env)
			ns$set_package_threads(budget)
			on.exit({
				assign("num_cores_override", prev_override, envir = edi_env)
				ns$set_package_threads(prev_threads)
			}, add = TRUE)
			force(expr)
		},
		with_param_bootstrap_seed = function(seed, expr){
			if (is.null(seed) || !is.finite(seed)) return(force(expr))
			old_rng_kind = RNGkind()
			had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
			if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
			on.exit({
				do.call(RNGkind, as.list(old_rng_kind))
				if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv)
			}, add = TRUE)
			RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
			set.seed(as.integer(seed)[1L])
			force(expr)
		},
		param_boot_failure_result = function(reason, attempts = 1L, details = NULL){
			list(
				success = FALSE,
				lr = NA_real_,
				reason = as.character(reason)[1L],
				attempts = as.integer(attempts)[1L],
				details = details
			)
		},
		param_boot_success_result = function(lr, attempts = 1L){
			list(
				success = TRUE,
				lr = as.numeric(lr)[1L],
				reason = "success",
				attempts = as.integer(attempts)[1L],
				details = NULL
			)
		},
		validate_param_bootstrap_spec = function(boot_spec){
			if (is.null(boot_spec) || !is.list(boot_spec)) return(FALSE)
			if (is.null(boot_spec$full_fit)) return(FALSE)
			if (!is.function(boot_spec$fit_null)) return(FALSE)
			if (!is.function(boot_spec$neg_loglik)) return(FALSE)
			TRUE
		},
		validate_param_bootstrap_worker_data = function(worker_state, worker_data){
			if (is.null(worker_state) || is.null(worker_data) || !is.list(worker_data)) return(FALSE)
			y = worker_data$y
			if (is.null(y) || !is.numeric(y) || length(y) != as.integer(worker_state$n)) return(FALSE)
			if (!all(is.finite(y))) return(FALSE)
			if (!is.null(worker_data$dead)) {
				dead = as.numeric(worker_data$dead)
				if (length(dead) != as.integer(worker_state$n) || any(!is.finite(dead))) return(FALSE)
			}
			TRUE
		},
		extract_param_bootstrap_failure_reason = function(boot_spec, default = "simulated_data_failure"){
			if (is.null(boot_spec) || !is.list(boot_spec)) return(default)
			reason = boot_spec$failure_reason %||% attr(boot_spec, "edi_param_boot_failure_reason", exact = TRUE)
			if (is.null(reason) || !length(reason) || !is.character(reason)) default else reason[[1L]]
		},
		compute_param_bootstrap_lr_from_boot_spec = function(boot_spec, delta){
			if (is.null(boot_spec)) {
				return(private$param_boot_failure_result("simulated_data_failure"))
			}
			if (!isTRUE(private$validate_param_bootstrap_spec(boot_spec))) {
				return(private$param_boot_failure_result(private$extract_param_bootstrap_failure_reason(boot_spec, default = "full_refit_failure")))
			}
			full_nll_boot = tryCatch(
				boot_spec$neg_loglik(boot_spec$full_fit),
				error = function(e) NA_real_
			)
			if (!is.finite(full_nll_boot)) return(private$param_boot_failure_result("full_refit_failure"))
			null_fit_boot = tryCatch(
				boot_spec$fit_null(delta),
				error = function(e) NULL
			)
			if (is.null(null_fit_boot)) return(private$param_boot_failure_result("null_refit_failure"))
			null_nll_boot = tryCatch(
				boot_spec$neg_loglik(null_fit_boot),
				error = function(e) NA_real_
			)
			if (!is.finite(null_nll_boot)) return(private$param_boot_failure_result("non_finite_lr"))
			lr_boot = 2 * (null_nll_boot - full_nll_boot)
			if (is.finite(lr_boot)) private$param_boot_success_result(lr_boot) else private$param_boot_failure_result("non_finite_lr")
		},
		compute_param_bootstrap_lr_impl = function(spec, delta, null_fit, seed = NULL, max_attempts_per_replicate = 1L, worker_priv = NULL){
			runner_priv = worker_priv %||% private
			private$with_param_bootstrap_seed(seed, {
				last_result = private$param_boot_failure_result("simulated_data_failure", attempts = 0L)
				for (attempt in seq_len(max(1L, as.integer(max_attempts_per_replicate)))) {
					boot_spec = tryCatch(
						runner_priv$simulate_under_lik_null(spec, delta, null_fit),
						error = function(e) NULL
					)
					res = runner_priv$compute_param_bootstrap_lr_from_boot_spec(boot_spec, delta)
					res$attempts = as.integer(attempt)
					last_result = res
					if (isTRUE(res$success) && is.finite(res$lr)) return(res)
				}
				last_result
			})
		},
		compute_param_bootstrap_lr_deterministic = function(spec, delta, null_fit, seed = NULL, max_attempts_per_replicate = 1L){
			worker_state = private$create_param_bootstrap_worker_state(spec, delta, null_fit)
			if (!is.null(worker_state)) {
				return(private$compute_param_bootstrap_worker_lrt(
					worker_state = worker_state,
					delta = delta,
					seed = seed,
					max_attempts_per_replicate = max_attempts_per_replicate
				))
			}
			worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			worker$num_cores = 1L
			worker_priv = worker$.__enclos_env__$private
			worker_priv$cached_values = list()
			worker_priv$clear_likelihood_null_warm_cache()
			worker_priv$clear_fit_warm_start()
			worker_spec = worker_priv$get_likelihood_test_spec()
			if (is.null(worker_spec)) {
				return(private$param_boot_failure_result("simulated_data_failure"))
			}
			worker_eval = worker_priv$get_memoized_likelihood_test_eval(
				delta = delta,
				testing_type = "lik_ratio",
				spec = worker_spec,
				include_full_negloglik = TRUE,
				include_null_negloglik = TRUE
			)
			if (isTRUE(worker_eval$invalid) || is.null(worker_eval$null_fit)) {
				return(private$param_boot_failure_result("null_refit_failure"))
			}
			worker_priv$compute_param_bootstrap_lr_impl(
				spec = worker_spec,
				delta = delta,
				null_fit = worker_eval$null_fit,
				seed = seed,
				max_attempts_per_replicate = max_attempts_per_replicate,
				worker_priv = worker_priv
			)
		},
		load_param_bootstrap_draw_into_worker = function(worker_state, sim_data){
			if (is.null(worker_state) || is.null(worker_state$worker_priv) || is.null(worker_state$state_env)) {
				return(FALSE)
			}
			worker_priv = worker_state$worker_priv
			worker_data = sim_data$worker_data %||% NULL
			if (!private$validate_param_bootstrap_worker_data(worker_state, worker_data)) {
				worker_state$state_env$current_param_bootstrap_draw = sim_data
				return(FALSE)
			}
			worker_des_priv = worker_state$worker_des_priv
			worker_priv$X = worker_state$base_X
			worker_priv$w = worker_state$base_w
			worker_priv$m = worker_state$base_m
			worker_priv$n = as.integer(worker_state$n)
			worker_priv$y = as.numeric(worker_data$y)
			worker_priv$dead = if (is.null(worker_data$dead)) NULL else as.numeric(worker_data$dead)
			worker_priv$any_censoring = !is.null(worker_priv$dead) && any(worker_priv$dead == 0)
			worker_priv$y_temp = worker_priv$y
			extra_names = setdiff(names(worker_data), c("y", "dead"))
			for (nm in extra_names) {
				worker_priv[[nm]] = worker_data[[nm]]
			}
			worker_priv$cached_values = list()
			worker_priv$clear_likelihood_null_warm_cache()
			worker_priv$clear_fit_warm_start()
			# Parametric bootstrap holds X and w fixed across draws (lines above restore both
			# to base_X / base_w).  cached_design_matrix = [1|w|X] and cached_hardened_X_cov
			# therefore remain valid and are deliberately preserved, avoiding an O(Np^2)
			# drop_highly_correlated_cols recomputation on every draw.
			# Compare: randomization preserves cached_hardened_X_cov for the same reason
			# (X_cov fixed) but must clear cached_design_matrix because w changes.
			worker_priv$reduced_design_keep_cache = NULL
			worker_priv$fixed_covariate_keep_cache = NULL
			worker_priv$best_X_colnames = NULL
			worker_priv$best_Xmm_colnames = NULL
			worker_priv$cached_mod = NULL
			worker_priv$fit_warm_start = worker_state$base_fit_warm_start
			worker_priv$fit_warm_start_type = worker_state$base_fit_warm_start_type
			worker_priv$fit_warm_start_fisher = worker_state$base_fit_warm_start_fisher
			if (!is.null(worker_des_priv)) {
				worker_des_priv$Xraw = worker_state$base_Xraw
				worker_des_priv$Ximp = worker_state$base_Ximp
				worker_des_priv$X = worker_state$base_X
				worker_des_priv$w = worker_state$base_w
				worker_des_priv$y = worker_priv$y
				worker_des_priv$dead = worker_priv$dead
				if (!is.null(worker_state$base_m)) worker_des_priv$m = worker_state$base_m
				worker_des_priv$n = as.integer(worker_state$n)
				worker_des_priv$t = as.integer(worker_state$n)
				worker_des_priv$all_subject_data_cache = list()
				worker_des_priv$lin_centered_covariates = NULL
				if (is.function(worker_des_priv$reset_matching_caches)) worker_des_priv$reset_matching_caches()
			}
			worker_state$state_env$current_param_bootstrap_draw = sim_data
			TRUE
		},
		create_param_bootstrap_worker_state = function(spec, delta, null_fit){
			worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			worker$num_cores = 1L
			worker_priv = worker$.__enclos_env__$private
			worker_des = if (!is.null(worker_priv$des_obj)) worker_priv$des_obj$duplicate(verbose = FALSE) else NULL
			worker_des_priv = if (!is.null(worker_des)) worker_des$.__enclos_env__$private else NULL
			source_des_priv = private$des_obj_priv_int
			if (!is.null(worker_des)) {
				worker_priv$des_obj = worker_des
				worker_priv$des_obj_priv_int = worker_des_priv
			}
			worker_priv$cached_values = list()
			worker_priv$clear_likelihood_null_warm_cache()
			worker_priv$clear_fit_warm_start()
			worker_priv$X = if (!is.null(private$X)) private$X else private$get_X()
			worker_spec = worker_priv$get_likelihood_test_spec()
			if (is.null(worker_spec)) return(NULL)
			worker_eval = worker_priv$get_memoized_likelihood_test_eval(
				delta = delta,
				testing_type = "lik_ratio",
				spec = worker_spec,
				include_full_negloglik = TRUE,
				include_null_negloglik = TRUE
			)
			if (isTRUE(worker_eval$invalid) || is.null(worker_eval$null_fit)) return(NULL)
			list(
				worker = worker,
				worker_priv = worker_priv,
				worker_des_priv = worker_des_priv,
				spec = worker_spec,
				null_fit = worker_eval$null_fit,
				base_fit_warm_start = private$fit_warm_start,
				base_fit_warm_start_type = private$fit_warm_start_type,
				base_fit_warm_start_fisher = private$fit_warm_start_fisher,
				base_Xraw = if (!is.null(source_des_priv$Xraw)) source_des_priv$Xraw else NULL,
				base_Ximp = if (!is.null(source_des_priv$Ximp)) source_des_priv$Ximp else NULL,
				base_X = if (!is.null(private$X)) private$X else private$get_X(),
				base_w = if (!is.null(source_des_priv$w)) as.numeric(source_des_priv$w) else NULL,
				base_y = if (!is.null(source_des_priv$y)) source_des_priv$y else NULL,
				base_dead = if (!is.null(source_des_priv$dead)) as.numeric(source_des_priv$dead) else NULL,
				base_m = if (!is.null(source_des_priv$m)) source_des_priv$m else NULL,
				n = private$n,
				state_env = list2env(list(
					current_param_bootstrap_draw = NULL
				), parent = emptyenv())
			)
		},
		compute_param_bootstrap_worker_lrt = function(worker_state, delta, seed = NULL, max_attempts_per_replicate = 1L){
			if (is.null(worker_state) || is.null(worker_state$worker_priv)) {
				return(private$param_boot_failure_result("simulated_data_failure"))
			}
			private$with_param_bootstrap_seed(seed, {
				last_result = private$param_boot_failure_result("simulated_data_failure", attempts = 0L)
				for (attempt in seq_len(max(1L, as.integer(max_attempts_per_replicate)))) {
					boot_spec = tryCatch(
						worker_state$worker_priv$simulate_under_lik_null(worker_state$spec, delta, worker_state$null_fit),
						error = function(e) NULL
					)
					loaded_into_worker = private$load_param_bootstrap_draw_into_worker(worker_state, boot_spec)
					if (isTRUE(loaded_into_worker)) {
						worker_spec = tryCatch(worker_state$worker_priv$get_likelihood_test_spec(), error = function(e) NULL)
						if (is.null(worker_spec)) {
							res = private$param_boot_failure_result("full_refit_failure")
						} else {
							worker_eval = tryCatch(
								worker_state$worker_priv$get_memoized_likelihood_test_eval(
									delta = delta,
									testing_type = "lik_ratio",
									spec = worker_spec,
									include_full_negloglik = TRUE,
									include_null_negloglik = TRUE
								),
								error = function(e) NULL
							)
							if (is.null(worker_eval) || isTRUE(worker_eval$invalid)) {
								res = private$param_boot_failure_result("null_refit_failure")
							} else if (!is.finite(worker_eval$full_negloglik)) {
								res = private$param_boot_failure_result("full_refit_failure")
							} else if (!is.finite(worker_eval$null_negloglik)) {
								res = private$param_boot_failure_result("null_refit_failure")
							} else {
								lr_boot = 2 * (worker_eval$null_negloglik - worker_eval$full_negloglik)
								res = if (is.finite(lr_boot)) private$param_boot_success_result(lr_boot) else private$param_boot_failure_result("non_finite_lr")
							}
						}
					} else {
						res = worker_state$worker_priv$compute_param_bootstrap_lr_from_boot_spec(
							worker_state$state_env$current_param_bootstrap_draw,
							delta
						)
					}
					res$attempts = as.integer(attempt)
					last_result = res
					if (isTRUE(res$success) && is.finite(res$lr)) return(res)
				}
				last_result
			})
		},
		summarize_param_bootstrap_diagnostics = function(results, B, min_number_usable_samples, max_attempts_per_replicate, used_reusable_worker, used_deterministic_mode = FALSE){
			if (is.null(results)) results = list()
			reasons = vapply(results, function(res) as.character(res$reason %||% "unknown_failure")[1L], character(1))
			attempts = vapply(results, function(res) as.integer(res$attempts %||% NA_integer_)[1L], integer(1))
			lrs = vapply(results, function(res) as.numeric(res$lr %||% NA_real_)[1L], numeric(1))
			success = is.finite(lrs)
			reason_levels = c("success", "simulated_data_failure", "full_refit_failure", "null_refit_failure", "non_finite_lr", "extreme_lr", "unknown_failure")
			reason_factor = factor(ifelse(reasons %in% reason_levels, reasons, "unknown_failure"), levels = reason_levels)
			reason_counts = as.list(as.integer(table(reason_factor)))
			names(reason_counts) = reason_levels
			list(
				B = as.integer(B),
				n_success = sum(success),
				n_failure = sum(!success),
				success_fraction = mean(success),
				min_number_usable_samples = as.integer(min_number_usable_samples),
				max_attempts_per_replicate = as.integer(max_attempts_per_replicate),
				used_reusable_worker = isTRUE(used_reusable_worker),
				used_deterministic_mode = isTRUE(used_deterministic_mode),
				reason_counts = reason_counts,
				prop_simulated_data_failure = unname(reason_counts$simulated_data_failure) / max(1L, as.integer(B)),
				prop_full_refit_failure = unname(reason_counts$full_refit_failure) / max(1L, as.integer(B)),
				prop_null_refit_failure = unname(reason_counts$null_refit_failure) / max(1L, as.integer(B)),
				prop_non_finite_lr = unname(reason_counts$non_finite_lr) / max(1L, as.integer(B)),
				mean_attempts = mean(attempts, na.rm = TRUE),
				replicate_results = results
			)
		},
		supports_lik_ratio_param_bootstrap = function() FALSE,
		supports_lik_ratio_param_bootstrap_confidence_interval = function(){
			isTRUE(private$supports_lik_ratio_param_bootstrap())
		},
		#' Approximate (Monte-Carlo) Bartlett support automatically follows
		#' parametric-bootstrap LR support: any family that already implements
		#' simulate_under_lik_null() gets the shared Monte-Carlo factor below for
		#' free. Families that need to withhold Bartlett support for a reason that
		#' doesn't affect parametric-bootstrap LR itself (e.g. a known raw-LR
		#' miscalibration) should override this method directly.
		supports_bartlett_likelihood_ratio_approx = function(){
			isTRUE(private$supports_lik_ratio_param_bootstrap())
		},
		bartlett_factor_mc_min_usable_fraction = 0.2,
		bartlett_factor_mc_max_attempts_per_replicate = 2L,
		#' Generic Monte-Carlo Bartlett correction factor, shared by every
		#' InferenceParamBootstrap family that already implements
		#' simulate_under_lik_null() for parametric-bootstrap LR calibration.
		#'
		#' Simulates B datasets under H0: (tested parameter) = delta using the same
		#' simulate_under_lik_null()/refit machinery as
		#' compute_lik_ratio_bootstrap_two_sided_pval(), and sets
		#'   c(delta) = mean(simulated LR statistics)
		#' which is a Monte-Carlo estimate of E[LR | H0], the quantity a classical
		#' analytic Bartlett correction targets exactly (E[LR]/df = 1 under a
		#' chi-square(df) reference). This is deliberately model-agnostic: any
		#' family wanting a faster closed-form factor instead can override
		#' get_bartlett_factor_exact()/supports_bartlett_likelihood_ratio_exact()
		#' (a subclass override of this method wins too, if preferred).
		get_bartlett_factor_approx = function(spec, delta, full_fit, null_fit, B = 99){
			if (!isTRUE(private$supports_lik_ratio_param_bootstrap())) return(NULL)
			if (is.null(null_fit)) return(NULL)

			B = as.integer(B)
			# Scales with B rather than a fixed count, so small-B calls (e.g. quick
			# smoke tests) aren't structurally impossible to satisfy.
			min_usable = max(2L, as.integer(ceiling(private$bartlett_factor_mc_min_usable_fraction * B)))
			max_attempts = as.integer(private$bartlett_factor_mc_max_attempts_per_replicate)

			if (!is.null(private$seed)) {
				had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				on.exit(
					if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv),
					add = TRUE
				)
				set.seed(private$seed)
			}
			replicate_seeds = sample.int(.Machine$integer.max, B, replace = TRUE)

			lr_boots = vapply(seq_len(B), function(b){
				res = tryCatch(
					private$compute_param_bootstrap_lr_impl(
						spec = spec,
						delta = delta,
						null_fit = null_fit,
						seed = replicate_seeds[[b]],
						max_attempts_per_replicate = max_attempts
					),
					error = function(e) NULL
				)
				as.numeric(res$lr %||% NA_real_)[1L]
			}, numeric(1))

			extreme = is.finite(lr_boots) & private$param_bootstrap_lr_extreme(lr_boots)
			lr_boots[extreme] = NA_real_
			finite_lr = lr_boots[is.finite(lr_boots)]
			if (length(finite_lr) < min_usable) return(NULL)

			factor = mean(finite_lr)
			if (!is.finite(factor) || factor <= 0) return(NULL)
			factor
		},
		#' Simulate a bootstrap dataset under the fitted null likelihood and return
		#' a minimal spec list for refitting.
		#'
		#' Must be overridden by families that set supports_lik_ratio_param_bootstrap()
		#' to TRUE.  The returned list must contain at least:
		#'   - full_fit:  unrestricted fit on the simulated data
		#'   - fit_null:  function(delta, start) returning a constrained fit
		#'   - neg_loglik: function(fit) returning the neg-log-likelihood
		simulate_under_lik_null = function(spec, delta, null_fit){
			NULL
		}
	)
)
