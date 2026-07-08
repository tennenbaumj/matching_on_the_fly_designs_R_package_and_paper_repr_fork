#' Randomization-based Inference
#'
#' Abstract class for randomization-based inference.
#'
#' @keywords internal
InferenceRand = R6::R6Class("InferenceRand",
	inherit = Inference,
	lock_objects = FALSE,
	public = list(
		#' @description Set Custom Randomization Statistic Computation
		#' @param custom_randomization_statistic_function  A function that returns a scalar value.
		set_custom_randomization_statistic_function = function(custom_randomization_statistic_function){
			if (!is.null(custom_randomization_statistic_function) && !is.null(private[["compiled_cpp_stat_fn"]])) {
				stop("Cannot specify both custom_randomization_statistic_function and custom_randomization_statistic_cpp.")
			}
			if (should_run_asserts()) {
				assertFunction(custom_randomization_statistic_function, null.ok = TRUE)
			}
			private[["custom_randomization_statistic_function"]] = custom_randomization_statistic_function
			private$cached_values$t0s_rand = NULL
			private$cached_values$rand_distr_cache = list()
			private$cached_values$custom_stat_analysis = NULL
		},
		#' @description Set Custom Randomization Statistic as a Pre-Compiled Rcpp Function
		#' @param fn  A compiled Rcpp function returning a scalar \code{double}. Must accept
		#'   either \code{(NumericVector y, IntegerVector w)} or
		#'   \code{(NumericVector y, IntegerVector w, IntegerVector dead)} as arguments.
		#'   Compile with \code{Rcpp::cppFunction()} before passing. Pass \code{NULL} to clear.
		#'   Cannot be combined with \code{set_custom_randomization_statistic_function}.
		set_custom_randomization_statistic_cpp = function(fn){
			if (!is.null(fn) && !is.null(private[["custom_randomization_statistic_function"]])) {
				stop("Cannot specify both custom_randomization_statistic_function and custom_randomization_statistic_cpp.")
			}
			if (!is.null(fn)) {
				if (!is.function(fn)) stop("custom_randomization_statistic_cpp must be a compiled Rcpp function, not a ", class(fn)[1], ".")
				arity = length(formals(fn))
				if (!arity %in% c(2L, 3L)) stop("custom_randomization_statistic_cpp must accept 2 arguments (y, w) or 3 arguments (y, w, dead); got ", arity, ".")
			}
			private[["compiled_cpp_stat_fn"]] = fn
			private$cached_values$t0s_rand = NULL
			private$cached_values$rand_distr_cache = list()
			private$cached_values$custom_stat_analysis = NULL
		},
		#' @description Computes the randomization distribution of the treatment effect estimate under the sharp null.
		#'
		#' @param r  					Number of randomization vectors. Default 501.
		#' @param delta  				The null difference. Default 0.
		#' @param transform_responses  Type of transformation. Default "none".
		#' @param show_progress  		Show progress bar. Default TRUE.
		#' @param permutations  		Pre-computed permutations. Default NULL.
		#' @param debug  				If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics including error messages, warning messages, counts of each,
		#'   and summary proportions for iterations with errors, warnings, and illegal (non-finite)
		#'   values. Runs serially. Default \code{FALSE}.
		#' @return 	When \code{debug = FALSE} (default), a numeric vector of length \code{r}. When
		#'   \code{debug = TRUE}, a list with: \code{values}, \code{errors} (list of character
		#'   vectors, one per iteration), \code{warnings} (list of character vectors, one per
		#'   iteration), \code{num_errors}, \code{num_warnings},
		#'   \code{prop_iterations_with_errors}, \code{prop_iterations_with_warnings}, and
		#'   \code{prop_illegal_values}.
		#' @param zero_one_logit_clamp The clamping amount for exact 0 and 1 values when logging
		approximate_randomization_distribution_beta_hat_T = function(r = 501, delta = 0, transform_responses = "none", show_progress = TRUE, permutations = NULL, debug = FALSE, zero_one_logit_clamp = .Machine$double.eps){
			private$active_resampling_operation = "rand"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Randomization inference")
				assertNumeric(delta); assertCount(r, positive = TRUE); assertFlag(debug)
			}
			if (is.null(permutations)) permutations = private$generate_permutations(r)
			setup = private$setup_randomization_template_and_shifts(delta, transform_responses, zero_one_logit_clamp)
			if (!isTRUE(debug) && !is.null(permutations) && private$has_private_method("compute_fast_randomization_distr")) {
				fast_distr = tryCatch(
					private$compute_fast_randomization_distr(setup$template$.__enclos_env__$private$y, permutations, delta, transform_responses, zero_one_logit_clamp),
					error = function(e) NULL
				)
				if (!is.null(fast_distr)) return(fast_distr)
				# If fast path threw, fall through to the standard reusable-worker path below.
			}
			if (!isTRUE(debug) && !is.null(permutations) &&
				isTRUE(private$use_reusable_bootstrap_worker()) &&
				is.null(private[["custom_randomization_statistic_function"]]) &&
				is.null(private[["compiled_cpp_stat_fn"]])) {
				actual_rand_cores = private$effective_parallel_cores("rand_pval", self$num_cores)
				return(private$compute_randomization_distr_via_reused_worker_states(
					permutations = permutations,
					delta = delta,
					transform_responses = transform_responses,
					actual_rand_cores = actual_rand_cores,
					show_progress = show_progress,
					setup = setup,
					zero_one_logit_clamp = zero_one_logit_clamp
				))
			}
			custom_stat_analysis = private$analyze_custom_randomization_statistic()
			use_lightweight_custom_stat = isTRUE(custom_stat_analysis$can_use_lightweight_yw_only)
			use_perms = !is.null(permutations) && (!is.null(permutations$w_mat) || length(permutations) >= r)
			need_thread_objs = !(use_lightweight_custom_stat && use_perms)
			inf_template = if (need_thread_objs) self$duplicate() else NULL
			des_template = if (need_thread_objs) setup$template$duplicate() else NULL
			# Warm up the design template cache if it is a sequential design that uses covariates.
			if (!is.null(des_template) && isTRUE(des_template$.__enclos_env__$private$uses_covariates)) {
				is_verbose = isTRUE(private$verbose)
				if (is_verbose) cat("Warming up design cache... ")
				tryCatch({
					priv = des_template$.__enclos_env__$private
					old_t = priv$t
					if (is.null(priv$all_subject_data_cache)) priv$all_subject_data_cache = list()
					n_subjects = des_template$get_n()
					for (t_temp in 1 : n_subjects) {
						priv$t = t_temp
						priv$compute_all_subject_data()
					}
					priv$t = old_t
					if (is_verbose) cat("done.\n")
				}, error = function(e) {
					if (exists("old_t")) des_template$.__enclos_env__$private$t = old_t
					if (is_verbose) cat("failed.\n")
				})
			}
			if (!is.null(inf_template) && private$has_match_structure && private$object_has_private_method(inf_template, "compute_basic_match_data"))
				inf_template$.__enclos_env__$private$compute_basic_match_data()
			if (isTRUE(debug)) {
				debug_results = if (isTRUE(private$use_reusable_bootstrap_worker()) && is.null(private$custom_randomization_statistic_function) && is.null(private[["compiled_cpp_stat_fn"]])){
					# Fast path: use reused workers
					worker_state = private$create_bootstrap_worker_state()
					cleanup_worker = private$cleanup_bootstrap_worker_state
					if (is.function(cleanup_worker)) on.exit(cleanup_worker(worker_state), add = TRUE)
					lapply(seq_len(r), function(idx){
						iter_warns = character(0)
						iter_result = withCallingHandlers(
							tryCatch({
								# Generate permuted weights
								perm_w = if (use_perms) {
									if (!is.null(permutations$w_mat)) {
										j = ((idx - 1L) %% ncol(permutations$w_mat)) + 1L
										permutations$w_mat[, j]
									} else {
										perm_data = permutations[[idx]]
										if (is.list(perm_data) && !is.null(perm_data$w)) perm_data$w else perm_data
									}
								} else {
									sample(private$w)
								}
								# Load into worker
								private$load_randomization_perm_into_worker(worker_state, perm_w, delta, transform_responses, setup$y_delta, setup$base_template_y, setup$base_template_dead, zero_one_logit_clamp)
								# Compute estimate
								list(val = private$compute_bootstrap_worker_estimate(worker_state))
							}, error = function(e) list(val = NA_real_, error = conditionMessage(e))),
							warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
						)
						list(
							val = as.numeric(iter_result$val)[1L],
							errors = if (!is.null(iter_result$error)) iter_result$error else character(0),
							warnings = iter_warns
						)
					})
				} else {
					# Standard path: duplicate objects (slow)
					lapply(seq_len(r), function(idx) {
						iter_warns = character(0)
						iter_result = withCallingHandlers(
							tryCatch({
								worker_des = if (!is.null(des_template)) setup$template$duplicate() else NULL
								worker_inf = if (!is.null(inf_template)) self$duplicate(verbose = FALSE, make_fork_cluster = FALSE) else NULL
								private$run_randomization_iteration(worker_des, worker_inf, if (use_perms) idx else NULL, permutations, delta, transform_responses, setup$y_delta, setup$base_template_y, setup$base_template_dead, custom_stat_analysis, setup$lightweight_custom_context, debug = TRUE, zero_one_logit_clamp = zero_one_logit_clamp)
							}, error = function(e) list(val = NA_real_, error = conditionMessage(e))),
							warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
						)
						list(
							val = as.numeric(iter_result$val)[1L],
							errors = if (!is.null(iter_result$error)) iter_result$error else character(0),
							warnings = iter_warns
						)
					})
				}
				debug_results = debug_results[!vapply(debug_results, is.null, logical(1))]
				if (length(debug_results) == 0L) {
					stop("All randomization iterations failed or returned invalid results. Check for worker crashes or out-of-memory issues.")
				}
				values = sapply(debug_results, `[[`, "val")
				errors_list = lapply(debug_results, `[[`, "errors")
				warnings_list = lapply(debug_results, `[[`, "warnings")
				num_errors_vec = lengths(errors_list)
				num_warnings_vec = lengths(warnings_list)
				return(list(
					values = values,
					errors = errors_list,
					warnings = warnings_list,
					num_errors = num_errors_vec,
					num_warnings = num_warnings_vec,
					prop_iterations_with_errors = mean(num_errors_vec > 0),
					prop_iterations_with_warnings = mean(num_warnings_vec > 0),
					prop_illegal_values = mean(!is.finite(values))
				))
			}
			actual_rand_cores = private$effective_parallel_cores("rand_pval", self$num_cores)
			if (actual_rand_cores > 1L && need_thread_objs) {
				do_warmup_iter = function() {
					w_des = if (!is.null(des_template)) des_template$duplicate() else NULL
					w_inf = if (!is.null(inf_template)) inf_template$duplicate(make_fork_cluster = FALSE) else NULL
					private$run_randomization_iteration(w_des, w_inf, if(use_perms) 1L else NULL, permutations, delta, transform_responses, setup$y_delta, setup$base_template_y, setup$base_template_dead, custom_stat_analysis, setup$lightweight_custom_context, zero_one_logit_clamp = zero_one_logit_clamp)
				}
				# Run warmup TWICE and use the second timing. The first call often pays
				# cold-start penalties (C++ JIT, OS page-cache misses, R bytecode compilation)
				# that inflate the estimate 5–15× vs steady-state cost, causing the guard to
				# wrongly choose parallel for small r values like r = 19.
				system.time(do_warmup_iter())  # First call: discarded (cold-start overhead)
				t_rand_warmup = system.time(do_warmup_iter())[[3]]  # Second call: representative cost
				# Existing cluster: ~10ms round-trip overhead. No cluster yet: ~300ms lazy creation.
				fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.3
				if (t_rand_warmup * r < fork_overhead_estimate * actual_rand_cores * 2.0) {
					actual_rand_cores = 1L
				}
			} else if (actual_rand_cores > 1L && !need_thread_objs) {
				# Use warmup timing for the lightweight path, same guard as the thread-obj path above.
				do_warmup_iter_lw = function() {
					private$run_randomization_iteration(
						NULL, NULL,
						if (use_perms) 1L else NULL,
						permutations, delta, transform_responses,
						setup$y_delta, setup$base_template_y, setup$base_template_dead,
						custom_stat_analysis, setup$lightweight_custom_context,
						zero_one_logit_clamp = zero_one_logit_clamp
					)
				}
				system.time(do_warmup_iter_lw())
				t_lw_warmup = system.time(do_warmup_iter_lw())[[3]]
				fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.3
				if (t_lw_warmup * r < fork_overhead_estimate * actual_rand_cores * 2.0) {
					actual_rand_cores = 1L
				}
			}
			beta_hat_T_diff_ws = unlist(private$par_lapply(1:r, function(idx) {
				suppressWarnings({
					worker_des = if (!is.null(des_template)) des_template$duplicate() else NULL
					worker_inf = if (!is.null(inf_template)) inf_template$duplicate(make_fork_cluster = FALSE) else NULL
					private$run_randomization_iteration(worker_des, worker_inf, if(use_perms) idx else NULL, permutations, delta, transform_responses, setup$y_delta, setup$base_template_y, setup$base_template_dead, custom_stat_analysis, setup$lightweight_custom_context, zero_one_logit_clamp = zero_one_logit_clamp)
				})
			}, n_cores = actual_rand_cores, show_progress = show_progress,
			export_list = list(
				des_template = des_template,
				inf_template = inf_template,
				permutations = permutations,
				delta = delta,
				setup = setup,
				custom_stat_analysis = custom_stat_analysis,
				use_perms = use_perms,
				zero_one_logit_clamp = zero_one_logit_clamp
			)))
			if (!is.numeric(beta_hat_T_diff_ws)) beta_hat_T_diff_ws = as.numeric(beta_hat_T_diff_ws)
			beta_hat_T_diff_ws
		},
		#' @description Computes a randomization-based p-value.
		#' @param r  	Number of randomization vectors.
		#' @param delta  				Null difference.
		#' @param transform_responses  Transformation.
		#' @param na.rm 				Remove NAs.
		#' @param show_progress  	Show progress.
		#' @param permutations  	Pre-computed permutations.
		#' @param zero_one_logit_clamp The clamping amount for exact 0 and 1 values when logging
		#' @return 	Randomization p-value.
		compute_rand_two_sided_pval = function(r = 501, delta = 0, transform_responses = "none", na.rm = TRUE, show_progress = TRUE, permutations = NULL, zero_one_logit_clamp = .Machine$double.eps){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Randomization inference")
				assertLogical(na.rm)
				if (private$des_obj_priv_int$response_type == "incidence" && is.null(private$custom_randomization_statistic_function)) stop("Randomization tests are not supported for incidence. Use Zhang method.")
			}
			if (is.null(permutations)) permutations = private$generate_permutations(r)
			if (identical(transform_responses, "none")) {
				transform_responses = switch(
					private$des_obj_priv_int$response_type,
					continuous = "none",
					proportion = "logit",
					count = "log",
					survival = "log",
					"none"
				)
			}
			cache_key = private$build_randomization_distribution_cache_key(r, delta, transform_responses, permutations)
			if (transform_responses == "none" && is.null(private[["custom_randomization_statistic_function"]]) && !is.null(private$cached_values$t0s_rand) && length(private$cached_values$t0s_rand) >= r) {
				t0s = private$cached_values$t0s_rand[seq_len(r)] + delta
				t = private$compute_treatment_estimate_during_randomization_inference()
				if (is.function(self$is_nonestimable) && isTRUE(self$is_nonestimable("estimate"))) return(NA_real_)
				if (length(t) != 1 || !is.finite(t)) return(NA_real_)
				na_t0s = !is.finite(t0s)
				nsim_adj = sum(!na_t0s)
				if (nsim_adj == 0L) return(NA_real_)
				return(min(1, max(2 / nsim_adj, 2 * min(sum(t0s >= t, na.rm = TRUE) / nsim_adj, sum(t0s <= t, na.rm = TRUE) / nsim_adj))))
			}
			if (is.null(private$cached_values$rand_distr_cache)) private$cached_values$rand_distr_cache = list()
			t = if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) {
				custom_stat_analysis = private$analyze_custom_randomization_statistic()
				if (isTRUE(custom_stat_analysis$can_use_lightweight_yw_only)) {
					private$evaluate_lightweight_custom_randomization_statistic(
						private$des_obj_priv_int,
						private$y,
						private$w,
						private$dead
					)
				} else {
					private$custom_randomization_statistic_function()
				}
			} else {
				private$compute_treatment_estimate_during_randomization_inference()
			}
			if (is.function(self$is_nonestimable) && isTRUE(self$is_nonestimable("estimate"))) return(NA_real_)
			if (length(t) != 1 || !is.finite(t)) return(NA_real_)
			mc_pval = private$compute_two_sided_pval_with_sequential_mc(
				t = t,
				r = r,
				delta = delta,
				transform_responses = transform_responses,
				show_progress = show_progress,
				permutations = permutations,
				cache_key = cache_key,
				zero_one_logit_clamp = zero_one_logit_clamp
			)
			if (!is.null(mc_pval)) return(mc_pval)
			t0s = private$get_randomization_distribution_prefix(
				r = r,
				delta = delta,
				transform_responses = transform_responses,
				show_progress = show_progress,
				permutations = permutations,
				cache_key = cache_key,
				zero_one_logit_clamp = zero_one_logit_clamp
			)
			private$compute_two_sided_randomization_pval_from_t0s(t0s, t)
		}
	),
	private = list(
		custom_randomization_statistic_function = NULL,
		compiled_cpp_stat_fn = NULL,
		randomization_mc_control = NULL,
		normalize_delta_for_cache = function(delta, resolution = NULL){
			if (!is.finite(delta)) return("NA")
			if (!is.null(resolution) && is.finite(resolution) && resolution > 0) {
				delta = round(as.numeric(delta) / resolution) * resolution
			}
			format(as.numeric(delta), scientific = TRUE, digits = 17)
		},
		compute_randomization_distr_via_reused_worker_states = function(permutations, delta, transform_responses, actual_rand_cores, show_progress, setup, zero_one_logit_clamp) {
			nsim = if (!is.null(permutations$w_mat)) ncol(permutations$w_mat) else length(permutations)
			if (!isTRUE(nsim > 0L)) return(numeric(0))
			get_perm_w = if (!is.null(permutations$w_mat)) {
				w_mat_local = permutations$w_mat
				function(i) w_mat_local[, i]
			} else {
				function(i) {
					p = permutations[[i]]
					if (is.list(p) && !is.null(p$w)) p$w else p
				}
			}
			chunk_n = max(1L, min(as.integer(actual_rand_cores), nsim))
			chunk_id = ceiling(seq_len(nsim) / ceiling(nsim / chunk_n))
			chunks = split(seq_len(nsim), chunk_id)
			run_chunk = function(idxs) {
				worker_state = private$create_bootstrap_worker_state()
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					perm_w = get_perm_w(idxs[k])
					out[k] = tryCatch({
						private$load_randomization_perm_into_worker(
							worker_state, perm_w, delta, transform_responses,
							setup$y_delta, setup$base_template_y, setup$base_template_dead,
							zero_one_logit_clamp
						)
						as.numeric(private$compute_bootstrap_worker_estimate(worker_state))[1L]
					}, error = function(e) NA_real_)
					# Sequential null anchoring: after each successful permutation, update
					# base_fit_warm_start to the converged parameters so the next permutation
					# starts from the current null-distribution point rather than the MLE.
					# Iterative models (logistic, Poisson, NegBin, ordinal, survival, etc.) call
					# set_fit_warm_start() inside their fitting function, which writes the converged
					# params to inf_priv$fit_warm_start.  We copy that here.
					# Cold objects (fit_warm_start_enabled = FALSE) have set_fit_warm_start() as a
					# no-op, so inf_priv$fit_warm_start stays NULL and no update occurs.
					if (is.finite(out[k])) {
						inf_priv_seq = if (!is.null(worker_state$worker_inf)) {
							worker_state$worker_inf$.__enclos_env__$private
						} else if (!is.null(worker_state$worker_priv)) {
							worker_state$worker_priv
						} else if (!is.null(worker_state$worker)) {
							worker_state$worker$.__enclos_env__$private
						} else NULL
						if (!is.null(inf_priv_seq)) {
							new_ws = inf_priv_seq$fit_warm_start
							if (!is.null(new_ws) && length(new_ws) > 0L && all(is.finite(new_ws))) {
								worker_state$base_fit_warm_start        = new_ws
								worker_state$base_fit_warm_start_type   = inf_priv_seq$fit_warm_start_type
								# Do NOT carry fit_warm_start_fisher: the Fisher information is
								# X_full'WX_full where X_full = [1 | w | X_cov].  The treatment
								# column w changes every permutation, invalidating all cross-terms
								# involving w.  Force a fresh recompute from the new design matrix.
								worker_state$base_fit_warm_start_fisher = NULL
							}
						}
					}
				}
				out
			}
			if (actual_rand_cores <= 1L) return(as.numeric(run_chunk(seq_len(nsim))))
			as.numeric(unlist(private$par_lapply(
				chunks,
				run_chunk,
				n_cores = actual_rand_cores,
				budget = 1L,
				show_progress = show_progress
			), use.names = FALSE))
		},
		build_fast_randomization_worker_cache = function(prev_cache = NULL, preserve_cache_keys = character()){
			cache = list()
			if (is.null(prev_cache)) {
				cache$rand_distr_cache = list()
				return(cache)
			}
			always_keep = c("m_cache", "t0s_rand", "custom_stat_analysis")
			for (nm in unique(c(always_keep, preserve_cache_keys))) {
				if (!is.null(prev_cache[[nm]])) cache[[nm]] = prev_cache[[nm]]
			}
			cache$rand_distr_cache = list()
			cache
		},
		compute_fast_randomization_distr_via_reused_worker = function(y, permutations, delta, transform_responses, preserve_cache_keys = character(), zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (is.null(permutations)) return(NULL)
			nsim = if (!is.null(permutations$w_mat)) ncol(permutations$w_mat) else length(permutations)
			if (!isTRUE(nsim > 0L)) return(numeric(0))
			get_perm_data = if (!is.null(permutations$w_mat)) {
				w_mat = permutations$w_mat
				m_mat = permutations$m_mat
				function(i) {
					list(
						w = w_mat[, i],
						m_vec = if (!is.null(m_mat)) m_mat[, i] else NULL
					)
				}
			} else {
				function(i) permutations[[i]]
			}
			actual_rand_cores = min(private$effective_parallel_cores("rand_pval", self$num_cores), nsim)
			chunk_n = max(1L, min(as.integer(actual_rand_cores), nsim))
			chunk_id = ceiling(seq_len(nsim) / ceiling(nsim / chunk_n))
			chunks = split(seq_len(nsim), chunk_id)
			run_chunk = function(idxs) {
				worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
				worker$num_cores = 1L
				w_priv = worker$.__enclos_env__$private
				worker_des = if (!is.null(w_priv$des_obj)) w_priv$des_obj$duplicate(verbose = FALSE) else NULL
				if (!is.null(worker_des)) private$sync_randomization_worker_state(worker_des, worker)
				worker_des_priv = if (!is.null(worker_des)) worker_des$.__enclos_env__$private else NULL
				base_m = w_priv$m
				base_cache = w_priv$cached_values
				w_priv$y = as.numeric(y)
				w_priv$y_temp = w_priv$y
				if (!is.null(worker_des_priv)) {
					worker_des_priv$y = w_priv$y
					worker_des_priv$dead = w_priv$dead
					if (!is.null(base_m)) worker_des_priv$m = base_m
					private$sync_randomization_worker_state(worker_des, worker)
				}
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					perm_data = get_perm_data(idxs[k])
					if (!is.null(worker_des_priv)) {
						worker_des_priv$w = as.integer(perm_data$w)
						worker_des_priv$m = if (!is.null(perm_data$m_vec)) perm_data$m_vec else base_m
						y_sim = w_priv$y_temp
						if (delta != 0) {
							resp_type = worker_des_priv$response_type
							if (transform_responses == "logit") {
								y_sim[perm_data$w == 1] = inv_logit(logit(y_sim[perm_data$w == 1], zero_one_logit_clamp) + delta, zero_one_logit_clamp)
							} else if (transform_responses == "log" && resp_type == "survival") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else if (transform_responses == "log" && resp_type == "count") {
								y_sim[perm_data$w == 1] = as.integer(round(y_sim[perm_data$w == 1] * exp(delta)))
							} else if (transform_responses == "log" && resp_type != "count") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] + delta
							}
						}
						worker_des_priv$y = y_sim
						worker_des_priv$dead = w_priv$dead
						private$sync_randomization_worker_state(worker_des, worker)
					} else {
						w_priv$w = as.integer(perm_data$w)
						w_priv$m = if (!is.null(perm_data$m_vec)) perm_data$m_vec else base_m
						y_sim = w_priv$y_temp
						if (delta != 0) {
							resp_type = w_priv$des_obj_priv_int$response_type
							if (transform_responses == "logit") {
								y_sim[perm_data$w == 1] = inv_logit(logit(y_sim[perm_data$w == 1], zero_one_logit_clamp) + delta, zero_one_logit_clamp)
							} else if (transform_responses == "log" && resp_type == "survival") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else if (transform_responses == "log" && resp_type == "count") {
								y_sim[perm_data$w == 1] = as.integer(round(y_sim[perm_data$w == 1] * exp(delta)))
							} else if (transform_responses == "log" && resp_type != "count") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] + delta
							}
						}
						w_priv$y = y_sim
					}
					w_priv$cached_values = private$build_fast_randomization_worker_cache(
						if (k == 1L) base_cache else w_priv$cached_values,
						preserve_cache_keys = preserve_cache_keys
					)
					est = tryCatch(
						w_priv$compute_treatment_estimate_during_randomization_inference(estimate_only = TRUE),
						error = function(e) NA_real_
					)
					if (is.function(worker$is_nonestimable) &&
					    isTRUE(worker$is_nonestimable("estimate"))) {
						est = NA_real_
					}
					if (is.list(est) && "b" %in% names(est)) est = est$b[1]
					out[k] = as.numeric(est)[1]
				}
				out
			}
			as.numeric(unlist(private$par_lapply(
				chunks,
				run_chunk,
				n_cores = actual_rand_cores,
				budget = 1L,
				show_progress = FALSE,
				export_list = list(
					permutations = permutations,
					y = y,
					transform_responses = transform_responses,
					preserve_cache_keys = preserve_cache_keys
				)
			), use.names = FALSE))
		},
		compute_two_sided_randomization_pval_from_t0s = function(t0s, t){
			na_t0s = !is.finite(t0s)
			nsim_adj = sum(!na_t0s)
			if (nsim_adj == 0L) return(NA_real_)
			min(1, max(2 / nsim_adj, 2 * min(sum(t0s >= t, na.rm = TRUE) / nsim_adj, sum(t0s <= t, na.rm = TRUE) / nsim_adj)))
		},
		compute_two_sided_randomization_pval_band = function(t0s, t, conf_level){
			valid = is.finite(t0s)
			n = sum(valid)
			if (n == 0L) return(c(NA_real_, NA_real_))
			x_ge = sum(t0s[valid] >= t)
			x_le = sum(t0s[valid] <= t)
			binom_band = function(x){
				alpha_band = 1 - conf_level
				lower = if (x <= 0L) 0 else stats::qbeta(alpha_band / 2, x, n - x + 1)
				upper = if (x >= n) 1 else stats::qbeta(1 - alpha_band / 2, x + 1, n - x)
				c(lower, upper)
			}
			band_ge = binom_band(x_ge)
			band_le = binom_band(x_le)
			band = c(2 * min(band_ge[1], band_le[1]), 2 * min(band_ge[2], band_le[2]))
			pmin(1, pmax(0, band))
		},
		subset_permutations = function(permutations, indices){
			if (is.null(permutations)) return(NULL)
			if (!is.null(permutations$w_mat)) {
				list(
					w_mat = permutations$w_mat[, indices, drop = FALSE],
					m_mat = if (!is.null(permutations$m_mat)) permutations$m_mat[, indices, drop = FALSE] else NULL
				)
			} else {
				permutations[indices]
			}
		},
		get_randomization_distribution_prefix = function(r, delta, transform_responses, show_progress, permutations, cache_key, batch_size = NULL, zero_one_logit_clamp = .Machine$double.eps){
			if (is.null(private$cached_values$rand_distr_cache)) private$cached_values$rand_distr_cache = list()
			cached = if (!is.null(cache_key)) private$cached_values$rand_distr_cache[[cache_key]] else NULL
			if (length(cached) > 0L && !any(is.finite(cached))) {
				cached = NULL
				if (!is.null(cache_key)) private$cached_values$rand_distr_cache[[cache_key]] = NULL
			}
			have = length(cached)
			target = if (is.null(batch_size)) as.integer(r) else min(as.integer(r), have + as.integer(batch_size))
			if (have < target) {
				idx = seq.int(have + 1L, target)
				new_t0s = self$approximate_randomization_distribution_beta_hat_T(
					r = length(idx),
					delta = delta,
					transform_responses = transform_responses,
					show_progress = isTRUE(show_progress) && target >= r && have == 0L,
					permutations = private$subset_permutations(permutations, idx),
					zero_one_logit_clamp = zero_one_logit_clamp
				)
				cached = c(cached, new_t0s)
				if (!is.null(cache_key)) private$cached_values$rand_distr_cache[[cache_key]] = cached
			}
			cached[seq_len(target)]
		},
		compute_two_sided_pval_with_sequential_mc = function(t, r, delta, transform_responses, show_progress, permutations, cache_key, zero_one_logit_clamp = .Machine$double.eps){
			mc_control = private$randomization_mc_control
			if (is.null(mc_control) || !isTRUE(mc_control$mc_enable) || !is.finite(mc_control$mc_stop_threshold)) return(NULL)
			batch_size = min(as.integer(r), as.integer(mc_control$mc_batch_size))
			min_draws = min(as.integer(r), as.integer(mc_control$mc_min_draws))
			if (batch_size <= 0L || min_draws <= 0L || batch_size >= as.integer(r)) return(NULL)
			conf_level = mc_control$mc_conf_level
			threshold = mc_control$mc_stop_threshold
			repeat {
				t0s = private$get_randomization_distribution_prefix(
					r = r,
					delta = delta,
					transform_responses = transform_responses,
					show_progress = FALSE,
					permutations = permutations,
					cache_key = cache_key,
					batch_size = batch_size,
					zero_one_logit_clamp = zero_one_logit_clamp
				)
				n_valid = sum(is.finite(t0s))
				p_hat = private$compute_two_sided_randomization_pval_from_t0s(t0s, t)
				if (length(t0s) >= as.integer(r) || n_valid < min_draws || !is.finite(p_hat)) {
					if (length(t0s) >= as.integer(r) || !is.finite(p_hat)) return(p_hat)
				} else {
					band = private$compute_two_sided_randomization_pval_band(t0s, t, conf_level)
					if (is.finite(band[1]) && is.finite(band[2]) && (band[2] < threshold || band[1] > threshold)) return(p_hat)
				}
				if (length(t0s) >= as.integer(r)) return(p_hat)
			}
		},
		generate_permutations = function(r){
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
			}
			design_sig = private$stable_signature(list(
				class = class(private$des_obj),
				n = private$n,
				prob_T = private$prob_T,
				m = private$des_obj_priv_int$m,
				strata_cols = private$des_obj_priv_int$strata_cols
			))
			cache_key = paste0(as.integer(r), "|", design_sig)
			cached = private$des_obj_priv_int$permutations_cache[[cache_key]]
			if (!is.null(cached)) return(cached)
			des_template = private$des_obj$duplicate()
			w_mat = des_template$draw_ws_according_to_design(as.integer(r))
			if (!is.matrix(w_mat)) {
				w_mat = matrix(as.numeric(w_mat), nrow = private$n)
			}
			# draw_ws_according_to_design returns {-1,+1}; convert to internal {0,1} for injection.
			w_mat = (w_mat + 1L) / 2L
			storage.mode(w_mat) = "numeric"
			permutations = list(
				w_mat = w_mat,
				m_mat = NULL
			)
			private$des_obj_priv_int$permutations_cache[[cache_key]] = permutations
			permutations
		},
		build_randomization_distribution_cache_key = function(r, delta, transform_responses, permutations){
			delta_key = formatC(as.numeric(delta), digits = 17, format = "fg", flag = "#")
			perm_sig = private$stable_signature(permutations)
			paste(as.integer(r), delta_key, transform_responses, perm_sig, sep = "|")
		},
		shift_randomization_responses = function(y, w, delta, transform_responses, response_type, inverse = FALSE, zero_one_logit_clamp = .Machine$double.eps){
			if (delta == 0) return(y)
			y_shifted = y
			idx_treated = which(w == 1)
			if (length(idx_treated) == 0L) return(y_shifted)
			signed_delta = if (isTRUE(inverse)) -delta else delta
			if (transform_responses == "logit") {
				y_shifted[idx_treated] = inv_logit(logit(y_shifted[idx_treated], zero_one_logit_clamp) + signed_delta, zero_one_logit_clamp)
				return(y_shifted)
			}
			if (transform_responses == "log" && response_type == "survival") {
				y_shifted[idx_treated] = y_shifted[idx_treated] * exp(signed_delta)
				return(y_shifted)
			}
			if (transform_responses == "log" && response_type == "count") {
				y_shifted[idx_treated] = as.integer(round(y_shifted[idx_treated] * exp(signed_delta)))
				return(y_shifted)
			}
			if (transform_responses == "log" && response_type != "count") {
				y_shifted[idx_treated] = y_shifted[idx_treated] * exp(signed_delta)
				return(y_shifted)
			}
			y_shifted[idx_treated] = y_shifted[idx_treated] + signed_delta
			y_shifted
		},
		setup_randomization_template_and_shifts = function(delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			# Use the design matrix and response vector from the design object.
			template = private$des_obj$duplicate()
			y_delta = template$.__enclos_env__$private$y
			if (delta != 0){
				if (should_run_asserts()) {
					if (private$des_obj_priv_int$response_type == "incidence" && is.null(private$custom_randomization_statistic_function)) stop("randomization tests with delta nonzero not supported for incidence")
				}
				template$.__enclos_env__$private$y = private$shift_randomization_responses(
					y = template$.__enclos_env__$private$y,
					w = private$w,
					delta = delta,
					transform_responses = transform_responses,
					response_type = private$des_obj_priv_int$response_type,
					inverse = TRUE,
					zero_one_logit_clamp = zero_one_logit_clamp
				)
				y_delta = template$.__enclos_env__$private$y
			}
			list(template = template, y_delta = y_delta, base_template_y = private$y, base_template_dead = private$dead, lightweight_custom_context = private$des_obj_priv_int)
		},
		load_randomization_perm_into_worker = function(worker_state, perm_w, delta, transform_responses, y_delta, base_template_y, base_template_dead, zero_one_logit_clamp = .Machine$double.eps){
			inf_priv = if (!is.null(worker_state$worker_inf)) {
				worker_state$worker_inf$.__enclos_env__$private
			} else if (!is.null(worker_state$worker_priv)) {
				worker_state$worker_priv
			} else {
				worker_state$worker$.__enclos_env__$private
			}
			des_priv = if (!is.null(worker_state$worker_des)) {
				worker_state$worker_des$.__enclos_env__$private
			} else {
				worker_state$worker_des_priv
			}
			
			# Update design private state
			if (!is.null(des_priv)) des_priv$w = perm_w
			if (delta != 0) {
				y_sim = private$shift_randomization_responses(
					y = y_delta,
					w = perm_w,
					delta = delta,
					transform_responses = transform_responses,
					response_type = if (!is.null(des_priv)) des_priv$response_type else private$des_obj_priv_int$response_type,
					inverse = FALSE,
					zero_one_logit_clamp = zero_one_logit_clamp
				)
				if (!is.null(des_priv)) des_priv$y = y_sim
			} else {
				y_sim = base_template_y
				if (!is.null(des_priv)) des_priv$y = y_sim
			}
			
			# Sync to inference private state
			inf_priv$w = perm_w
			inf_priv$y = y_sim
			inf_priv$y_temp = y_sim
			inf_priv$dead = if (!is.null(des_priv)) des_priv$dead else base_template_dead
			inf_priv$cached_values$KKstats = NULL # reset
			inf_priv$cached_values$beta_hat_T = NULL
			inf_priv$cached_values$s_beta_hat_T = NULL
			inf_priv$likelihood_null_warm_cache = list()
			
			# Reset all private design matrix and covariate caches
			inf_priv$cached_design_matrix = NULL
			inf_priv$cached_w_for_design_matrix = NULL
			inf_priv$cached_harden_for_design_matrix = NULL
			# inf_priv$cached_hardened_X_cov = NULL  # Preserve covariate-only cache under randomization
			inf_priv$cached_reduced_X = NULL
			inf_priv$cached_X_full_for_reduced = NULL
			inf_priv$cached_keep_for_reduced = NULL
			inf_priv$cached_j_treat_for_reduced = NULL
			
			inf_priv$fit_warm_start = worker_state$base_fit_warm_start
			inf_priv$fit_warm_start_type = worker_state$base_fit_warm_start_type
			inf_priv$fit_warm_start_fisher = worker_state$base_fit_warm_start_fisher
			
			if (!is.null(inf_priv$compute_basic_match_data)) inf_priv$compute_basic_match_data()
			invisible(NULL)
		},
		sync_randomization_worker_state = function(thread_des_obj, thread_inf_obj){
			if (is.null(thread_des_obj) || is.null(thread_inf_obj)) return(invisible(NULL))
			des_priv = thread_des_obj$.__enclos_env__$private
			inf_priv = thread_inf_obj$.__enclos_env__$private
			inf_priv$des_obj = thread_des_obj
			inf_priv$des_obj_priv_int = des_priv
			inf_priv$w = des_priv$w
			inf_priv$y = des_priv$y
			inf_priv$y_temp = des_priv$y
			inf_priv$dead = des_priv$dead
			if (private$has_match_structure) inf_priv$m = des_priv$m
			if (!is.null(inf_priv$compute_basic_match_data)) inf_priv$compute_basic_match_data()
			
			# Reset all private design matrix and covariate caches
			inf_priv$cached_design_matrix = NULL
			inf_priv$cached_w_for_design_matrix = NULL
			inf_priv$cached_harden_for_design_matrix = NULL
			# inf_priv$cached_hardened_X_cov = NULL  # Preserve covariate-only cache under randomization
			inf_priv$cached_reduced_X = NULL
			inf_priv$cached_X_full_for_reduced = NULL
			inf_priv$cached_keep_for_reduced = NULL
			inf_priv$cached_j_treat_for_reduced = NULL
			
			invisible(NULL)
		},
		run_randomization_iteration = function(thread_des_obj, thread_inf_obj, perm_idx, permutations, delta, transform_responses, y_delta, base_template_y, base_template_dead, custom_stat_analysis, lightweight_custom_context, debug = FALSE, zero_one_logit_clamp = .Machine$double.eps){
			use_perms = !is.null(perm_idx)
			get_perm_data = if (use_perms) {
				if (!is.null(permutations$w_mat)) {
					n_avail = ncol(permutations$w_mat)
					function(i) { j = ((i - 1L) %% n_avail) + 1L; list(w = permutations$w_mat[, j], m_vec = if (!is.null(permutations$m_mat)) permutations$m_mat[, j] else NULL) }
				} else function(i) permutations[[i]]
			} else NULL
			if (isTRUE(custom_stat_analysis$can_use_lightweight_yw_only) && use_perms) {
				perm_data = get_perm_data(perm_idx); w_sim = perm_data$w; y_sim = y_delta
				if (delta != 0) {
					y_sim = private$shift_randomization_responses(
						y = y_sim,
						w = w_sim,
						delta = delta,
						transform_responses = transform_responses,
						response_type = lightweight_custom_context$response_type,
						inverse = FALSE,
						zero_one_logit_clamp = zero_one_logit_clamp
					)
				}
				val = private$evaluate_lightweight_custom_randomization_statistic(lightweight_custom_context, y_sim, w_sim, base_template_dead)
				if (isTRUE(debug)) return(list(val = val, error = NULL))
				return(val)
			}
			if (use_perms) {
				perm_data = get_perm_data(perm_idx)
				thread_des_obj$.__enclos_env__$private$w = perm_data$w
				if (!is.null(perm_data$m_vec)) thread_des_obj$.__enclos_env__$private$m = perm_data$m_vec
			} else {
				thread_des_obj$.__enclos_env__$private$resample_assignment()
			}
			if (delta != 0) {
				y_sim = private$shift_randomization_responses(
					y = y_delta,
					w = thread_des_obj$.__enclos_env__$private$w,
					delta = delta,
					transform_responses = transform_responses,
					response_type = thread_des_obj$.__enclos_env__$private$response_type,
					inverse = FALSE,
					zero_one_logit_clamp = zero_one_logit_clamp
				)
				thread_des_obj$.__enclos_env__$private$y = y_sim
			}
			private$sync_randomization_worker_state(thread_des_obj, thread_inf_obj)
			iter_error = NULL
			estimate = tryCatch(
				thread_inf_obj$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference(estimate_only = TRUE),
				error = function(e) { iter_error <<- conditionMessage(e); NA_real_ }
			)
			if (is.function(thread_inf_obj$is_nonestimable) &&
			    isTRUE(thread_inf_obj$is_nonestimable("estimate"))) {
				estimate = NA_real_
			}
			val = if (is.list(estimate) && "b" %in% names(estimate)) as.numeric(estimate$b[1]) else as.numeric(estimate)
			if (isTRUE(debug)) return(list(val = val, error = iter_error))
			val
		},
		get_compiled_cpp_stat = function() private[["compiled_cpp_stat_fn"]],
		analyze_custom_randomization_statistic = function(){
			if (!is.null(private$cached_values$custom_stat_analysis)) return(private$cached_values$custom_stat_analysis)
			if (is.null(private$custom_randomization_statistic_function) && is.null(private[["compiled_cpp_stat_fn"]])) {
				analysis = list(can_use_lightweight_yw_only = FALSE, needs_match_data = TRUE)
				private$cached_values$custom_stat_analysis = analysis; return(analysis)
			}
			# C++ stat is always lightweight: it only ever receives (y, w) or (y, w, dead).
			if (!is.null(private[["compiled_cpp_stat_fn"]])) {
				analysis = list(can_use_lightweight_yw_only = TRUE, needs_match_data = FALSE)
				private$cached_values$custom_stat_analysis = analysis; return(analysis)
			}
			# Basic analysis: does it only use y and w?
			body_str = paste(deparse(body(private$custom_randomization_statistic_function)), collapse = " ")
			# Look for access to other members of private$des_obj_priv_int
			can_use_lightweight = !grepl("private\\$des_obj_priv_int\\$(?!y|w|dead)", body_str, perl = TRUE)
			analysis = list(can_use_lightweight_yw_only = can_use_lightweight, needs_match_data = FALSE)
			private$cached_values$custom_stat_analysis = analysis
			analysis
		},
		evaluate_lightweight_custom_randomization_statistic = function(lightweight_custom_context, y, w, dead){
			# Fast path: compiled C++ function — no R interpreter overhead per permutation.
			cpp_fn = private$get_compiled_cpp_stat()
			if (!is.null(cpp_fn)) {
				arity = length(formals(cpp_fn))
				return(as.numeric(
					if (arity >= 3L) cpp_fn(y, as.integer(w), as.integer(dead))
					else cpp_fn(y, as.integer(w))
				)[1L])
			}
			# We simulate the environment for the custom statistic
			fn = private$custom_randomization_statistic_function
			old_env = environment(fn)
			on.exit(environment(fn) <- old_env, add = TRUE)
			eval_env = new.env(parent = environment(fn))
			
			private_proxy = new.env(parent = emptyenv())
			seq_priv_proxy = new.env(parent = emptyenv())
			seq_priv_proxy$y = y; seq_priv_proxy$w = w; seq_priv_proxy$dead = dead
			private_proxy$des_obj_priv_int = seq_priv_proxy
			
			eval_env$private = private_proxy
			eval_env$inf_priv = private_proxy
			eval_env$des_priv = seq_priv_proxy
			eval_env$des_obj_priv_int = seq_priv_proxy
			
			environment(fn) = eval_env
			eval_env$.custom_randomization_statistic_function = fn
			eval(quote(.custom_randomization_statistic_function()), envir = eval_env)
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (identical(private$des_obj_priv_int$response_type, "proportion") &&
			    (inherits(self, "InferenceAbstractKKQuantileRegrIVWC") || inherits(self, "InferenceAbstractKKQuantileRegrOneLik"))){
				private$y = .sanitize_proportion_response(private$y, interior = TRUE)
				private$cached_values$KKstats = NULL
				private$cached_values$beta_hat_T = NULL
				private$cached_values$s_beta_hat_T = NULL
				if (!is.null(private$compute_basic_match_data)) private$compute_basic_match_data()
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			if (!is.null(private[["compiled_cpp_stat_fn"]])) {
				cpp_fn = private$get_compiled_cpp_stat()
				arity = length(formals(cpp_fn))
				return(as.numeric(
					if (arity >= 3L) cpp_fn(private$y, as.integer(private$w), as.integer(private$dead))
					else cpp_fn(private$y, as.integer(private$w))
				)[1L])
			}
			if (is.null(private$custom_randomization_statistic_function)) self$compute_estimate(estimate_only = estimate_only)
			else private$custom_randomization_statistic_function()
		}
	)
)
