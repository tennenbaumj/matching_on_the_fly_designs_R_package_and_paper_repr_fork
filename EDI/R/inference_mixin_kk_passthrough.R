#' Mixin for KK Matching-on-the-Fly Pass-Through Inference
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' providing KK design validation, bootstrap distribution logic, and match-data
#' helpers. Splice into a daughter class via
#' \code{public = c(InferenceMixinKKPassThrough$public, list(...))} and
#' \code{private = c(InferenceMixinKKPassThrough$private, list(...))}.
#'
#' Capability flag: \code{private$kk_passthrough == TRUE}.
#'
#' @keywords internal
InferenceMixinKKPassThrough = list(
	public = list(
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect
		#'
		#' @param B  					Number of bootstrap samples. The default is 501.
		#'
		#' @return A vector of length \code{B} with the bootstrap values of the estimates of the
		#'   treatment effect
		#'
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics. Default \code{FALSE}.
		#' @param bootstrap_type Optional bootstrap-resampling scheme. Legal public values are:
		#'   \describe{
		#'     \item{\code{NULL}}{Use the design's default row-resampling bootstrap.}
		#'     \item{\code{"within_blocks"}}{Only legal for blocking-style designs.}
		#'     \item{\code{"resample_blocks"}}{Only legal for the same blocking-style designs.}
		#'   }
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
				if (should_run_asserts()) {
					private$assert_valid_bootstrap_type(bootstrap_type)
				}
				if (!private$has_match_structure){
					super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug = debug, bootstrap_type = bootstrap_type)
				} else {
					if (should_run_asserts()) {
						assertCount(B, positive = TRUE)
					}
					if (is.null(private$cached_values$KKstats)){
						private$compute_basic_match_data()
					}
					n = private$n
					y = private$y
					dead = private$dead
					w = private$w
					X = private$get_X()
					# Let Design initialise and own the bootstrap pair structure
					des_priv = private$des_obj_priv_int
					des_priv$init_matching_bootstrap_structure()
					n_reservoir = des_priv$boot_n_reservoir
					m            = nrow(des_priv$boot_pair_rows)
					# For the C++ fast-path we still need i_reservoir / m_vec in Inference scope
					i_reservoir  = des_priv$boot_i_reservoir
					m_vec = private$m
					if (is.null(m_vec)) m_vec = rep(0L, n)
					m_vec = as.integer(m_vec)
					m_vec[is.na(m_vec)] = 0L
					# Check if subclass provides a C++ OpenMP dispatcher to bypass the slow R loop
					if (!isTRUE(debug) && private$has_private_method("compute_fast_bootstrap_distr")) {
						fast_distr = private$compute_fast_bootstrap_distr(B, i_reservoir, n_reservoir, m, y, w, m_vec)
						if (!is.null(fast_distr)) {
							return(fast_distr)
						}
					}
					kk_boot_context = private$create_kk_bootstrap_context(
						y = y, dead = dead, w = w, X = X,
						m = m, n_reservoir = n_reservoir
					)
					kk_boot_draws = replicate(
						as.integer(B),
						des_priv$draw_bootstrap_indices(),
						simplify = FALSE
					)
					if (isTRUE(private$use_reusable_kk_bootstrap_worker())) {
						if (isTRUE(debug)) {
							return(private$compute_kk_bootstrap_debug_with_reused_worker(kk_boot_draws, kk_boot_context))
						}
						actual_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
						if (actual_cores > 1L) {
							do_warmup_iter = function() {
								worker_state = private$create_kk_bootstrap_worker_state(kk_boot_context)
								sample_info = kk_boot_draws[[1L]]
								private$load_kk_bootstrap_sample_into_worker(worker_state, sample_info)
								tryCatch(private$compute_kk_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
							}
							system.time(do_warmup_iter())
							t_boot_warmup = system.time(do_warmup_iter())[[3]]
							fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.5
							if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores)) {
								actual_cores = 1L
							}
						}
						return(private$compute_kk_bootstrap_distribution_with_reused_workers(
							kk_boot_draws = kk_boot_draws,
							kk_boot_context = kk_boot_context,
							actual_cores = actual_cores,
							show_progress = show_progress
						))
					}
					if (isTRUE(debug)) {
						debug_results = vector("list", B)
						has_res_stat_debug = private$has_private_method("compute_reservoir_and_match_statistics")
						for (b in seq_len(B)) {
							boot_sample = kk_boot_draws[[b]]
							i_b   = boot_sample$i_b
							m_vec_b = boot_sample$m_vec_b
							iter_warns = character(0)
							iter_val = withCallingHandlers(
								tryCatch({
									boot_inf_obj = self$duplicate()
									boot_inf_obj$.__enclos_env__$private$y = y[i_b]
									boot_inf_obj$.__enclos_env__$private$y_temp = boot_inf_obj$.__enclos_env__$private$y
									boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
										boot_inf_obj$.__enclos_env__$private$w = w[i_b]
										boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
										private$clear_kk_bootstrap_worker_design_caches(boot_inf_obj$.__enclos_env__$private)
										boot_inf_obj$.__enclos_env__$private$m = m_vec_b
										boot_inf_obj$.__enclos_env__$private$cached_values = list(
										KKstats = compute_bootstrap_matching_stats_cpp(
											X = X, y = y, w = w, i_b = i_b, n_reservoir = n_reservoir
										)
									)
									if (has_res_stat_debug) {
										boot_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
									}
									as.numeric(boot_inf_obj$compute_estimate(estimate_only = TRUE))[1L]
								}, error = function(e) list(val = NA_real_, error = conditionMessage(e))),
								warning = function(w) {
									iter_warns <<- c(iter_warns, conditionMessage(w))
									invokeRestart("muffleWarning")
								}
							)
							res_val = if (is.list(iter_val)) iter_val$val else iter_val
							res_err = if (is.list(iter_val)) iter_val$error else NULL
							debug_results[[b]] = list(
								val = res_val,
								error = res_err,
								warnings = iter_warns
							)
						}
						return(debug_results)
					}
					# Standard Loop
					cores_to_use = private$effective_parallel_cores("bootstrap", self$num_cores)
					res = if (cores_to_use > 1L) {
						unlist(private$par_lapply(
							as.list(seq_len(B)),
							function(b) {
								tryCatch({
									boot_sample = kk_boot_draws[[b]]
									i_b = boot_sample$i_b
									m_vec_b = boot_sample$m_vec_b
									boot_inf_obj = self$duplicate()
									boot_inf_obj$.__enclos_env__$private$y = y[i_b]
									boot_inf_obj$.__enclos_env__$private$y_temp = boot_inf_obj$.__enclos_env__$private$y
									boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
										boot_inf_obj$.__enclos_env__$private$w = w[i_b]
										boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
										private$clear_kk_bootstrap_worker_design_caches(boot_inf_obj$.__enclos_env__$private)
										boot_inf_obj$.__enclos_env__$private$m = m_vec_b
										boot_inf_obj$.__enclos_env__$private$cached_values = list(
										KKstats = compute_bootstrap_matching_stats_cpp(
											X = X, y = y, w = w, i_b = i_b, n_reservoir = n_reservoir
										)
									)
									as.numeric(boot_inf_obj$compute_estimate(estimate_only = TRUE))[1L]
								}, error = function(e) { NA_real_ })
							},
							n_cores = cores_to_use,
							show_progress = show_progress
						))
					} else {
						vapply(seq_len(B), function(b) {
							tryCatch({
								boot_sample = kk_boot_draws[[b]]
								i_b = boot_sample$i_b
								m_vec_b = boot_sample$m_vec_b
								boot_inf_obj = self$duplicate()
								boot_inf_obj$.__enclos_env__$private$y = y[i_b]
								boot_inf_obj$.__enclos_env__$private$y_temp = boot_inf_obj$.__enclos_env__$private$y
								boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
									boot_inf_obj$.__enclos_env__$private$w = w[i_b]
									boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
									private$clear_kk_bootstrap_worker_design_caches(boot_inf_obj$.__enclos_env__$private)
									boot_inf_obj$.__enclos_env__$private$m = m_vec_b
									boot_inf_obj$.__enclos_env__$private$cached_values = list(
									KKstats = compute_bootstrap_matching_stats_cpp(
										X = X, y = y, w = w, i_b = i_b, n_reservoir = n_reservoir
									)
								)
								as.numeric(boot_inf_obj$compute_estimate(estimate_only = TRUE))[1L]
							}, error = function(e) { NA_real_ })
						}, numeric(1L))
					}
					res
				}
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE) {
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			
			# If the class provides its own specialized weighted estimator (e.g. via backend), use it.
			if (is.function(private$compute_weighted_estimate_ivwc)) {
				w_info = kk_pair_and_reservoir_bootstrap_weights(private, row_weights)
				return(private$compute_weighted_estimate_ivwc(w_info, estimate_only))
			}
			
			# Fallback for KK designs: weighted combination of matches and reservoir.
			if (is.null(private$cached_values$KKstats)) {
				private$compute_basic_match_data()
			}
			stats = private$cached_values$KKstats
			if (is.null(stats)) return(NA_real_)
			
			w_info = kk_pair_and_reservoir_bootstrap_weights(private, row_weights)
			
			# 1. Matched differences mean
			d_bar_w = if (length(w_info$pair_weights) > 0) {
				sum(stats$y_matched_diffs * w_info$pair_weights, na.rm = TRUE) / sum(w_info$pair_weights, na.rm = TRUE)
			} else NA_real_
			
			# 2. Reservoir mean diff
			r_idx = w_info$reservoir_idx
			r_w = w_info$reservoir_weights
			if (length(r_idx) > 0 && !is.null(stats$y_reservoir)) {
				y_r = stats$y_reservoir; w_r = stats$w_reservoir
				num_t = sum(y_r[w_r == 1] * r_w[w_r == 1], na.rm = TRUE)
				den_t = sum(r_w[w_r == 1], na.rm = TRUE)
				num_c = sum(y_r[w_r == 0] * r_w[w_r == 0], na.rm = TRUE)
				den_c = sum(r_w[w_r == 0], na.rm = TRUE)
				r_bar_t = if (is.finite(den_t) && den_t > 0) num_t / den_t else NA_real_
				r_bar_c = if (is.finite(den_c) && den_c > 0) num_c / den_c else NA_real_
				r_bar_w = r_bar_t - r_bar_c
			} else r_bar_w = NA_real_
			
			# Combine using fixed weights from observed fit
			w_star = stats$w_star
			if (is.null(w_star) || is.na(w_star)) {
				if (!is.na(d_bar_w)) return(d_bar_w)
				return(r_bar_w)
			}
			if (is.na(d_bar_w)) return(r_bar_w)
			if (is.na(r_bar_w)) return(d_bar_w)
			
			as.numeric(w_star * d_bar_w + (1 - w_star) * r_bar_w)
		}
	),
	private = list(
		m = NULL,
		kk_passthrough = TRUE,
		y_temp = NULL,
		dead = NULL,
		w = NULL,
		X = NULL,
		any_censoring = NULL,
		best_par = NULL,
		optimization_alg = "lbfgs",
		cached_mod = NULL,
		best_X_colnames = NULL,
		best_Xmm_colnames = NULL,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_information_preference = function(){
			FALSE
		},
		supports_observed_information = function(){
			FALSE
		},
		get_supported_testing_types_impl = function(){
			"wald"
		},
		get_supported_information_preferences_impl = function(){
			"auto"
		},
		use_reusable_kk_bootstrap_worker = function(){
			FALSE
		},
		init_kk_passthrough = function(des_obj){
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")) {
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14) or DesignFixedBinaryMatch.")
				}
			}
			if (private$has_match_structure){
				if (inherits(des_obj, "DesignFixedBinaryMatch")){
					des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
				}
				private$m = des_obj$.__enclos_env__$private$m
				private$compute_basic_match_data()
			}
		},
		create_kk_bootstrap_context = function(y, dead, w, X, m, n_reservoir){
			X_mat = if (is.null(X)) {
				matrix(numeric(0), nrow = length(y), ncol = 0L)
			} else {
				as.matrix(X)
			}
			list(
				y = as.numeric(y),
				dead = dead,
				w = as.integer(w),
				X = X_mat,
				m = as.integer(m),
				n_reservoir = as.integer(n_reservoir)
			)
		},
		create_kk_bootstrap_worker_state = function(kk_boot_context){
			worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			worker$num_cores = 1L
			worker_priv = worker$.__enclos_env__$private
			list(
				worker = worker,
				worker_priv = worker_priv,
				base_y = kk_boot_context$y,
				base_dead = kk_boot_context$dead,
				base_w = kk_boot_context$w,
				base_X = kk_boot_context$X,
				n_reservoir = kk_boot_context$n_reservoir,
				has_res_stat = private$object_has_private_method(worker, "compute_reservoir_and_match_statistics")
			)
		},
			load_kk_bootstrap_sample_into_worker = function(worker_state, sample_info){
				worker_priv = worker_state$worker_priv
				i_b = sample_info$i_b
			worker_priv$y = worker_state$base_y[i_b]
			worker_priv$y_temp = worker_priv$y
				worker_priv$dead = worker_state$base_dead[i_b]
				worker_priv$w = worker_state$base_w[i_b]
				worker_priv$X = worker_state$base_X[i_b, , drop = FALSE]
				private$clear_kk_bootstrap_worker_design_caches(worker_priv)
				worker_priv$cached_values = list(
				KKstats = compute_bootstrap_matching_stats_cpp(
					X = worker_state$base_X,
					y = worker_state$base_y,
					w = worker_state$base_w,
					i_b = i_b,
					n_reservoir = worker_state$n_reservoir
				)
			)
			worker_priv$best_X_colnames = NULL
			worker_priv$best_Xmm_colnames = NULL
			worker_priv$fit_warm_coefficients = NULL
			worker_priv$cached_mod = NULL
			worker_priv$m = sample_info$m_vec_b
				if (isTRUE(worker_state$has_res_stat)) {
					worker_priv$compute_reservoir_and_match_statistics()
				}
			},
			clear_kk_bootstrap_worker_design_caches = function(worker_priv){
				worker_priv$cached_design_matrix = NULL
				worker_priv$cached_w_for_design_matrix = NULL
				worker_priv$cached_harden_for_design_matrix = NULL
				worker_priv$cached_hardened_X_cov = NULL
				worker_priv$cached_reduced_X = NULL
				worker_priv$cached_X_full_for_reduced = NULL
				worker_priv$cached_keep_for_reduced = NULL
				worker_priv$cached_j_treat_for_reduced = NULL
				worker_priv$reduced_design_keep_cache = NULL
				worker_priv$fixed_covariate_keep_cache = NULL
				invisible(NULL)
			},
			compute_kk_bootstrap_worker_estimate = function(worker_state){
				as.numeric(worker_state$worker$compute_estimate(estimate_only = TRUE))[1L]
			},
		compute_kk_bootstrap_debug_with_reused_worker = function(kk_boot_draws, kk_boot_context){
			B = length(kk_boot_draws)
			worker_state = private$create_kk_bootstrap_worker_state(kk_boot_context)
			debug_results = vector("list", B)
			for (b in seq_len(B)) {
				sample_info = kk_boot_draws[[b]]
				iter_warns = character(0)
				iter_val = withCallingHandlers(
					tryCatch({
						private$load_kk_bootstrap_sample_into_worker(worker_state, sample_info)
						private$compute_kk_bootstrap_worker_estimate(worker_state)
					}, error = function(e) list(val = NA_real_, error = conditionMessage(e))),
					warning = function(w) {
						iter_warns <<- c(iter_warns, conditionMessage(w))
						invokeRestart("muffleWarning")
					}
				)
				res_val = if (is.list(iter_val)) iter_val$val else iter_val
				res_err = if (is.list(iter_val)) iter_val$error else NULL
				debug_results[[b]] = list(
					val = res_val,
					error = res_err,
					warnings = iter_warns
				)
			}
			debug_results
		},
		compute_kk_bootstrap_distribution_with_reused_workers = function(kk_boot_draws, kk_boot_context, actual_cores, show_progress){
			B = length(kk_boot_draws)
			if (actual_cores > 1L) {
				unlist(private$par_lapply(
					as.list(seq_len(B)),
					function(idx) {
						b = as.integer(idx)[1L]
						# In a fork, we can't easily reuse state across indices without complicated orchestration.
						# Re-creating a fresh state is safer and the worker state setup is cheap.
						worker_state = private$create_kk_bootstrap_worker_state(kk_boot_context)
						private$load_kk_bootstrap_sample_into_worker(worker_state, kk_boot_draws[[b]])
						tryCatch(private$compute_kk_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
					},
					n_cores = actual_cores,
					budget = 1L,
					show_progress = show_progress
				), recursive = FALSE, use.names = FALSE)
			} else {
				worker_state = private$create_kk_bootstrap_worker_state(kk_boot_context)
				pbar = if (show_progress) utils::txtProgressBar(min = 0, max = B, style = 3) else NULL
				res = vapply(seq_len(B), function(b) {
					private$load_kk_bootstrap_sample_into_worker(worker_state, kk_boot_draws[[b]])
					val = tryCatch(private$compute_kk_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
					if (show_progress) utils::setTxtProgressBar(pbar, b)
					val
				}, numeric(1L))
				if (show_progress) close(pbar)
				res
			}
		},
		compute_basic_kk_match_data_impl = function(){
			if (!isTRUE(private$has_match_structure)) {
				private$cache_nonestimable_estimate("kk_design_required")
				return(invisible(NULL))
			}
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv     = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
			private$cached_values$KKstats
		}
	)
)
