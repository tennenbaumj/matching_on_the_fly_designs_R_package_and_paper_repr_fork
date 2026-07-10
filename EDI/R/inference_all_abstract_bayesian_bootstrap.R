#' Bayesian Bootstrap-capable Inference
#'
#' Abstract class for Dirichlet-weight Bayesian bootstrap inference layered on
#' top of the existing nonparametric bootstrap infrastructure.
#'
#' @keywords internal
InferenceBayesianBootstrap = R6::R6Class("InferenceBayesianBootstrap",
	lock_objects = FALSE,
	inherit = InferenceNonParamBootstrap,
	public = list(
		#' @description Recomputes the treatment estimate under Bayesian-bootstrap
		#'   subject-, block-, cluster-, or matched-set weights.
		#'
		#' This is an abstract hook implemented by concrete inference families that
		#' support weighted re-estimation.
		#'
		#' @param subject_or_block_weights Numeric Bayesian-bootstrap weights at the
		#'   design's exchangeable resampling unit. For ordinary designs these are
		#'   subject-level weights. For blocking, clustering, or matching designs
		#'   these may instead be block-, cluster-, pair-, or matched-set-level
		#'   weights, depending on \code{weighting_unit_type}.
		#' @param estimate_only If \code{TRUE}, compute only the point estimate for
		#'   the weighted replicate.
		#'
		#' @return A numeric treatment-effect estimate for the weighted replicate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			stop(class(self)[1], " must implement weighted bootstrap estimation.")
		},
		#' @description Creates the Bayesian-bootstrap distribution of the treatment
		#'   estimate using Dirichlet weights.
		#'
		#' @param B Number of Bayesian-bootstrap replicates. The default is 501.
		#' @param show_progress A flag indicating whether a progress bar should be displayed.
		#' @param debug If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics including error messages, warning messages,
		#'   counts of each, and summary proportions for iterations with errors,
		#'   warnings, and illegal (non-finite) values. Runs serially. Default
		#'   \code{FALSE}.
		#' @param weighting_unit_type Optional Bayesian-bootstrap weighting-unit
		#'   scheme. Legal public values are:
		#'   \describe{
		#'     \item{\code{NULL}}{Use the design's default weighting-unit logic. For
		#'       ordinary non-blocking designs this is the usual subject-level
		#'       Bayesian bootstrap. For certain blocking designs, \code{NULL} maps
		#'       to the same behavior as \code{"within_blocks"}.}
		#'     \item{\code{"within_blocks"}}{Only legal for blocking-style designs
		#'       that support block-aware weighting:
		#'       \code{DesignFixedBlocking}, \code{DesignFixedOptimalBlocks},
		#'       \code{DesignSeqOneByOneSPBR}, and
		#'       \code{DesignFixedBlockedCluster}. Draws Dirichlet weights on
		#'       observational units within each observed block/stratum. For blocked
		#'       cluster designs this means cluster-within-stratum weights.}
		#'     \item{\code{"resample_blocks"}}{Only legal for the same
		#'       blocking-style designs as \code{"within_blocks"}. Draws Dirichlet
		#'       weights on whole observed blocks/strata rather than on units within
		#'       each block.}
		#'   }
		#'   Any non-\code{NULL} value is rejected for designs outside that blocking
		#'   family.
		#'
		#' @return When \code{debug = FALSE} (default), a numeric vector of length
		#'   \code{B} containing the Bayesian-bootstrap estimates. When
		#'   \code{debug = TRUE}, a list with: \code{values}, \code{errors} (list of
		#'   character vectors, one per iteration), \code{warnings} (list of
		#'   character vectors, one per iteration), \code{num_errors},
		#'   \code{num_warnings}, \code{prop_iterations_with_errors},
		#'   \code{prop_iterations_with_warnings}, and
		#'   \code{prop_illegal_values}.
		approximate_bayesian_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, weighting_unit_type = NULL){
			private$active_resampling_operation = "bayesian_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Bayesian bootstrap inference")
				private$assert_valid_bootstrap_type(weighting_unit_type)
				assertCount(B, positive = TRUE)
				assertFlag(debug)
			}
			cache_key = private$bayesian_bootstrap_cache_key(B = B, weighting_unit_type = weighting_unit_type)
			if (!isTRUE(debug) && !is.null(private$cached_values$bayes_boot_distr_cache[[cache_key]])) {
				return(private$cached_values$bayes_boot_distr_cache[[cache_key]])
			}
			inf_template = self$duplicate()
			if (!is.null(private$seed)) {
				had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				on.exit(
					if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv),
					add = TRUE
				)
				set.seed(private$seed)
			}
			draws = replicate(
				as.integer(B),
				private$bayesian_bootstrap_sample_weights(weighting_unit_type = weighting_unit_type),
				simplify = FALSE
			)
			run_one_iter = function(worker_inf) {
				draw = draws[[1L]]
				worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
				as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
					subject_or_block_weights = draw$subject_or_block_weights,
					estimate_only = TRUE
				))[1L]
			}
			if (isTRUE(debug)) {
				run_debug_iter = function(draw, worker_inf = NULL, worker_state = NULL) {
					iter_warns = character(0)
					iter_errs = character(0)
					iter_val = withCallingHandlers(
						tryCatch({
							if (!is.null(worker_state)) {
								private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
								private$compute_bayesian_bootstrap_worker_estimate(worker_state)
							} else {
								worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
								as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
									subject_or_block_weights = draw$subject_or_block_weights,
									estimate_only = TRUE
								))[1L]
							}
						}, error = function(e) { iter_errs <<- c(iter_errs, conditionMessage(e)); NA_real_ }),
						warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
					)
					list(val = as.numeric(iter_val)[1L], errors = iter_errs, warnings = iter_warns)
				}
				actual_debug_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
				chunk_n = max(1L, min(as.integer(actual_debug_cores), as.integer(B)))
				chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
				chunks = split(seq_len(B), chunk_id)
				run_debug_chunk = if (isTRUE(private$use_reusable_bootstrap_worker())) {
					function(idxs) {
						worker_state = private$create_bootstrap_worker_state()
						lapply(idxs, function(idx) run_debug_iter(draw = draws[[idx]], worker_state = worker_state))
					}
				} else {
					function(idxs) {
						lapply(idxs, function(idx) {
							worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
							run_debug_iter(draw = draws[[idx]], worker_inf = worker_inf)
						})
					}
				}
				debug_results = if (actual_debug_cores <= 1L) {
					run_debug_chunk(seq_len(B))
				} else {
					unlist(private$par_lapply(
						chunks,
						run_debug_chunk,
						n_cores = actual_debug_cores,
						budget = 1L,
						show_progress = show_progress
					), recursive = FALSE, use.names = FALSE)
				}
				debug_results = Filter(function(x) is.list(x) && !is.null(x$val), debug_results)
				values = vapply(debug_results, function(x) as.numeric(x$val)[1L], numeric(1))
				errors_list = lapply(debug_results, `[[`, "errors")
				warnings_list = lapply(debug_results, `[[`, "warnings")
				num_errors_vec = lengths(errors_list)
				num_warnings_vec = lengths(warnings_list)
				if (is.null(private$cached_values$bayes_boot_distr_cache)) private$cached_values$bayes_boot_distr_cache = list()
				private$cached_values$bayes_boot_distr_cache[[cache_key]] = values
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
			actual_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
			if (actual_cores > 1L) {
				do_warmup_iter = if (isTRUE(private$use_reusable_bootstrap_worker())) {
					function() {
						worker_state = private$create_bootstrap_worker_state()
						draw = draws[[1L]]
						private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
						tryCatch(private$compute_bayesian_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
					}
				} else {
					function() {
						worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
						draw = draws[[1L]]
						tryCatch({
							worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
							as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
								subject_or_block_weights = draw$subject_or_block_weights,
								estimate_only = TRUE
							))[1L]
						}, error = function(e) NA_real_)
					}
				}
				system.time(do_warmup_iter())
				t_boot_warmup = system.time(do_warmup_iter())[[3]]
				fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.5
				if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores)) actual_cores = 1L
			}
			boot_distr = if (isTRUE(private$use_reusable_bootstrap_worker())) {
				private$compute_bayesian_bootstrap_distribution_with_reused_workers(
					draws = draws,
					actual_cores = actual_cores,
					show_progress = show_progress
				)
			} else {
				unlist(private$par_lapply(
					seq_along(draws),
					function(idx) {
						worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
						draw = draws[[idx]]
						tryCatch({
							worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
							as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
								subject_or_block_weights = draw$subject_or_block_weights,
								estimate_only = TRUE
							))[1L]
						}, error = function(e) NA_real_)
					},
					n_cores = actual_cores,
					show_progress = show_progress,
					export_list = list(inf_template = inf_template, draws = draws)
				))
			}
			boot_distr = as.numeric(boot_distr)
			if (is.null(private$cached_values$bayes_boot_distr_cache)) private$cached_values$bayes_boot_distr_cache = list()
			private$cached_values$bayes_boot_distr_cache[[cache_key]] = boot_distr
			boot_distr
		},
		#' @description Computes a Bayesian-bootstrap-based two-sided p-value for
		#'   the treatment effect.
		#'
		#' @param delta Null hypothesis value. Default 0.
		#' @param B Number of Bayesian-bootstrap replicates. Default 501.
		#' @param type Type of Bayesian-bootstrap p-value. Supported values are
		#'   \code{"percentile"} (default), \code{"symmetric"}, \code{"wald"},
		#'   \code{"studentized"} / \code{"bootstrap-t"} (pivots by replicate SE from
		#'   \code{compute_estimate_with_bootstrap_weights(..., estimate_only = FALSE)}),
		#'   and \code{"bca"} (bias-corrected and accelerated via leave-one-unit-out
		#'   Bayesian jackknife).
		#' @param na.rm If \code{TRUE}, discard non-finite bootstrap replicates before
		#'   computing the p-value. Otherwise, any non-finite replicate returns
		#'   \code{NA}.
		#' @param show_progress A flag indicating whether a progress bar should be displayed.
		#' @param min_number_usable_samples Minimum number of finite Bayesian-bootstrap
		#'   replicates required after filtering. Default 5.
		#' @param weighting_unit_type Optional Bayesian-bootstrap weighting-unit
		#'   scheme. See
		#'   \code{\link{InferenceBayesianBootstrap$approximate_bayesian_bootstrap_distribution_beta_hat_T}()}.
		#'
		#' @return A numeric two-sided p-value, or \code{NA_real_} if too few usable
		#'   replicates remain or the estimate is non-finite.
		compute_bayesian_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = FALSE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertFlag(na.rm)
			}
			type = tolower(type %||% "percentile")
			if (should_run_asserts()) {
				assertChoice(type, c("percentile", "symmetric", "wald", "studentized", "bootstrap-t", "bca"))
			}
			est = as.numeric(self$compute_estimate())[1L]
			if (!is.finite(est)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate("bayesian_bootstrap_original_estimate_unavailable")
				return(NA_real_)
			}
			if (type == "bca" && private$mark_jackknife_nonestimable_if_block_unsupported(unit = "auto")) {
				return(NA_real_)
			}
			if (type %in% c("studentized", "bootstrap-t")) {
				boot_stats = private$approximate_bayesian_bootstrap_statistics_beta_hat_T(
					B = B, show_progress = show_progress, na.rm = isTRUE(na.rm),
					require_se = TRUE, weighting_unit_type = weighting_unit_type
				)
				se_hat = tryCatch(private$infer_original_se(), error = function(e) NA_real_)
				if (!is.finite(se_hat) || se_hat <= 0) {
					if (isTRUE(private$harden)) private$cache_nonestimable_se("bayesian_bootstrap_original_standard_error_unavailable")
					return(NA_real_)
				}
				t_boot = tryCatch(
					private$studentized_bootstrap_pivots(
						theta = boot_stats$theta,
						se = boot_stats$se,
						est = est,
						se_hat = se_hat,
						min_number_usable_samples = min_number_usable_samples,
						symmetric = TRUE
					),
					error = function(e) numeric(0)
				)
				if (length(t_boot) < as.integer(min_number_usable_samples)) {
					if (isTRUE(private$harden)) private$cache_nonestimable_se("bayesian_bootstrap_unstable_studentized_standard_errors")
					return(NA_real_)
				}
				t_obs = abs(est - delta) / se_hat
				return(min(1, max(1 / length(t_boot), mean(t_boot >= t_obs))))
			}
			boot_distr = self$approximate_bayesian_bootstrap_distribution_beta_hat_T(
				B = B,
				show_progress = show_progress,
				weighting_unit_type = weighting_unit_type
			)
			if (isTRUE(na.rm)) boot_distr = boot_distr[is.finite(boot_distr)]
			else if (any(!is.finite(boot_distr))) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate("bayesian_bootstrap_nonfinite_estimates")
				return(NA_real_)
			}
			if (length(boot_distr) < as.integer(min_number_usable_samples)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate("bayesian_bootstrap_too_few_finite_estimates")
				return(NA_real_)
			}
			if (type == "percentile") {
				boot_null = boot_distr - mean(boot_distr) + delta
				n_bs = length(boot_null)
				return(min(1, max(2 / n_bs, 2 * min(
					sum(boot_null >= est) / n_bs,
					sum(boot_null <= est) / n_bs
				))))
			}
			if (type == "symmetric") {
				D_obs = abs(est - delta)
				D_boot = abs(boot_distr - mean(boot_distr))
				n_bs = length(D_boot)
				return(min(1, max(1 / n_bs, mean(D_boot >= D_obs))))
			}
			if (type == "bca") {
				p_value = tryCatch(
					private$pval_bayesian_bca(boot_distr, est, delta, weighting_unit_type = weighting_unit_type),
					error = function(e) {
						if (isTRUE(private$harden)) private$cache_nonestimable_se("bayesian_bootstrap_bca_pvalue_unavailable")
						NA_real_
					}
				)
				if (!is.finite(p_value) && isTRUE(private$harden) && !isTRUE(self$is_nonestimable())) {
					private$cache_nonestimable_se("bayesian_bootstrap_bca_pvalue_unavailable")
				}
				return(p_value)
			}
			se_boot = stats::sd(boot_distr)
			if (!is.finite(se_boot) || se_boot <= 0) {
				if (isTRUE(private$harden)) private$cache_nonestimable_se("bayesian_bootstrap_standard_error_unavailable")
				return(NA_real_)
			}
			2 * stats::pnorm(-abs((est - delta) / se_boot))
		},
		#' @description Computes a Bayesian-bootstrap confidence interval for the
		#'   treatment effect.
		#'
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of Bayesian-bootstrap replicates. Default 501.
		#' @param type Type of Bayesian-bootstrap interval. Supported values are
		#'   \code{"percentile"} (default), \code{"basic"}, \code{"wald"},
		#'   \code{"studentized"} / \code{"bootstrap-t"} (pivots by replicate SE from
		#'   \code{compute_estimate_with_bootstrap_weights(..., estimate_only = FALSE)}),
		#'   and \code{"bca"} (bias-corrected and accelerated via leave-one-unit-out
		#'   Bayesian jackknife).
		#' @param na.rm If \code{TRUE}, discard non-finite bootstrap replicates before
		#'   constructing the interval.
		#' @param show_progress A flag indicating whether a progress bar should be displayed.
		#' @param min_number_usable_samples Minimum number of finite Bayesian-bootstrap
		#'   replicates required after filtering. Default 5.
		#' @param weighting_unit_type Optional Bayesian-bootstrap weighting-unit
		#'   scheme. See
		#'   \code{\link{InferenceBayesianBootstrap$approximate_bayesian_bootstrap_distribution_beta_hat_T}()}.
		#'
		#' @return A length-2 numeric confidence interval. Returns
		#'   \code{c(NA_real_, NA_real_)} when the estimate is non-finite or too few
		#'   usable replicates remain.
		compute_bayesian_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
			}
			type = tolower(type %||% "percentile")
			if (should_run_asserts()) {
				assertChoice(type, c("percentile", "basic", "wald", "studentized", "bootstrap-t", "bca"))
			}
			est = as.numeric(self$compute_estimate())[1L]
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			if (!is.finite(est)) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_original_estimate_unavailable", stage = "estimate"))
			}
			if (type == "bca" && private$mark_jackknife_nonestimable_if_block_unsupported(unit = "auto")) {
				return(private$missing_bootstrap_ci(alpha, "jackknife_block_size_gt_one_not_supported", stage = "se"))
			}
			if (type %in% c("studentized", "bootstrap-t")) {
				boot_stats = private$approximate_bayesian_bootstrap_statistics_beta_hat_T(
					B = B, show_progress = show_progress, na.rm = isTRUE(na.rm),
					require_se = TRUE, weighting_unit_type = weighting_unit_type
				)
				ok = is.finite(boot_stats$theta) & is.finite(boot_stats$se) & boot_stats$se > 0
				if (sum(ok) < as.integer(min_number_usable_samples)) {
					return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_too_few_finite_standard_errors", stage = "se"))
				}
				ci[] = tryCatch(
					private$ci_studentized(boot_stats, alpha, est, min_number_usable_samples = min_number_usable_samples),
					error = function(e) private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_standard_error_ci_unavailable", stage = "se")
				)
				if (length(ci) < 2L || !all(is.finite(ci[1:2]))) {
					return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_standard_error_ci_unavailable", stage = "se"))
				}
				return(ci)
			}
			boot_distr = self$approximate_bayesian_bootstrap_distribution_beta_hat_T(
				B = B,
				show_progress = show_progress,
				weighting_unit_type = weighting_unit_type
			)
			boot_distr = boot_distr[is.finite(boot_distr)]
			if (length(boot_distr) < as.integer(min_number_usable_samples)) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_too_few_finite_estimates", stage = "estimate"))
			}
			if (type == "wald") {
				se_boot = stats::sd(boot_distr)
				if (!is.finite(se_boot) || se_boot <= 0) {
					return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_standard_error_unavailable", stage = "se"))
				}
				z = stats::qnorm(1 - alpha / 2)
				ci[] = c(est - z * se_boot, est + z * se_boot)
				return(ci)
			}
			if (type == "bca") {
				ci[] = tryCatch(
					private$ci_bayesian_bca(boot_distr, alpha, est, weighting_unit_type = weighting_unit_type),
					error = function(e) private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_bca_ci_unavailable", stage = "se")
				)
				if (length(ci) < 2L || !all(is.finite(ci[1:2]))) {
					return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_bca_ci_unavailable", stage = "se"))
				}
				return(ci)
			}
			ci[] = private$ci_from_boot_distribution(boot_distr, alpha, type, est = est)
			if (length(ci) < 2L || !all(is.finite(ci[1:2]))) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_ci_unavailable", stage = "estimate"))
			}
			ci
		}
	),
	private = list(
		current_bayesian_bootstrap_context = NULL,
		current_bayesian_bootstrap_subject_or_block_weights = NULL,
		bayesian_bootstrap_cache_key = function(B, weighting_unit_type = NULL){
			paste0(as.integer(B), "::", weighting_unit_type %||% "default")
		},
		build_bayesian_bootstrap_context = function(weighting_unit_type = NULL){
			n = private$des_obj$get_n()
			design_obj = private$des_obj
			is_matching_design = is(design_obj, "DesignMatching") &&
				isTRUE(tryCatch(design_obj$is_matching_design(), error = function(e) FALSE))
			is_blocking_design = is(design_obj, "DesignBlocking") &&
				isTRUE(tryCatch(design_obj$is_blocking_design(), error = function(e) FALSE))
			if (is_matching_design) {
				private$des_obj_priv_int$ensure_matching_structure_computed()
				cluster_ids = as.integer(design_obj$get_matching_cluster_ids(private$m))
				row_to_unit = match(cluster_ids, unique(cluster_ids))
				unit_group_id = rep(1L, max(row_to_unit))
				return(list(
					row_to_unit = as.integer(row_to_unit),
					unit_group_id = as.integer(unit_group_id),
					n_units = length(unit_group_id)
				))
			}
			if (is_blocking_design) {
				block_ids = as.integer(design_obj$get_block_ids())
				if (identical(weighting_unit_type, "resample_blocks")) {
					row_to_unit = match(block_ids, unique(block_ids))
					unit_group_id = rep(1L, max(row_to_unit))
					return(list(
						row_to_unit = as.integer(row_to_unit),
						unit_group_id = as.integer(unit_group_id),
						n_units = length(unit_group_id)
					))
				}
				if (is(design_obj, "DesignFixedBlockedCluster")) {
					cluster_ids = as.character(private$des_obj_priv_int$Xraw[seq_len(n), ][[private$des_obj_priv_int$cluster_col]])
					unit_keys = paste(block_ids, cluster_ids, sep = "::")
					row_to_unit = match(unit_keys, unique(unit_keys))
					unit_group_id = block_ids[match(unique(unit_keys), unit_keys)]
					unit_group_id = match(unit_group_id, unique(unit_group_id))
					return(list(
						row_to_unit = as.integer(row_to_unit),
						unit_group_id = as.integer(unit_group_id),
						n_units = length(unit_group_id)
					))
				}
				row_to_unit = seq_len(n)
				unit_group_id = match(block_ids, unique(block_ids))
				return(list(
					row_to_unit = as.integer(row_to_unit),
					unit_group_id = as.integer(unit_group_id),
					n_units = n
				))
			}
			list(
				row_to_unit = seq_len(n),
				unit_group_id = rep(1L, n),
				n_units = n
			)
		},
		bayesian_bootstrap_sample_weights = function(weighting_unit_type = NULL){
			ctx = private$build_bayesian_bootstrap_context(weighting_unit_type = weighting_unit_type)
			subject_or_block_weights = numeric(ctx$n_units)
			for (group_id in unique(ctx$unit_group_id)) {
				idx = which(ctx$unit_group_id == group_id)
				draw = stats::rgamma(length(idx), shape = 1, rate = 1)
				subject_or_block_weights[idx] = draw / sum(draw) * length(idx)
			}
			list(
				subject_or_block_weights = as.numeric(subject_or_block_weights),
				context = ctx
			)
		},
		expand_subject_or_block_weights_to_row_weights = function(subject_or_block_weights){
			ctx = private$current_bayesian_bootstrap_context
			if (is.null(ctx)) {
				stop("No Bayesian-bootstrap context is installed on this inference object.", call. = FALSE)
			}
			if (should_run_asserts()) {
				assertNumeric(subject_or_block_weights, len = ctx$n_units, lower = 0, any.missing = FALSE)
			}
			as.numeric(subject_or_block_weights[ctx$row_to_unit])
		},
		load_bayesian_bootstrap_weights_into_worker = function(worker_state, draw){
			worker_priv = worker_state$worker$.__enclos_env__$private
			worker_priv$current_bayesian_bootstrap_subject_or_block_weights = as.numeric(draw$subject_or_block_weights)
			worker_priv$current_bayesian_bootstrap_context = draw$context
		},
		compute_bayesian_bootstrap_worker_estimate = function(worker_state){
			worker_priv = worker_state$worker$.__enclos_env__$private
			theta = as.numeric(worker_state$worker$compute_estimate_with_bootstrap_weights(
				subject_or_block_weights = worker_priv$current_bayesian_bootstrap_subject_or_block_weights,
				estimate_only = TRUE
			))[1L]
			if (is.function(worker_state$worker$is_nonestimable) &&
			    isTRUE(worker_state$worker$is_nonestimable("estimate"))) {
				return(NA_real_)
			}
			theta
		},
		approximate_bayesian_bootstrap_statistics_beta_hat_T = function(B = 501, show_progress = TRUE, na.rm = TRUE, require_se = FALSE, weighting_unit_type = NULL) {
			private$active_resampling_operation = "bayesian_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (!isTRUE(require_se)) {
				theta = self$approximate_bayesian_bootstrap_distribution_beta_hat_T(
					B = B, show_progress = show_progress, weighting_unit_type = weighting_unit_type
				)
				if (isTRUE(na.rm)) theta = theta[is.finite(theta)]
				return(list(theta = theta, se = rep(NA_real_, length(theta))))
			}
			stats_cache_key = paste0("stats::", private$bayesian_bootstrap_cache_key(B, weighting_unit_type))
			if (!is.null(private$cached_values$bayes_boot_stats_cache[[stats_cache_key]])) {
				cached = private$cached_values$bayes_boot_stats_cache[[stats_cache_key]]
				if (isTRUE(na.rm)) {
					ok = is.finite(cached$theta) & is.finite(cached$se) & cached$se > 0
					return(list(theta = cached$theta[ok], se = cached$se[ok]))
				}
				return(cached)
			}
			if (!is.null(private$seed)) {
				had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				on.exit(
					if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv),
					add = TRUE
				)
				set.seed(private$seed)
			}
			inf_work = self$duplicate(make_fork_cluster = FALSE)
			inf_work$num_cores = 1L
			stats_mat = matrix(NA_real_, nrow = B, ncol = 2L)
			pb = if (isTRUE(show_progress) && B > 1L) {
				pb_obj = utils::txtProgressBar(min = 0, max = B, style = 3)
				on.exit(try(close(pb_obj), silent = TRUE), add = TRUE)
				pb_obj
			} else NULL
			for (b in seq_len(B)) {
				draw = private$bayesian_bootstrap_sample_weights(weighting_unit_type = weighting_unit_type)
				tryCatch({
					inf_work$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
					theta_b = as.numeric(inf_work$compute_estimate_with_bootstrap_weights(
						subject_or_block_weights = draw$subject_or_block_weights,
						estimate_only = FALSE
					))[1L]
					se_b = as.numeric(inf_work$.__enclos_env__$private$cached_values$s_beta_hat_T)[1L]
					if (is.finite(theta_b)) stats_mat[b, 1L] = theta_b
					if (is.finite(se_b) && se_b > 0) stats_mat[b, 2L] = se_b
				}, error = function(e) NULL)
				if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
			}
			result = list(theta = stats_mat[, 1L], se = stats_mat[, 2L])
			if (is.null(private$cached_values$bayes_boot_stats_cache)) private$cached_values$bayes_boot_stats_cache = list()
			private$cached_values$bayes_boot_stats_cache[[stats_cache_key]] = result
			if (isTRUE(na.rm)) {
				ok = is.finite(result$theta) & is.finite(result$se) & result$se > 0
				return(list(theta = result$theta[ok], se = result$se[ok]))
			}
			result
		},
		approximate_bayesian_jackknife_distribution_beta_hat_T = function(weighting_unit_type = NULL) {
			private$active_resampling_operation = "bayesian_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (private$mark_jackknife_nonestimable_if_block_unsupported(unit = "auto")) {
				return(numeric(0))
			}
			jack_cache_key = weighting_unit_type %||% "default"
			if (!is.null(private$cached_values$bayes_jack_distr_cache[[jack_cache_key]])) {
				return(private$cached_values$bayes_jack_distr_cache[[jack_cache_key]])
			}
			ctx = private$build_bayesian_bootstrap_context(weighting_unit_type = weighting_unit_type)
			n_units = ctx$n_units
			if (n_units <= 1L) return(numeric(0))
			unit_group_id = ctx$unit_group_id
			group_sizes = tabulate(unit_group_id, nbins = max(unit_group_id))
			inf_work = self$duplicate(make_fork_cluster = FALSE)
			inf_work$num_cores = 1L
			inf_work$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx
			jack = vapply(seq_len(n_units), function(k) {
				g = unit_group_id[k]
				mg = group_sizes[g]
				w = rep(1, n_units)
				if (mg > 1L) {
					w[which(unit_group_id == g)] = mg / (mg - 1)
				}
				w[k] = 0
				tryCatch(
					as.numeric(inf_work$compute_estimate_with_bootstrap_weights(
						subject_or_block_weights = w,
						estimate_only = TRUE
					))[1L],
					error = function(e) NA_real_
				)
			}, numeric(1))
			if (is.null(private$cached_values$bayes_jack_distr_cache)) private$cached_values$bayes_jack_distr_cache = list()
			private$cached_values$bayes_jack_distr_cache[[jack_cache_key]] = jack
			jack
		},
		ci_bayesian_bca = function(boot_distr, alpha, est, weighting_unit_type = NULL) {
			jack = private$approximate_bayesian_jackknife_distribution_beta_hat_T(weighting_unit_type = weighting_unit_type)
			jack = jack[is.finite(jack)]
			if (length(jack) < 2L) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_bca_jackknife_unavailable", stage = "se"))
			}
			p_less = mean(boot_distr < est)
			p_less = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_less))
			z0 = stats::qnorm(p_less)
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			if (!is.finite(z0) || !is.finite(a) || abs(z0) > 2.5 || abs(a) > 1) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_bca_unstable_bias_or_acceleration", stage = "se"))
			}
			z_alpha = stats::qnorm(c(alpha / 2, 1 - alpha / 2))
			denom = 1 - a * (z0 + z_alpha)
			if (any(!is.finite(denom)) || any(abs(denom) < sqrt(.Machine$double.eps))) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_bca_unstable_bias_or_acceleration", stage = "se"))
			}
			adj = stats::pnorm(z0 + (z0 + z_alpha) / denom)
			prob_eps = 1 / (length(boot_distr) + 1)
			adj = pmin(1 - prob_eps, pmax(prob_eps, adj))
			adj = sort(adj)
			if (any(adj <= 2 * prob_eps) || any(adj >= 1 - 2 * prob_eps)) {
				return(private$missing_bootstrap_ci(alpha, "bayesian_bootstrap_bca_adjustment_on_boundary", stage = "se"))
			}
			if (diff(adj) < prob_eps) {
				return(stats::quantile(boot_distr, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8))
			}
			stats::quantile(boot_distr, probs = adj, names = FALSE, type = 8)
		},
		pval_bayesian_bca = function(boot_distr, est, delta, weighting_unit_type = NULL) {
			jack = private$approximate_bayesian_jackknife_distribution_beta_hat_T(weighting_unit_type = weighting_unit_type)
			jack = jack[is.finite(jack)]
			if (length(jack) < 2L) {
				private$cache_nonestimable_se("bayesian_bootstrap_bca_jackknife_unavailable")
				return(NA_real_)
			}
			p_less = mean(boot_distr < est)
			p_less = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_less))
			z0 = stats::qnorm(p_less)
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			if (!is.finite(z0) || !is.finite(a) || abs(z0) > 2.5 || abs(a) > 1) {
				private$cache_nonestimable_se("bayesian_bootstrap_bca_unstable_bias_or_acceleration")
				return(NA_real_)
			}
			z_alpha = stats::qnorm(c(0.025, 0.975))
			denom_ci = 1 - a * (z0 + z_alpha)
			if (any(!is.finite(denom_ci)) || any(abs(denom_ci) < sqrt(.Machine$double.eps))) {
				private$cache_nonestimable_se("bayesian_bootstrap_bca_unstable_bias_or_acceleration")
				return(NA_real_)
			}
			prob_eps = 1 / (length(boot_distr) + 1)
			adj_ci = sort(stats::pnorm(z0 + (z0 + z_alpha) / denom_ci))
			if (any(!is.finite(adj_ci)) || any(adj_ci <= 2 * prob_eps) || any(adj_ci >= 1 - 2 * prob_eps)) {
				private$cache_nonestimable_se("bayesian_bootstrap_bca_adjustment_on_boundary")
				return(NA_real_)
			}
			p_delta = mean(boot_distr < delta)
			p_delta = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_delta))
			z_delta = stats::qnorm(p_delta)
			s = z_delta - z0
			denom = 1 + a * s
			if (!is.finite(denom) || abs(denom) < sqrt(.Machine$double.eps)) {
				private$cache_nonestimable_se("bayesian_bootstrap_bca_unstable_bias_or_acceleration")
				return(NA_real_)
			}
			adj_z = s / denom - z0
			if (!is.finite(adj_z) || abs(adj_z) > 8) {
				private$cache_nonestimable_se("bayesian_bootstrap_bca_adjustment_on_boundary")
				return(NA_real_)
			}
			p_raw = min(1, 2 * min(stats::pnorm(adj_z), 1 - stats::pnorm(adj_z)))
			min(1, max(1 / length(boot_distr), p_raw))
		},
		compute_bayesian_bootstrap_distribution_with_reused_workers = function(draws, actual_cores, show_progress = FALSE){
			B = length(draws)
			chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
			chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
			chunks = split(seq_len(B), chunk_id)
			run_chunk = function(idxs) {
				worker_state = private$create_bootstrap_worker_state()
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					draw = draws[[idxs[[k]]]]
					out[k] = tryCatch({
						private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
						private$compute_bayesian_bootstrap_worker_estimate(worker_state)
					}, error = function(e) NA_real_)
				}
				out
			}
			if (actual_cores <= 1L) {
				return(as.numeric(run_chunk(seq_len(B))))
			}
			as.numeric(unlist(private$par_lapply(
				chunks,
				run_chunk,
				n_cores = actual_cores,
				budget = 1L,
				show_progress = show_progress
			), use.names = FALSE))
		}
	)
)
