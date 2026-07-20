#' Bootstrap-based Inference
#'
#' Abstract class for bootstrap-based inference.
#'
#' @section Design-specific validity caveats for the nonparametric bootstrap:
#' The nonparametric bootstrap resamples experimental units with replacement from their
#' empirical distribution, carrying each unit's realized \code{(x, w, y)} into the
#' replicate, and recomputes the estimator. Its validity rests on the resampled units
#' being (approximately) iid draws from the design's unit-level superpopulation. Since
#' covariate-adaptive designs induce dependence among the assignments \eqn{w_i} (and
#' between \eqn{w} and \eqn{X}), the appropriate resampling unit and the fidelity with
#' which the design's dependence is replicated differ by design. In all cases below the
#' inference is asymptotic, never finite-sample exact (for exact finite-sample inference
#' under the design's actual randomization mechanism, use the randomization tests and
#' randomization confidence intervals instead). Where a balance constraint of the design
#' is broken by resampling, the bootstrap variance is inflated relative to the design's
#' true sampling variance, so tests and intervals err \emph{conservative}
#' (over-coverage), not anti-conservative.
#'
#' \describe{
#'   \item{\code{DesignFixedBernoulli}, \code{DesignSeqOneByOneBernoulli}}{Assignments
#'     are iid coin flips independent of \eqn{X}, so rows genuinely are iid and
#'     row-level resampling is fully justified. No caveat.}
#'   \item{\code{DesignFixediBCRD}, \code{DesignSeqOneByOneiBCRD}}{Assignment depends
#'     only on the treatment counts (completely randomized / without-replacement urn),
#'     inducing negative correlation among the \eqn{w_i} through the fixed-margin
#'     constraint. Row-level iid resampling does not replicate this constraint:
#'     replicates have a random number of treated subjects. The extra variability is
#'     \eqn{O(1/n)}, so the bootstrap is conservative by an asymptotically negligible
#'     amount.}
#'   \item{\code{DesignSeqOneByOneEfron}, \code{DesignSeqOneByOneUrn}}{Assignment
#'     depends on the running treatment imbalance (not on \eqn{X}), inducing serial
#'     negative dependence among the \eqn{w_i}. Row-level resampling ignores this
#'     dependence; as with fixed margins the effect on smooth estimators is
#'     \eqn{O(1/n)}, so the bootstrap is conservative by a negligible amount.}
#'   \item{\code{DesignSeqOneByOneRandomBlockSize} (no strata)}{Permuted-block balance
#'     over entry order is broken by row-level resampling. Conservative, minor: only
#'     the counts constraint is lost since the design does not use \eqn{X}.}
#'   \item{\code{DesignSeqOneByOneRandomBlockSize} (with strata),
#'     \code{DesignSeqOneByOneSPBR}}{Resampling is within-strata, preserving stratum
#'     sizes and the stratum-covariate composition. The within-block time-order balance
#'     inside each stratum is still broken, so replicates have random within-stratum
#'     treatment counts. Conservative, minor.}
#'   \item{\code{DesignFixedBlocking}}{Resampling is within-strata by default
#'     (\code{bootstrap_type = "within_blocks"}), preserving stratum sizes; the exact
#'     within-stratum treatment/control split is not enforced in replicates, so the
#'     block-randomization variance reduction is partially unreplicated. Conservative,
#'     minor. \code{bootstrap_type = "resample_blocks"} instead resamples whole blocks,
#'     preserving within-block composition at the price of fewer resampling atoms.}
#'   \item{\code{DesignFixedOptimalBlocks}}{Same within-block resampling caveats as
#'     \code{DesignFixedBlocking}, plus the blocks themselves are computed from the
#'     realized covariate sample: the block structure is a global function of the data
#'     that the bootstrap conditions on rather than re-derives. The justification for
#'     this conditioning is asymptotic: as \eqn{n} grows the blocking depends on the
#'     sample only through the (convergent) empirical distribution of \eqn{X}, so
#'     between-block dependence vanishes. Conservative.}
#'   \item{\code{DesignFixedCluster}}{Assignment is at the cluster level and outcomes
#'     are correlated within clusters, so whole clusters are resampled with
#'     replacement. This is the correct exchangeable unit; with few clusters the
#'     bootstrap distribution rests on few resampling atoms and becomes unstable.
#'     Asymptotics are in the number of clusters, not the number of subjects.}
#'   \item{\code{DesignFixedBlockedCluster}}{Clusters are resampled within strata,
#'     matching both levels of the design's dependence (stratum and cluster). Sound,
#'     with the same small-sample caution: few clusters per stratum means few
#'     resampling atoms per stratum, and asymptotics are in the number of clusters.}
#'   \item{\code{DesignSeqOneByOnePocockSimon}}{Minimization makes each assignment a
#'     near-deterministic function of the running stratum-count imbalances. Row-level
#'     iid resampling does not replicate this balance-forcing, so the bootstrap
#'     variance corresponds to iid assignment rather than the (smaller)
#'     minimization-design variance (cf. Bugni, Canay & Shah 2018). Conservative, with
#'     the largest expected over-coverage among the sequential designs.}
#'   \item{\code{DesignSeqOneByOneAtkinson}}{The biased-coin \eqn{D_A}-optimal rule
#'     makes \eqn{w_i} depend on the full covariate and assignment history, and
#'     conditional assignment probabilities differ from 1/2. Row-level resampling does
#'     not replicate the covariate balance the rule enforces. Conservative, moderate.}
#'   \item{\code{DesignFixedAOptimal}, \code{DesignFixedDOptimal},
#'     \code{DesignFixedGreedy}, \code{DesignFixedRerandomization}}{The observed
#'     \eqn{w} vector is one draw from a tightly constrained (optimized or
#'     acceptance-sampled) set of allocations. Resampled replicates carry per-row
#'     assignments whose recombined \eqn{w} vector no longer satisfies the balance
#'     constraint, so the bootstrap reflects the variance of unconstrained assignment
#'     (cf. Li, Ding & Rubin 2018 for rerandomization). Conservative,
#'     moderate-to-large: the stronger the optimization, the greater the
#'     over-coverage.}
#'   \item{\code{DesignFixedMatchingGreedyPairSwitching}}{The greedy switching search
#'     only ever flips assignments \emph{within} binary-match pairs, so every pair has
#'     exactly one treated subject; the bootstrap resamples intact pairs, preserving
#'     the within-pair anticorrelation. Remaining caveats: the pairing is a global
#'     function of the sample (conditioned on, justified asymptotically as for the
#'     matched designs below), and the greedy choice of \emph{which} pair member is
#'     treated couples the pairs, which resampling does not replicate --- the residual
#'     effect errs conservative.}
#'   \item{\code{DesignFixedBinaryMatch}, \code{DesignFixedNaiveMatch}}{Matched pairs
#'     are resampled intact, preserving the within-pair anticorrelation of \eqn{w} and
#'     the pair-level variance reduction. The pairing itself is a global function of
#'     the covariate sample (an Abadie & Imbens 2008-type concern): pairs are
#'     exchangeable but not exactly independent. Validity is asymptotic --- as \eqn{n}
#'     grows the pairing depends on the sample only through the empirical distribution
#'     of \eqn{X} and between-pair dependence vanishes --- and the bootstrap conditions
#'     on the realized match structure.}
#'   \item{\code{DesignSeqOneByOneKK14}}{Matched pairs and reservoir subjects are
#'     resampled separately as intact units, preserving within-pair anticorrelation
#'     and the reservoir's Bernoulli assignments. Same asymptotic caveats as the fixed
#'     matched designs, plus the split between number of pairs and reservoir size is
#'     treated as fixed rather than random. Sound asymptotically.}
#'   \item{\code{DesignSeqOneByOneKK21}, \code{DesignSeqOneByOneKK21stepwise}}{All
#'     \code{DesignSeqOneByOneKK14} caveats apply, plus the matching weights are
#'     estimated from earlier \emph{responses}, so \eqn{W} depends on \eqn{y} as well
#'     as \eqn{X}. The bootstrap conditions on the realized response-adaptive weights
#'     and match structure rather than replicating their sampling variability; this
#'     extra conditioning is not quantified, and validity remains asymptotic.}
#'   \item{\code{DesignFixedFactorial}}{Row-level resampling does not replicate the
#'     balanced allocation across factor combinations. Conservative, minor.}
#'   \item{\code{DesignFixedCustom}, \code{DesignCustomSequential}}{Warning: iid
#'     row-level resampling is used because the package has no knowledge of the
#'     user-supplied assignment mechanism. If that mechanism balances on covariates,
#'     the bootstrap is likely conservative; if it induces clustering or other
#'     positive dependence, the bootstrap may not even be valid (anti-conservative).
#'     Use the randomization-based inference, which draws from the actual custom
#'     mechanism, whenever possible.}
#' }
#'
#' @keywords internal
InferenceNonParamBootstrap = R6::R6Class("InferenceNonParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceRandCI,
	public = list(
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#'   The resampling unit is design-specific (rows, within-strata rows, matched pairs plus
		#'   reservoir, or clusters); see the class-level section \emph{Design-specific validity
		#'   caveats for the nonparametric bootstrap} for the conservativeness and asymptotics of
		#'   each concrete design.
		#'
		#' @param B  					Number of bootstrap samples. The default is 501.
		#' @param show_progress  		A flag indicating whether a progress bar should be displayed.
		#' @param debug  				If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics including error messages, warning messages, counts of each,
		#'   and summary proportions for iterations with errors, warnings, and illegal (non-finite)
		#'   values. Runs serially. Default \code{FALSE}.
		#' @return 	When \code{debug = FALSE} (default), a numeric vector of length \code{B}
		#'   containing the bootstrap estimates. When \code{debug = TRUE}, a list with: \code{values},
		#'   \code{errors} (list of character vectors, one per iteration), \code{warnings} (list of
		#'   character vectors, one per iteration), \code{num_errors}, \code{num_warnings},
		#'   \code{prop_iterations_with_errors}, \code{prop_iterations_with_warnings}, and
		#'   \code{prop_illegal_values}.
		#' @param bootstrap_type Optional bootstrap-resampling scheme. Legal public values are:
		#'   \describe{
		#'     \item{\code{NULL}}{Use the design's default row-resampling bootstrap. For ordinary
		#'       non-blocking designs this is the usual subject-level resample-with-replacement
		#'       bootstrap. For certain blocking designs, \code{NULL} maps to the same behavior
		#'       as \code{"within_blocks"}.}
		#'     \item{\code{"within_blocks"}}{Only legal for blocking-style designs that support
		#'       block-aware bootstrap resampling:
		#'       \code{DesignFixedBlocking}, \code{DesignFixedOptimalBlocks},
		#'       \code{DesignSeqOneByOneSPBR}, and \code{DesignFixedBlockedCluster}.
		#'       Resamples observational units within each observed block/stratum. For blocked
		#'       cluster designs this means resampling clusters within strata.}
		#'     \item{\code{"resample_blocks"}}{Only legal for the same blocking-style designs as
		#'       \code{"within_blocks"}. Resamples entire observed blocks/strata with replacement
		#'       rather than resampling units within each block.}
		#'   }
		#'   Any non-\code{NULL} value is rejected for designs outside that blocking family.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			private$active_resampling_operation = "non_param_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Bootstrap inference")
				private$assert_valid_bootstrap_type(bootstrap_type)
				assertCount(B, positive = TRUE); assertFlag(debug)
			}
			# Check cache (skipped in debug mode to always get fresh diagnostic results)
			cache_key = as.character(B)
			if (!isTRUE(debug) && !is.null(private$cached_values$boot_distr_cache[[cache_key]])) {
				return(private$cached_values$boot_distr_cache[[cache_key]])
			}
			if (private$verbose) cat("Computing bootstrap distribution...\n")
			# Duplicate objects for thread safety
			inf_template = self$duplicate()
			des_template = private$des_obj$duplicate()
			has_match_structure_local = private$has_match_structure
			if (!is.null(private$seed)) {
				had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				on.exit(
					if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv),
					add = TRUE
				)
				set.seed(private$seed)
			}
			boot_draws = replicate(
				as.integer(B),
				private$bootstrap_sample_indices(private$n, bootstrap_type),
				simplify = FALSE
			)
			run_one_boot_iter = function(worker_des, worker_inf) {
				boot_draw = boot_draws[[1L]]
				sub_inf = private$bootstrap_subset_inference(boot_draw, smooth = FALSE)
				if (is.null(sub_inf)) return(NA_real_)
				as.numeric(sub_inf$compute_estimate(estimate_only = TRUE))[1L]
			}
			if (isTRUE(debug)) {
				run_debug_boot_iter = function(boot_draw, worker_des = NULL, worker_inf = NULL, worker_state = NULL) {
					iter_warns = character(0)
					iter_errs = character(0)
					iter_val = withCallingHandlers(
						tryCatch({
							if (!is.null(worker_state)) {
								private$load_bootstrap_sample_into_worker(worker_state, boot_draw)
								private$compute_bootstrap_worker_estimate(worker_state)
							} else {
								sub_inf = private$bootstrap_subset_inference(boot_draw, smooth = FALSE)
								if (is.null(sub_inf)) NA_real_ else as.numeric(sub_inf$compute_estimate(estimate_only = TRUE))[1L]
							}
						}, error = function(e) { iter_errs <<- c(iter_errs, conditionMessage(e)); NA_real_ }),
						warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
					)
					list(
						val = as.numeric(iter_val)[1L],
						errors = iter_errs,
						warnings = iter_warns
					)
				}
					actual_debug_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
					chunk_n = max(1L, min(as.integer(actual_debug_cores), as.integer(B)))
					chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
					chunks = split(seq_len(B), chunk_id)
					run_debug_chunk = if (isTRUE(private$use_reusable_bootstrap_worker())) {
						function(idxs) {
							worker_state = private$create_bootstrap_worker_state()
							lapply(idxs, function(idx) run_debug_boot_iter(boot_draw = boot_draws[[idx]], worker_state = worker_state))
						}
					} else {
						function(idxs) {
							lapply(idxs, function(idx) {
								worker_des = des_template$duplicate()
								worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
								run_debug_boot_iter(boot_draw = boot_draws[[idx]], worker_des = worker_des, worker_inf = worker_inf)
							})
						}
					}
					debug_results = if (actual_debug_cores <= 1L) {
						run_debug_chunk(seq_len(B))
					} else {
						# par_lapply flattens its internal chunking, so if run_debug_chunk returns a list,
						# par_lapply returns a list of lists.
						res_raw = private$par_lapply(
							chunks,
							run_debug_chunk,
							n_cores = actual_debug_cores,
							budget = 1L,
							show_progress = show_progress
						)
						# Flatten chunks into a single list of results
						unlist(res_raw, recursive = FALSE, use.names = FALSE)
					}
					# Harden debug_results: remove any NULLs or non-list results that might 
					# have come from worker crashes in par_lapply.
					debug_results_surviving = Filter(function(x) is.list(x) && !is.null(x$val), debug_results)
					if (length(debug_results_surviving) < length(debug_results)) {
						warning("Some bootstrap iterations (", length(debug_results) - length(debug_results_surviving), ") were lost due to worker crashes or invalid results.")
					}
					debug_results = debug_results_surviving
				
					if (length(debug_results) == 0L) {
						if (should_run_asserts()) {
							stop("All bootstrap iterations failed or returned invalid results. Check for worker crashes or out-of-memory issues.")
						}
					}
					values = vapply(debug_results, function(x) {
						if (is.list(x) && !is.null(x[["val"]])) as.numeric(x[["val"]])[1L] else NA_real_
					}, numeric(1))
					errors_list = lapply(debug_results, function(x) {
						if (is.list(x) && !is.null(x[["errors"]])) x[["errors"]] else character(0)
					})
					warnings_list = lapply(debug_results, function(x) {
						if (is.list(x) && !is.null(x[["warnings"]])) x[["warnings"]] else character(0)
					})
					num_errors_vec = lengths(errors_list)
					num_warnings_vec = lengths(warnings_list)
					# Populate the normal cache so subsequent non-debug calls can reuse the values
					if (is.null(private$cached_values$boot_distr_cache)) private$cached_values$boot_distr_cache = list()
					private$cached_values$boot_distr_cache[[cache_key]] = values
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
			# Determine cores — warm-up guard: run one iteration serially to estimate per-iteration
			# cost, then only parallelize if computation outweighs overhead per worker.
			# For an existing persistent fork cluster the per-call overhead is ~10ms (socket
			# round-trip); if the cluster doesn't exist yet it will be created lazily (~300ms).
			# We run the warmup iteration TWICE and use the second timing. The first call often
			# pays cold-start penalties (C++ JIT, OS page-cache misses, R bytecode compilation)
			# that can inflate the estimate 5-15x vs steady-state cost, causing the guard to
			# wrongly choose parallel for small B values like r = 19.
			actual_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
			if (actual_cores > 1L) {
				do_warmup_iter = if (isTRUE(private$use_reusable_bootstrap_worker())) {
					function() {
						worker_state = private$create_bootstrap_worker_state()
						boot_draw = boot_draws[[1L]]
						private$load_bootstrap_sample_into_worker(worker_state, boot_draw)
						tryCatch(private$compute_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
					}
				} else {
					function() {
						boot_draw = boot_draws[[1L]]
						tryCatch({
							sub_inf = private$bootstrap_subset_inference(boot_draw, smooth = FALSE)
							if (is.null(sub_inf)) NA_real_ else as.numeric(sub_inf$compute_estimate(estimate_only = TRUE))[1L]
						}, error = function(e) NA_real_)
					}
				}
				system.time(do_warmup_iter())  # First call: discarded (cold-start overhead)
				t_boot_warmup = system.time(do_warmup_iter())[[3]]  # Second call: representative cost
				# Existing cluster: ~10ms round-trip overhead. No cluster yet: ~300ms lazy creation.
				fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.3
				if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores))
					actual_cores = 1L
			}
			boot_distr = if (isTRUE(private$use_reusable_bootstrap_worker())) {
				private$compute_bootstrap_distribution_with_reused_workers(
					boot_draws = boot_draws,
					actual_cores = actual_cores,
					show_progress = show_progress,
					bootstrap_type = bootstrap_type
				)
			} else {
				unlist(private$par_lapply(seq_along(boot_draws), function(idx) {
					boot_draw = boot_draws[[idx]]
					tryCatch({
						sub_inf = private$bootstrap_subset_inference(boot_draw, smooth = FALSE)
						if (is.null(sub_inf)) NA_real_ else as.numeric(sub_inf$compute_estimate(estimate_only = TRUE))[1L]
					}, error = function(e) NA_real_)
				}, n_cores = actual_cores, show_progress = show_progress,
				export_list = list(
					des_template = des_template,
					inf_template = inf_template,
					has_match_structure_local = has_match_structure_local,
					boot_draws = boot_draws
				)))
			}
			if (!is.numeric(boot_distr)) boot_distr = as.numeric(boot_distr)
			if (is.null(private$cached_values$boot_distr_cache)) private$cached_values$boot_distr_cache = list()
			private$cached_values$boot_distr_cache[[cache_key]] = boot_distr
			boot_distr
		},
		#' @description Computes a bootstrap-based two-sided p-value for the treatment effect.
		#'   Validity is asymptotic and design-dependent; for most covariate-adaptive designs the
		#'   p-value errs conservative. See the class-level section \emph{Design-specific validity
		#'   caveats for the nonparametric bootstrap}.
		#'
		#' @param delta  				Null hypothesis value. Default 0.
		#' @param B  					Number of bootstrap samples. Default 501.
		#' @param min_number_usable_samples Minimum number of finite bootstrap samples
		#'   required after filtering. Default 50. Must be smaller than \code{B}.
		#' @param type  				Bootstrap p-value type. Supported values are
		#'   \code{"percentile"} (default), \code{"symmetric"}, \code{"studentized"},
		#'   \code{"bootstrap-t"}, and \code{"bca"}.
		#'   \code{"percentile"}: shifts the bootstrap distribution to be centred at
		#'   \code{delta} and counts the two-tail proportion (Hall 1992).
		#'   \code{"symmetric"}: uses \eqn{|T^* - \bar{T}^*| \ge |t_{\rm obs} - \delta|}
		#'   for a symmetric one-sample test; recommended by Hall & Wilson (1991) when the
		#'   null distribution may be skewed. This pooled-tail test is offered only as a
		#'   p-value here, not as a confidence-interval \code{type} in
		#'   \code{compute_bootstrap_confidence_interval}: pooling both tails via
		#'   \eqn{|\cdot|} improves testing power (Hall & Wilson's original use case), but
		#'   inverting it unstudentized would add no value as an interval. The unstudentized
		#'   pivot is not asymptotically pivotal, so the resulting interval would have the
		#'   same first-order \eqn{O(n^\{-1/2\})} coverage error as \code{"percentile"}/
		#'   \code{"basic"}, while forcing symmetric bounds around a possibly skewed
		#'   bootstrap distribution --- strictly worse than \code{"percentile"}/\code{"basic"}
		#'   for shape-adaptivity, and strictly worse than \code{"symmetric-percentile-t"} for
		#'   accuracy, since studentizing (not the absolute-value pooling) is what buys the
		#'   \eqn{O(n^\{-1\})} improvement. The CI-worthy symmetric variant is therefore
		#'   \code{"symmetric-percentile-t"} (studentized pivot), not a plain \code{"symmetric"}
		#'   CI type.
		#'   \code{"studentized"} / \code{"bootstrap-t"}: pivots by the per-replicate
		#'   standard error, giving O(n^\{-1\}) error versus O(n^\{-1/2\}) for the percentile
		#'   method (Hall 1992; Davidson & MacKinnon 1999).
		#'   \code{"bca"}: bias-corrected and accelerated p-value via closed-form CI
		#'   inversion using the jackknife acceleration and bias-correction constants;
		#'   second-order accurate (Efron 1987; Efron & Tibshirani 1993).
		#' @param na.rm  				Remove non-finite bootstrap replicates. Default FALSE.
		#' @param show_progress  		A flag indicating whether a progress bar should be displayed.
		#'
		#' @return 	A bootstrap two-sided p-value.
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = FALSE, show_progress = TRUE, min_number_usable_samples = 5L){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
			}
			if (should_run_asserts()) {
				if (as.integer(B) <= as.integer(min_number_usable_samples)) {
					stop("B must be greater than min_number_usable_samples.", call. = FALSE)
				}
			}
			if (should_run_asserts()) {
				assertFlag(na.rm)
				assertFlag(show_progress)
			}
			type = tolower(private$get_bootstrap_type(type))
			if (should_run_asserts()) {
				assertChoice(type, c("percentile", "symmetric", "studentized", "bootstrap-t", "bca"))
			}
			est = as.numeric(self$compute_estimate())
			if (length(est) == 0L || !is.finite(est[1])) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate("bootstrap_original_estimate_unavailable")
				return(NA_real_)
			}
			est = est[1]
			private$clear_nonestimable_state()
			if (type %in% c("studentized", "bootstrap-t")) {
				boot_stats = private$approximate_bootstrap_statistics_beta_hat_T(
					B = B,
					na.rm = isTRUE(na.rm),
					require_se = TRUE,
					show_progress = show_progress
				)
				boot_distr = boot_stats$theta
			} else {
				boot_distr = self$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = show_progress)
				boot_stats = NULL
			}
			if (isTRUE(na.rm)) boot_distr = boot_distr[is.finite(boot_distr)]
			else if (any(!is.finite(boot_distr))) return(NA_real_)
			if (length(boot_distr) < as.integer(min_number_usable_samples)) {
				if (isTRUE(private$harden)) {
					if (!is.null(boot_stats)) {
						# studentized path: empty because SEs are missing, not because estimates are missing
						private$cache_nonestimable_se("bootstrap_too_few_finite_standard_errors")
					} else {
						private$cache_nonestimable_estimate("bootstrap_too_few_finite_estimates")
					}
				}
				return(NA_real_)
			}
			if (length(boot_distr) == 0L) return(NA_real_)
			if (private$bootstrap_estimates_extreme(boot_distr, est = est)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate("bootstrap_extreme_finite_estimates")
				return(NA_real_)
			}
			if (type == "percentile") {
				# Shift bootstrap distribution to be centred at delta (null hypothesis)
				boot_null = boot_distr - mean(boot_distr) + delta
				n_bs = length(boot_null)
				pval = min(1, max(2 / n_bs, 2 * min(
					sum(boot_null >= est) / n_bs,
					sum(boot_null <= est) / n_bs
				)))
				if (!is.finite(pval)) {
					if (isTRUE(private$harden)) private$cache_nonestimable_estimate("bootstrap_pvalue_unavailable")
					return(NA_real_)
				}
				pval
			} else if (type == "symmetric") {
				# Hall & Wilson (1991) symmetric test: pool both tails via absolute deviations.
				# p = P(|T* - mean(T*)| >= |t_obs - delta|)
				D_obs = abs(est - delta)
				D_boot = abs(boot_distr - mean(boot_distr))
				n_bs = length(D_boot)
				pval = min(1, max(2 / n_bs, mean(D_boot >= D_obs)))
				if (!is.finite(pval)) {
					if (isTRUE(private$harden)) private$cache_nonestimable_estimate("bootstrap_pvalue_unavailable")
					return(NA_real_)
				}
				pval
			} else if (type %in% c("studentized", "bootstrap-t")) {
				# Studentized (bootstrap-t) p-value: pivot by standard error.
				# t_obs = (est - delta) / se_hat
				# t*_b  = (T*_b  - est)  / se*_b    [centred at original estimate, not null]
				# p = P(|t*_b| >= |t_obs|)
				se_hat = tryCatch(private$infer_original_se(), error = function(e) NA_real_)
				if (!is.finite(se_hat) || se_hat <= 0) {
					if (isTRUE(private$harden)) private$cache_nonestimable_se("bootstrap_original_standard_error_unavailable")
					return(NA_real_)
				}
				t_obs = abs(est - delta) / se_hat
				se_boot = boot_stats$se
				t_boot = tryCatch(
					private$studentized_bootstrap_pivots(
						theta = boot_distr,
						se = se_boot,
						est = est,
						se_hat = se_hat,
						min_number_usable_samples = min_number_usable_samples,
						symmetric = TRUE
					),
					error = function(e) numeric(0)
				)
					if (length(t_boot) < as.integer(min_number_usable_samples)) {
						if (isTRUE(private$harden)) private$cache_nonestimable_se("bootstrap_unstable_studentized_standard_errors")
						return(NA_real_)
					}
					if (private$studentized_interval_scale_unstable(
						theta = boot_distr,
						se_hat = se_hat,
						pivots = t_boot,
						est = est,
						alpha = 0.05
					)) {
						if (isTRUE(private$harden)) private$cache_nonestimable_se("bootstrap_unstable_studentized_standard_errors")
						return(NA_real_)
					}
					min(1, max(2 / length(t_boot), mean(t_boot >= t_obs)))
				} else if (type == "bca") {
				# BCa p-value via closed-form CI inversion (Efron 1987; Efron & Tibshirani 1993).
				pval = tryCatch(
					private$pval_bca(boot_distr, est, delta),
					error = function(e) {
						if (isTRUE(private$harden)) {
							private$cache_nonestimable_estimate("bootstrap_pvalue_unavailable")
							return(NA_real_)
						}
						stop(e)
					}
				)
				if (!is.finite(pval) && isTRUE(private$harden)) {
					if (!isTRUE(self$is_nonestimable())) {
						private$cache_nonestimable_estimate("bootstrap_pvalue_unavailable")
					}
					return(NA_real_)
				}
				pval
			}
		},
		#' @description Computes a bootstrap-based confidence interval.
		#'   Coverage is asymptotic and design-dependent; for most covariate-adaptive designs the
		#'   interval errs conservative (over-coverage). See the class-level section
		#'   \emph{Design-specific validity caveats for the nonparametric bootstrap}.
		#'
		#' @param alpha  				The confidence level 1 - \code{alpha}. Default 0.05.
		#' @param B  					Number of bootstrap samples. Default 501.
		#' @param min_number_usable_samples Minimum number of finite bootstrap samples
		#'   required after filtering. Default 50. Must be smaller than \code{B}.
		#' @param type  				Bootstrap CI type. Supported values are
		#'   \code{"percentile"}, \code{"basic"}, \code{"studentized"},
		#'   \code{"bootstrap-t"}, \code{"symmetric-percentile-t"},
		#'   \code{"bca"}, \code{"prepivoted"}, \code{"double-bootstrap"},
		#'   \code{"calibrated"}, and \code{"smoothed"}.
		#'   There is no plain \code{"symmetric"} CI type (contrast with the \code{"symmetric"}
		#'   p-value type in \code{compute_bootstrap_two_sided_pval}): inverting the unstudentized
		#'   Hall & Wilson pooled-tail statistic would add no value as an interval, since it is not
		#'   asymptotically pivotal and so has the same first-order \eqn{O(n^\{-1/2\})} coverage
		#'   error as \code{"percentile"}/\code{"basic"}, while forcing symmetric bounds around a
		#'   possibly skewed bootstrap distribution --- strictly worse than \code{"percentile"}/
		#'   \code{"basic"} for shape-adaptivity, and strictly worse than
		#'   \code{"symmetric-percentile-t"} for accuracy, since studentizing (not the
		#'   absolute-value pooling) is what buys the \eqn{O(n^\{-1\})} improvement.
		#'   \code{"symmetric-percentile-t"} is the CI-worthy symmetric variant.
		#' @param na.rm                                   Remove non-finite bootstrap replicates.
		#'   Default TRUE. Non-finite replicates are always removed internally.
		#' @param show_progress  		Show progress bar.
		#'
		#' @return 	A bootstrap confidence interval.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Bootstrap inference")
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
			}
			if (should_run_asserts()) {
				if (as.integer(B) <= as.integer(min_number_usable_samples)) {
					stop("B must be greater than min_number_usable_samples.", call. = FALSE)
				}
			}
			type = tolower(private$get_bootstrap_type(type))
			if (should_run_asserts()) {
				assertChoice(type, c(
					"percentile", "basic", "studentized", "bootstrap-t",
					"symmetric-percentile-t", "bca", "prepivoted",
					"double-bootstrap", "calibrated", "smoothed"
				))
			}
			est = as.numeric(self$compute_estimate(estimate_only = FALSE))
			if (length(est) == 0L || !is.finite(est[1])) {
				if (isTRUE(private$harden)) {
					return(private$missing_bootstrap_ci(alpha, "bootstrap_original_estimate_unavailable", stage = "estimate"))
				}
				stop("Bootstrap confidence interval returned NA bounds")
			}
			est = est[1]
			# BCa only needs theta (not se), so route it through the parallel bootstrap path.
			# studentized/symmetric-percentile-t need se per replicate → serial statistics path.
			boot_stats = if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t", "prepivoted", "double-bootstrap", "calibrated", "smoothed")) {
				private$approximate_bootstrap_statistics_beta_hat_T(
					B = B,
					show_progress = show_progress,
					na.rm = na.rm,
					smooth = identical(type, "smoothed"),
					require_se = type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t")
				)
			} else {
				list(theta = self$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = show_progress), se = NULL)
			}
			boot_distr = boot_stats$theta
			boot_distr = boot_distr[is.finite(boot_distr)]
			if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t")) {
				se_boot = boot_stats$se
				ok = is.finite(boot_stats$theta) & is.finite(se_boot) & se_boot > 0
				if (sum(ok) < as.integer(min_number_usable_samples)) {
					if (isTRUE(private$harden)) {
						return(private$missing_bootstrap_ci(alpha, "bootstrap_too_few_finite_standard_errors", stage = "se"))
					}
					stop("Bootstrap confidence interval returned NA bounds")
				}
			} else if (length(boot_distr) < as.integer(min_number_usable_samples)) {
				if (isTRUE(private$harden)) {
					return(private$missing_bootstrap_ci(alpha, "bootstrap_too_few_finite_estimates", stage = "estimate"))
				}
				stop("Bootstrap confidence interval returned NA bounds")
			}
			if (private$bootstrap_estimates_extreme(boot_distr, est = est)) {
				if (isTRUE(private$harden)) {
					return(private$missing_bootstrap_ci(alpha, "bootstrap_extreme_finite_estimates", stage = "estimate"))
				}
				stop("Bootstrap estimates are numerically unstable.")
			}
			ci = tryCatch({
				if (type == "percentile") {
					private$ci_from_boot_distribution(boot_distr, alpha, "percentile")
				} else if (type == "basic") {
					private$ci_from_boot_distribution(boot_distr, alpha, "basic", est = est)
				} else if (type %in% c("studentized", "bootstrap-t")) {
					private$ci_studentized(boot_stats, alpha, est, min_number_usable_samples = min_number_usable_samples)
				} else if (type == "symmetric-percentile-t") {
					private$ci_symmetric_studentized(boot_stats, alpha, est, min_number_usable_samples = min_number_usable_samples)
				} else if (type == "bca") {
					private$ci_bca(boot_distr, alpha, est)
				} else if (type %in% c("prepivoted", "double-bootstrap", "calibrated")) {
					private$ci_calibrated_bootstrap(alpha, B, type, est, show_progress = show_progress, na.rm = na.rm)
				} else if (type == "smoothed") {
					private$ci_smoothed_bootstrap(alpha, B, est, show_progress = show_progress, na.rm = na.rm)
				} else {
					if (should_run_asserts()) {
						stop("Unsupported bootstrap CI type: ", type)
					}
				}
			}, error = function(e) {
				if (isTRUE(private$harden)) {
					stage = if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t")) "se" else "estimate"
					reason = if (identical(stage, "se")) "bootstrap_standard_error_ci_unavailable" else "bootstrap_ci_unavailable"
					return(private$missing_bootstrap_ci(alpha, reason, stage = stage))
				}
				stop(e)
			})
				if (length(ci) != 2L || !all(is.finite(ci[1:2]))) {
					if (isTRUE(private$harden)) {
						stage = if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t", "bca")) "se" else "estimate"
					reason = if (identical(stage, "se")) {
						if (type == "bca" && private$jackknife_block_size_gt_one_unsupported(unit = "auto")) {
							"jackknife_block_size_gt_one_not_supported"
						} else {
							"bootstrap_standard_error_ci_unavailable"
						}
					} else {
						"bootstrap_ci_unavailable"
					}
					return(private$missing_bootstrap_ci(alpha, reason, stage = stage))
					}
					stop("Bootstrap confidence interval returned NA bounds")
				}
				if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t") &&
				    private$studentized_interval_scale_unstable(theta = boot_stats$theta, ci = ci, est = est, alpha = alpha)) {
					if (isTRUE(private$harden)) {
						return(private$missing_bootstrap_ci(alpha, "bootstrap_unstable_studentized_interval", stage = "se"))
					}
					stop("Studentized bootstrap interval is numerically unstable.")
				}
				if (private$bootstrap_confidence_interval_extreme(ci, est = est)) {
					if (isTRUE(private$harden)) {
						stage = if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t")) "se" else "estimate"
						return(private$missing_bootstrap_ci(alpha, "bootstrap_extreme_confidence_interval", stage = stage))
					}
					stop("Bootstrap confidence interval is numerically unstable.")
				}
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				ci
			}
	),
	private = list(
		# Cache for bootstrap distributions
		boot_distr_cache = list(),
		jack_distr_cache = list(),
		bootstrap_extreme_estimate_threshold = 1e6,
		bootstrap_extreme_ci_width_threshold = 5,
		assert_valid_bootstrap_type = function(bootstrap_type){
			if (is.null(bootstrap_type)) return(invisible(NULL))
			if (should_run_asserts()) {
				assertChoice(bootstrap_type, c("within_blocks", "resample_blocks"))
			}
			valid_blocking_classes = c("DesignFixedBlocking", "DesignFixedOptimalBlocks", "DesignSeqOneByOneSPBR", "DesignFixedBlockedCluster")
			if (should_run_asserts()) {
				if (!any(vapply(valid_blocking_classes, function(cls) is(private$des_obj, cls), logical(1)))){
					stop("bootstrap_type can only be set for blocking designs: ", paste(valid_blocking_classes, collapse = ", "))
				}
			}
			invisible(NULL)
		},
		get_bootstrap_type = function(type) {
			if (!is.null(type)) return(type)
			edi_bootstrap_dispatch_policy(class(self), object = self)
		},
		resolve_jackknife_unit = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			if (unit != "auto") return(unit)
			design_obj = private$des_obj
			is_matching_design = is(design_obj, "DesignMatching") &&
				isTRUE(tryCatch(design_obj$is_matching_design(), error = function(e) FALSE))
			is_cluster_design = is(design_obj, "DesignFixedCluster") || is(design_obj, "DesignFixedBlockedCluster")
			is_blocking_design = is(design_obj, "DesignBlocking") &&
				isTRUE(tryCatch(design_obj$is_blocking_design(), error = function(e) FALSE))
			if (is_matching_design) {
				if (isTRUE(private$is_KK)) "matched_set" else "pair"
			} else if (is_cluster_design) {
				"cluster"
			} else if (is_blocking_design) {
				"block"
			} else {
				"observation"
			}
		},
		jackknife_block_size_gt_one_unsupported = function(unit = "auto"){
			resolved_unit = private$resolve_jackknife_unit(unit)
			design_obj = private$des_obj
			is_blocking_design = is(design_obj, "DesignBlocking") &&
				isTRUE(tryCatch(design_obj$is_blocking_design(), error = function(e) FALSE))
			if (resolved_unit == "block") return(FALSE)
			if (resolved_unit %in% c("pair", "matched_set")) return(FALSE)
			if (resolved_unit == "cluster" && !is_blocking_design) return(FALSE)
			if (!is_blocking_design) return(FALSE)
			block_ids = tryCatch(design_obj$get_block_ids(), error = function(e) NULL)
			if (is.null(block_ids)) return(FALSE)
			block_ids = as.integer(block_ids)
			block_ids = block_ids[is.finite(block_ids) & block_ids > 0L]
			if (!length(block_ids)) return(FALSE)
			any(as.integer(table(block_ids)) > 1L)
		},
		mark_jackknife_nonestimable_if_block_unsupported = function(unit = "auto"){
			if (!private$jackknife_block_size_gt_one_unsupported(unit = unit)) return(FALSE)
			private$cache_nonestimable_se("jackknife_block_size_gt_one_not_supported")
			TRUE
		},
		missing_bootstrap_ci = function(alpha, reason, stage = c("estimate", "se")){
			stage = match.arg(stage)
			if (identical(stage, "se")) private$cache_nonestimable_se(reason)
			else private$cache_nonestimable_estimate(reason)
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		bootstrap_estimates_extreme = function(theta, est = NA_real_, max_abs = private$bootstrap_extreme_estimate_threshold){
			theta = as.numeric(theta)
			theta = theta[is.finite(theta)]
			if (!length(theta)) return(FALSE)
			max_abs = as.numeric(max_abs)[1L]
			if (!is.finite(max_abs) || max_abs <= 0) max_abs = 1e6
			if (any(abs(theta) > max_abs)) return(TRUE)
			scale_ref = max(1, abs(as.numeric(est)[1L]), stats::median(abs(theta)), na.rm = TRUE)
			if (!is.finite(scale_ref) || scale_ref <= 0) scale_ref = 1
			theta_width = diff(stats::quantile(theta, probs = c(0.025, 0.975), names = FALSE, type = 8))
			is.finite(theta_width) && theta_width > max_abs * scale_ref
		},
		bootstrap_confidence_interval_extreme = function(ci, est = NA_real_, max_abs = private$bootstrap_extreme_estimate_threshold){
			ci = as.numeric(ci)
			if (length(ci) < 2L || !all(is.finite(ci[1:2]))) return(FALSE)
			max_abs = as.numeric(max_abs)[1L]
			if (!is.finite(max_abs) || max_abs <= 0) max_abs = 1e6
			if (any(abs(ci[1:2]) > max_abs)) return(TRUE)
			scale_ref = max(1, abs(as.numeric(est)[1L]), na.rm = TRUE)
			width = abs(diff(ci[1:2]))
			scaled_max_width = max_abs * scale_ref
			absolute_max_width = as.numeric(private$bootstrap_extreme_ci_width_threshold)[1L]
			if (!is.finite(absolute_max_width) || absolute_max_width <= 0) absolute_max_width = Inf
			max_width = min(scaled_max_width, absolute_max_width, na.rm = TRUE)
			if (!is.finite(max_width) || max_width <= 0) max_width = scaled_max_width
			is.finite(width) && width > max_width
		},
		supports_reusable_bootstrap_worker = function(){
			FALSE
		},
		create_bootstrap_worker_state = function(){
			NULL
		},
		create_design_backed_bootstrap_worker_state = function(){
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
			worker_priv$X = private$get_X()
			list(
				worker = worker,
				worker_priv = worker_priv,
				worker_des_priv = worker_des_priv,
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
				base_y_i_t_i = if (!is.null(source_des_priv$y_i_t_i)) source_des_priv$y_i_t_i else NULL,
				base_za_X_cov_all = if (!is.null(private$za_X_cov_all)) private$za_X_cov_all else NULL,
				base_za_Xzi_cov_all = if (!is.null(private$za_Xzi_cov_all)) private$za_Xzi_cov_all else NULL,
				n = private$n
			)
		},
		load_bootstrap_sample_into_worker = function(worker_state, indices){
			stop("Reusable bootstrap workers are not implemented for this class.")
		},
		load_bootstrap_sample_into_design_backed_worker = function(worker_state, indices){
			if (is.list(indices)) {
				m_vec_override = indices$m_vec_b
				indices = indices$i_b
			} else {
				m_vec_override = NULL
			}
			indices = as.integer(indices)
			w_priv = worker_state$worker_priv
			w_priv$X = if (!is.null(worker_state$base_X)) {
				worker_state$base_X[indices, , drop = FALSE]
			} else {
				NULL
			}
			w_priv$w = if (!is.null(worker_state$base_w)) as.numeric(worker_state$base_w[indices]) else NULL
			w_priv$y = if (!is.null(worker_state$base_y)) as.numeric(worker_state$base_y[indices]) else NULL
			w_priv$dead = if (!is.null(worker_state$base_dead)) as.numeric(worker_state$base_dead[indices]) else NULL
			w_priv$any_censoring = !is.null(w_priv$dead) && any(w_priv$dead == 0)
			w_priv$y_temp = w_priv$y
			if (!is.null(worker_state$base_m)) {
				w_priv$m = if (!is.null(m_vec_override)) as.integer(m_vec_override) else worker_state$base_m[indices]
			}
			w_priv$n = length(indices)
			w_priv$cached_values = list()
			w_priv$likelihood_null_warm_cache = list()
			w_priv$reduced_design_keep_cache = NULL
			w_priv$fixed_covariate_keep_cache = NULL
			w_priv$best_X_colnames = NULL
			w_priv$best_Xmm_colnames = NULL
			w_priv$fit_warm_start = worker_state$base_fit_warm_start
			w_priv$fit_warm_start_type = worker_state$base_fit_warm_start_type
			w_priv$fit_warm_start_fisher = worker_state$base_fit_warm_start_fisher
			w_priv$cached_mod = NULL
			
			w_priv$za_X_cov_all = if (!is.null(worker_state$base_za_X_cov_all)) {
				worker_state$base_za_X_cov_all[indices, , drop = FALSE]
			} else {
				NULL
			}
			w_priv$za_Xzi_cov_all = if (!is.null(worker_state$base_za_Xzi_cov_all)) {
				worker_state$base_za_Xzi_cov_all[indices, , drop = FALSE]
			} else {
				NULL
			}
			
			# Reset all private design matrix and covariate caches
			w_priv$cached_design_matrix = NULL
			w_priv$cached_w_for_design_matrix = NULL
			w_priv$cached_harden_for_design_matrix = NULL
			w_priv$cached_hardened_X_cov = NULL
			w_priv$cached_reduced_X = NULL
			w_priv$cached_X_full_for_reduced = NULL
			w_priv$cached_keep_for_reduced = NULL
			w_priv$cached_j_treat_for_reduced = NULL
			des_priv = worker_state$worker_des_priv
			if (!is.null(des_priv)) {
				des_priv$X = w_priv$X
				des_priv$t = length(indices)
				des_priv$n = length(indices)
				des_priv$w = w_priv$w
				des_priv$y = w_priv$y
				des_priv$dead = w_priv$dead
				if (!is.null(worker_state$base_m)) des_priv$m = w_priv$m
				
				# Subset Xraw, Ximp, and y_i_t_i in the worker's design
				subset_field = function(x) {
					if (is.null(x)) return(NULL)
					if (is.data.frame(x) || is.matrix(x)) return(x[indices, , drop = FALSE])
					if (is.list(x) && !is.data.frame(x)) return(x[indices])
					if (is.atomic(x) && length(x) >= max(indices)) return(x[indices])
					x
				}
				des_priv$Xraw = subset_field(worker_state$base_Xraw)
				des_priv$Ximp = subset_field(worker_state$base_Ximp)
				if (!is.null(worker_state$base_y_i_t_i)) des_priv$y_i_t_i = worker_state$base_y_i_t_i[indices]

				des_priv$all_subject_data_cache = list()
				des_priv$lin_centered_covariates = NULL
				if (is.function(des_priv$reset_matching_caches)) des_priv$reset_matching_caches()
			}
		},
		compute_bootstrap_worker_estimate = function(worker_state){
			stop("Reusable bootstrap workers are not implemented for this class.")
		},
		compute_bootstrap_worker_estimate_via_compute_treatment_estimate = function(worker_state){
			theta = as.numeric(worker_state$worker$compute_estimate(estimate_only = TRUE))[1L]
			if (is.function(worker_state$worker$is_nonestimable) &&
			    isTRUE(worker_state$worker$is_nonestimable("estimate"))){
				return(NA_real_)
			}
			theta
		},
		compute_bootstrap_distribution_with_reused_workers = function(boot_draws, actual_cores, show_progress = FALSE, bootstrap_type = NULL){
			B = length(boot_draws)
			chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
			chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
			chunks = split(seq_len(B), chunk_id)
			run_chunk = function(idxs) {
				worker_state = private$create_bootstrap_worker_state()
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					boot_draw = boot_draws[[idxs[[k]]]]
					out[k] = tryCatch({
						private$load_bootstrap_sample_into_worker(worker_state, boot_draw)
						private$compute_bootstrap_worker_estimate(worker_state)
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
		},
		compute_jackknife_distribution_with_reused_workers = function(deletion_draws, actual_cores, show_progress = FALSE){
			n_draws = length(deletion_draws)
			chunk_n = max(1L, min(as.integer(actual_cores), as.integer(n_draws)))
			chunk_id = ceiling(seq_len(n_draws) / ceiling(n_draws / chunk_n))
			chunks = split(seq_len(n_draws), chunk_id)
			run_chunk = function(idxs) {
				worker_state = private$create_bootstrap_worker_state()
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					out[k] = tryCatch({
						private$load_bootstrap_sample_into_worker(worker_state, deletion_draws[[idxs[k]]])
						private$compute_bootstrap_worker_estimate(worker_state)
					}, error = function(e) NA_real_)
				}
				out
			}
			if (actual_cores <= 1L) {
				return(as.numeric(run_chunk(seq_len(n_draws))))
			}
			as.numeric(unlist(private$par_lapply(
				chunks,
				run_chunk,
				n_cores = actual_cores,
				budget = 1L,
				show_progress = show_progress
			), use.names = FALSE))
		},
		bootstrap_sample_indices = function(n, bootstrap_type = NULL){
			if (!is.null(private$des_obj)){
				return(private$des_obj_priv_int$draw_bootstrap_indices(bootstrap_type))
			}
			list(i_b = sample.int(n, n, replace = TRUE), m_vec_b = NULL)
		},
		renumber_match_ids = function(m_vec){
			if (is.null(m_vec)) return(NULL)
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			pos = sort(unique(m_vec[m_vec > 0L]))
			if (!length(pos)) return(m_vec)
			map = seq_along(pos)
			names(map) = as.character(pos)
			out = integer(length(m_vec))
			is_match = m_vec > 0L
			out[is_match] = unname(map[as.character(m_vec[is_match])])
			out
		},
		get_cluster_jackknife_ids = function(design_obj){
			des_priv = private$des_obj_priv_int
			if (is(design_obj, "DesignFixedBlockedCluster") || is(design_obj, "DesignFixedCluster")) {
				cluster_col = des_priv$cluster_col
				Xraw = des_priv$Xraw
				n = private$des_obj$get_n()
				if (!is.null(cluster_col) && !is.null(Xraw) && cluster_col %in% names(Xraw)) {
					return(as.character(Xraw[seq_len(n), ][[cluster_col]]))
				}
			}
			NULL
		},
		build_jackknife_deletion_draws = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			n = private$des_obj$get_n()
			all_idx = seq_len(n)
			design_obj = private$des_obj
			is_blocking_design = is(design_obj, "DesignBlocking") &&
				isTRUE(tryCatch(design_obj$is_blocking_design(), error = function(e) FALSE))
			is_matching_design = is(design_obj, "DesignMatching") &&
				isTRUE(tryCatch(design_obj$is_matching_design(), error = function(e) FALSE))
			is_cluster_design = is(design_obj, "DesignFixedCluster") || is(design_obj, "DesignFixedBlockedCluster")
			if (unit == "auto") {
				unit = if (is_matching_design) {
					if (isTRUE(private$is_KK)) "matched_set" else "pair"
				} else if (is_cluster_design) {
					"cluster"
				} else if (is_blocking_design) {
					"block"
				} else {
					"observation"
				}
			}
			private$assert_jackknife_supported(unit = unit)
			if (unit == "cluster") {
				cluster_ids = private$get_cluster_jackknife_ids(design_obj)
				if (!is.null(cluster_ids)) {
					unique_clusters = unique(cluster_ids[!is.na(cluster_ids)])
					if (length(unique_clusters) > 0L) {
						return(lapply(unique_clusters, function(cluster_id) {
							keep_idx = all_idx[cluster_ids != cluster_id]
							list(i_b = keep_idx, m_vec_b = NULL)
						}))
					}
				}
			}
			if (unit == "block") {
				block_ids = tryCatch(design_obj$get_block_ids(), error = function(e) NULL)
				if (!is.null(block_ids)) {
					block_ids = as.integer(block_ids)
					unique_blocks = sort(unique(block_ids[is.finite(block_ids) & block_ids > 0L]))
					if (length(unique_blocks) >= 2L) {
						return(lapply(unique_blocks, function(block_id) {
							keep_idx = all_idx[block_ids != block_id]
							list(i_b = keep_idx, m_vec_b = NULL)
						}))
					}
				}
				unit = "observation"
			}
			if (unit == "observation") {
				return(lapply(all_idx, function(i) all_idx[all_idx != i]))
			}
			if (unit %in% c("pair", "matched_set") && !is_matching_design) {
				return(lapply(all_idx, function(i) all_idx[all_idx != i]))
			}
			des_priv = private$des_obj_priv_int
			des_priv$init_matching_bootstrap_structure()
			pair_rows = des_priv$boot_pair_rows
			i_reservoir = des_priv$boot_i_reservoir
			m_vec_full = if (!is.null(private$m)) private$m else des_priv$m
			if (is.null(m_vec_full)) {
				return(lapply(all_idx, function(i) all_idx[all_idx != i]))
			}
			deletion_draws = list()
			if (isTRUE(private$is_KK) && length(i_reservoir) > 0L) {
				for (i in seq_along(i_reservoir)) {
					keep_idx = all_idx[all_idx != i_reservoir[i]]
					deletion_draws[[length(deletion_draws) + 1L]] = list(
						i_b = keep_idx,
						m_vec_b = private$renumber_match_ids(m_vec_full[keep_idx])
					)
				}
			}
			if (!is.null(pair_rows) && nrow(pair_rows) > 0L) {
				for (pair_idx in seq_len(nrow(pair_rows))) {
					drop_idx = as.integer(pair_rows[pair_idx, ])
					keep_idx = all_idx[!(all_idx %in% drop_idx)]
					deletion_draws[[length(deletion_draws) + 1L]] = list(
						i_b = keep_idx,
						m_vec_b = private$renumber_match_ids(m_vec_full[keep_idx])
					)
				}
			}
			if (!length(deletion_draws)) {
				return(lapply(all_idx, function(i) all_idx[all_idx != i]))
			}
			deletion_draws
		},
		bootstrap_subset_inference = function(boot_draw, smooth = FALSE){
			# boot_draw is list(i_b, m_vec_b) from des_obj$draw_bootstrap_indices()
			# For backward compatibility also accept a bare integer vector.
			if (is.list(boot_draw)){
				indices = as.integer(boot_draw$i_b)
				m_vec_b = boot_draw$m_vec_b
			} else {
				indices = as.integer(boot_draw)
				m_vec_b = NULL
			}
			if (length(indices) == 0L) return(NULL)
			orig_des = private$des_obj
			orig_des_priv = private$des_obj_priv_int
			sub_des = orig_des$duplicate(verbose = FALSE)
			sub_des_priv = sub_des$.__enclos_env__$private
			subset_field = function(x){
				if (is.null(x)) return(NULL)
				if (is.data.frame(x) || is.matrix(x)) return(x[indices, , drop = FALSE])
				if (is.list(x) && !is.data.frame(x)) return(x[indices])
				if (is.atomic(x) && length(x) >= max(indices)) return(x[indices])
				x
			}
			if (!is.null(orig_des_priv$Xraw)) sub_des_priv$Xraw = subset_field(orig_des_priv$Xraw)
			if (!is.null(orig_des_priv$Ximp)) sub_des_priv$Ximp = subset_field(orig_des_priv$Ximp)
			if (!is.null(orig_des_priv$X)) sub_des_priv$X = subset_field(orig_des_priv$X)
			if (!is.null(orig_des_priv$w)) sub_des_priv$w = as.numeric(orig_des_priv$w[indices])
			if (!is.null(orig_des_priv$y)) sub_des_priv$y = as.numeric(orig_des_priv$y[indices])
			if (!is.null(orig_des_priv$dead)) sub_des_priv$dead = as.numeric(orig_des_priv$dead[indices])
			# Use Design-provided m_vec_b (pair-aware) if available; otherwise subset original m
			if (!is.null(orig_des_priv$m)) {
				sub_des_priv$m = if (!is.null(m_vec_b)) as.integer(m_vec_b) else as.integer(orig_des_priv$m[indices])
			}
			if (is.function(sub_des_priv$reset_matching_caches)) sub_des_priv$reset_matching_caches()
			if (!is.null(orig_des_priv$y_i_t_i)) sub_des_priv$y_i_t_i = orig_des_priv$y_i_t_i[indices]
			sub_des_priv$all_subject_data_cache = list()
			sub_des_priv$p_raw_t = if (!is.null(sub_des_priv$Xraw)) ncol(sub_des_priv$Xraw) else NULL
			sub_des_priv$t = length(indices)
			sub_des_priv$n = length(indices)
			sub_des_priv$fixed_sample = TRUE
			# Reset bootstrap structure cache — indices changed, pair structure no longer valid
			sub_des_priv$boot_pair_rows   = NULL
			sub_des_priv$boot_i_reservoir = NULL
			sub_des_priv$boot_n_reservoir = NULL
			# Reset model-specific design caches that hold n-row matrices — must be recomputed
			# on the subset (e.g. Lin centered covariates, which are cached as an n-row Xc).
			sub_des_priv$lin_centered_covariates = NULL
			if (smooth && !is.null(sub_des_priv$y) && private$des_obj_priv_int$response_type == "continuous") {
				sd_y = stats::sd(as.numeric(sub_des_priv$y), na.rm = TRUE)
				if (is.finite(sd_y) && sd_y > 0) {
					sub_des_priv$y = as.numeric(sub_des_priv$y) + stats::rnorm(length(sub_des_priv$y), 0, sd_y / sqrt(max(1, length(indices))))
				}
			}
			sub_inf = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			sub_inf_priv = sub_inf$.__enclos_env__$private
			sub_inf_priv$des_obj = sub_des
			sub_inf_priv$des_obj_priv_int = sub_des_priv
			sub_inf_priv$X = if (!is.null(private$X)) subset_field(private$X) else NULL
			sub_inf_priv$y = sub_des_priv$y
			sub_inf_priv$y_temp = sub_des_priv$y
			sub_inf_priv$w = sub_des_priv$w
			sub_inf_priv$dead = sub_des_priv$dead
			sub_inf_priv$n = length(indices)
			sub_inf_priv$cached_values = list()
			sub_inf_priv$cached_values$rand_distr_cache = list()
			sub_inf_priv$cached_values$m_cache = list()
			sub_inf_priv$reduced_design_keep_cache = NULL
			sub_inf_priv$fixed_covariate_keep_cache = NULL
			sub_inf_priv$cached_mod = NULL
			
			# Reset all private design matrix and covariate caches
			sub_inf_priv$cached_design_matrix = NULL
			sub_inf_priv$cached_w_for_design_matrix = NULL
			sub_inf_priv$cached_harden_for_design_matrix = NULL
			sub_inf_priv$cached_hardened_X_cov = NULL
			sub_inf_priv$cached_reduced_X = NULL
			sub_inf_priv$cached_X_full_for_reduced = NULL
			sub_inf_priv$cached_keep_for_reduced = NULL
			sub_inf_priv$cached_j_treat_for_reduced = NULL

			if (!is.null(sub_des_priv$m)) sub_inf_priv$m = sub_des_priv$m
			sub_inf
		},
		bootstrap_replication_stats = function(boot_draw, smooth = FALSE, require_se = FALSE){
			sub_inf = private$bootstrap_subset_inference(boot_draw, smooth = smooth)
			if (is.null(sub_inf)) return(c(theta = NA_real_, se = NA_real_))
			tryCatch({
				# Use the return value of compute_estimate() for theta.
				# Reading cached_values$beta_hat_T is unreliable: some classes (e.g.
				# GComp, KMDiff) return estimates directly without storing beta_hat_T.
				theta = as.numeric(sub_inf$compute_estimate(estimate_only = !isTRUE(require_se)))[1L]
				se = if (isTRUE(require_se)) as.numeric(sub_inf$.__enclos_env__$private$cached_values$s_beta_hat_T)[1L] else NA_real_
				if (is.function(sub_inf$is_nonestimable) && isTRUE(sub_inf$is_nonestimable("estimate"))) {
					return(c(theta = NA_real_, se = NA_real_))
				}
				if (!is.finite(theta)) theta = NA_real_
				if (isTRUE(require_se) && is.function(sub_inf$is_nonestimable) && isTRUE(sub_inf$is_nonestimable("se"))) {
					se = NA_real_
				}
				if (isTRUE(require_se) && !is.finite(se)) se = NA_real_
				c(theta = theta, se = se)
			}, error = function(e) c(theta = NA_real_, se = NA_real_))
		},
		approximate_bootstrap_statistics_beta_hat_T = function(B = 501, show_progress = TRUE, na.rm = TRUE, smooth = FALSE, require_se = FALSE){
			private$active_resampling_operation = "non_param_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (should_run_asserts()) {
				assertCount(B, positive = TRUE)
				assertFlag(require_se)
			}
			if (!isTRUE(require_se) && !isTRUE(smooth)) {
				theta = self$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = show_progress)
				if (isTRUE(na.rm)) theta = theta[is.finite(theta)]
				if (length(theta) == 0L) {
					return(list(theta = numeric(0), se = numeric(0)))
				}
				return(list(theta = theta, se = rep(NA_real_, length(theta))))
			}
			n = private$des_obj$get_n()
			stats_mat = matrix(NA_real_, nrow = B, ncol = 2L)
			pb = NULL
			if (isTRUE(show_progress) && B > 1L) {
				pb = utils::txtProgressBar(min = 0, max = B, style = 3)
				on.exit(try(close(pb), silent = TRUE), add = TRUE)
			}
			for (b in seq_len(B)) {
				idx = private$bootstrap_sample_indices(n)
				stats_mat[b, ] = private$bootstrap_replication_stats(idx, smooth = smooth, require_se = require_se)  # idx is a boot_draw list
				if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
			}
			if (isTRUE(na.rm)) {
				ok = is.finite(stats_mat[, 1L])
				if (isTRUE(require_se)) ok = ok & is.finite(stats_mat[, 2L]) & stats_mat[, 2L] > 0
				stats_mat = stats_mat[ok, , drop = FALSE]
			}
			if (nrow(stats_mat) == 0L) {
				return(list(theta = numeric(0), se = numeric(0)))
			}
			list(theta = stats_mat[, 1L], se = stats_mat[, 2L])
		},
		approximate_jackknife_distribution_beta_hat_T_private = function(unit = "auto"){
			private$active_resampling_operation = "jackknife"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			unit = private$normalize_jackknife_unit(unit)
			if (private$mark_jackknife_nonestimable_if_block_unsupported(unit = unit)) {
				return(numeric(0))
			}
			deletion_draws = private$build_jackknife_deletion_draws(unit = unit)
			n_draws = length(deletion_draws)
			if (n_draws <= 1L) return(numeric(0))
			actual_cores = private$effective_parallel_cores("jackknife", self$num_cores)
			cache_key = paste0(unit, "::", as.character(n_draws))
			if (!is.null(private$cached_values$jack_distr_cache[[cache_key]])) {
				return(private$cached_values$jack_distr_cache[[cache_key]])
			}
			jack = if (isTRUE(private$use_reusable_bootstrap_worker())) {
				private$compute_jackknife_distribution_with_reused_workers(
					deletion_draws = deletion_draws,
					actual_cores = actual_cores,
					show_progress = FALSE
				)
			} else {
				unlist(private$par_lapply(seq_len(n_draws), function(i) {
					sub_inf = private$bootstrap_subset_inference(deletion_draws[[i]], smooth = FALSE)
					if (is.null(sub_inf)) return(NA_real_)
					tryCatch({
						theta = as.numeric(sub_inf$compute_estimate(estimate_only = TRUE))[1L]
						if (is.finite(theta)) theta else NA_real_
					}, error = function(e) {
						message("Jackknife error at iteration ", i, ": ", e$message)
						NA_real_
					})
				}, n_cores = actual_cores, show_progress = FALSE), use.names = FALSE)
			}
			jack = as.numeric(jack)
			private$cached_values$jack_distr_cache[[cache_key]] = jack
			as.numeric(jack)
		},
			ci_from_boot_distribution = function(boot_distr, alpha, type, est = NULL){
			type = tolower(type)
			if (should_run_asserts()) {
				if (length(boot_distr) == 0L) stop("Bootstrap confidence interval returned NA bounds")
			}
			if (is.null(est)) est = as.numeric(self$compute_estimate())[1]
			if (type == "percentile") {
				stats::quantile(boot_distr, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8)
			} else {
				2 * est - stats::quantile(boot_distr, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
				}
			},
			studentized_interval_scale_unstable = function(theta, ci = NULL, se_hat = NULL, pivots = NULL, est = 0, alpha = 0.05, max_width_ratio = 5){
				theta = as.numeric(theta)
				theta = theta[is.finite(theta)]
				if (length(theta) < 5L) return(FALSE)
				theta_width = diff(stats::quantile(theta, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8))
				scale_ref = max(as.numeric(theta_width), 1e-8, 1e-3 * max(1, abs(as.numeric(est)[1L])), na.rm = TRUE)
				if (!is.finite(scale_ref) || scale_ref <= 0) return(FALSE)
				if (!is.null(ci)) {
					ci = as.numeric(ci)
					if (length(ci) < 2L || !all(is.finite(ci[1:2]))) return(FALSE)
					width = abs(diff(ci[1:2]))
				} else {
					if (is.null(pivots) || !is.finite(se_hat) || se_hat <= 0) return(FALSE)
					pivots = as.numeric(pivots)
					pivots = pivots[is.finite(pivots)]
					if (length(pivots) < 5L) return(FALSE)
					width = 2 * stats::quantile(abs(pivots), probs = 1 - alpha / 2, names = FALSE, type = 8) * as.numeric(se_hat)[1L]
				}
				is.finite(width) && width > max_width_ratio * scale_ref
			},
			studentized_bootstrap_pivots = function(theta, se, est, se_hat, min_number_usable_samples = 10L, symmetric = FALSE){
			se = as.numeric(se)
			theta = as.numeric(theta)
			min_number_usable_samples = as.integer(min_number_usable_samples)
			se_pos = se[is.finite(se) & se > 0]
			if (!length(se_pos) || !is.finite(se_hat) || se_hat <= 0) {
				stop("Studentized bootstrap requires finite positive standard errors.")
			}
			se_ref = stats::median(se_pos)
			se_floor = max(.Machine$double.eps, 1e-6 * as.numeric(se_hat), 1e-6 * as.numeric(se_ref))
			ok = is.finite(theta) & is.finite(se) & se > se_floor
			pivots = (theta[ok] - est) / se[ok]
			pivots = pivots[is.finite(pivots)]
			if (isTRUE(symmetric)) pivots = abs(pivots)
			if (length(pivots) < min_number_usable_samples) {
				stop("Studentized bootstrap returned too few stable standard errors.")
			}
			if (stats::quantile(abs(pivots), probs = 0.975, names = FALSE, type = 8) > 50) {
				stop("Studentized bootstrap pivots are numerically unstable.")
			}
			pivots
		},
		ci_studentized = function(boot_stats, alpha, est, min_number_usable_samples = 10L){
			se_hat = private$infer_original_se()
			if (should_run_asserts()) {
				if (!is.finite(se_hat) || se_hat <= 0) stop("Studentized bootstrap CI requires a finite standard error.")
			}
			t_vals = private$studentized_bootstrap_pivots(
				theta = boot_stats$theta,
				se = boot_stats$se,
				est = est,
				se_hat = se_hat,
				min_number_usable_samples = min_number_usable_samples
			)
			q = stats::quantile(t_vals, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
			c(est - q[1L] * se_hat, est - q[2L] * se_hat)
		},
		ci_symmetric_studentized = function(boot_stats, alpha, est, min_number_usable_samples = 10L){
			se_hat = private$infer_original_se()
			if (should_run_asserts()) {
				if (!is.finite(se_hat) || se_hat <= 0) stop("Symmetric percentile-t bootstrap CI requires a finite standard error.")
			}
			t_vals = private$studentized_bootstrap_pivots(
				theta = boot_stats$theta,
				se = boot_stats$se,
				est = est,
				se_hat = se_hat,
				min_number_usable_samples = min_number_usable_samples,
				symmetric = TRUE
			)
			q = stats::quantile(t_vals, probs = 1 - alpha, names = FALSE, type = 8)
			c(est - q[1L] * se_hat, est + q[1L] * se_hat)
		},
		ci_bca = function(boot_distr, alpha, est){
			jack = private$approximate_jackknife_distribution_beta_hat_T_private(unit = "auto")
			jack = jack[is.finite(jack)]
			if (length(jack) < 2L) {
				reason = if (private$jackknife_block_size_gt_one_unsupported(unit = "auto")) {
					"jackknife_block_size_gt_one_not_supported"
				} else {
					"bootstrap_bca_jackknife_unavailable"
				}
				private$cache_nonestimable_se(reason)
				return(c(NA_real_, NA_real_))
			}
			if (should_run_asserts()) {
				if (length(jack) < 2L) stop("BCa interval requires jackknife estimates.")
			}
			p_less = mean(boot_distr < est)
			p_less = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_less))
			z0 = stats::qnorm(p_less)
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			if (!is.finite(z0) || !is.finite(a) || abs(z0) > 2.5 || abs(a) > 1) {
				private$cache_nonestimable_se("bootstrap_bca_unstable_bias_or_acceleration")
				return(c(NA_real_, NA_real_))
			}
			alpha_vec = c(alpha / 2, 1 - alpha / 2)
			z_alpha = stats::qnorm(alpha_vec)
			denom = 1 - a * (z0 + z_alpha)
			if (any(!is.finite(denom)) || any(abs(denom) < sqrt(.Machine$double.eps))) {
				private$cache_nonestimable_se("bootstrap_bca_unstable_bias_or_acceleration")
				return(c(NA_real_, NA_real_))
			}
			adj = stats::pnorm(z0 + (z0 + z_alpha) / denom)
			prob_eps = 1 / (length(boot_distr) + 1)
			adj = pmin(1 - prob_eps, pmax(prob_eps, adj))
			adj = sort(adj)
			if (any(adj <= 2 * prob_eps) || any(adj >= 1 - 2 * prob_eps)) {
				private$cache_nonestimable_se("bootstrap_bca_adjustment_on_boundary")
				return(c(NA_real_, NA_real_))
			}
			if (diff(adj) < prob_eps) {
				return(stats::quantile(boot_distr, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8))
			}
			stats::quantile(boot_distr, probs = adj, names = FALSE, type = 8)
		},
		ci_calibrated_bootstrap = function(alpha, B, type, est, show_progress = TRUE, na.rm = TRUE){
			n_outer = max(25L, min(as.integer(B), 101L))
			n_inner = max(25L, min(as.integer(B), 51L))
			alpha_grid = unique(pmin(0.49, pmax(.Machine$double.eps, c(alpha / 4, alpha / 2, alpha, min(0.25, alpha * 1.5), min(0.4, alpha * 2)))))
			coverage = rep(NA_real_, length(alpha_grid))
			for (j in seq_along(alpha_grid)) {
				a = alpha_grid[j]
				covered = logical(n_outer)
				for (b in seq_len(n_outer)) {
					idx = private$bootstrap_sample_indices(private$des_obj$get_n())  # returns boot_draw list
					outer_inf = private$bootstrap_subset_inference(idx, smooth = identical(type, "smoothed"))
					if (is.null(outer_inf)) next
					inner_boot = outer_inf$approximate_bootstrap_distribution_beta_hat_T(B = n_inner, show_progress = FALSE, debug = FALSE)
					inner_boot = inner_boot[is.finite(inner_boot)]
					if (length(inner_boot) < 10L) next
					ci_inner = private$ci_from_boot_distribution(inner_boot, a, "percentile")
					covered[b] = is.finite(ci_inner[1]) && is.finite(ci_inner[2]) && est >= ci_inner[1] && est <= ci_inner[2]
				}
				coverage[j] = mean(covered, na.rm = TRUE)
			}
			target_coverage = 1 - alpha
			idx_best = which.min(abs(coverage - target_coverage))
			alpha_adj = alpha_grid[idx_best]
			outer_boot = private$approximate_bootstrap_statistics_beta_hat_T(
				B = B,
				show_progress = show_progress,
				na.rm = na.rm,
				smooth = identical(type, "smoothed")
			)$theta
			outer_boot = outer_boot[is.finite(outer_boot)]
			if (should_run_asserts()) {
				if (length(outer_boot) < 10L) stop("Calibrated bootstrap CI returned too few finite bootstrap draws.")
			}
			ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "percentile")
			if (identical(type, "double-bootstrap")) {
				ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "basic", est = est)
			}
			if (identical(type, "prepivoted")) {
				ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "percentile")
			}
			if (identical(type, "calibrated")) {
				ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "percentile")
			}
			ci
		},
		ci_smoothed_bootstrap = function(alpha, B, est, show_progress = TRUE, na.rm = TRUE){
			boot_stats = private$approximate_bootstrap_statistics_beta_hat_T(B = B, show_progress = show_progress, na.rm = na.rm, smooth = TRUE)
			boot_distr = boot_stats$theta[is.finite(boot_stats$theta)]
			if (should_run_asserts()) {
				if (length(boot_distr) < 10L) stop("Smoothed bootstrap CI returned too few finite bootstrap draws.")
			}
			private$ci_from_boot_distribution(boot_distr, alpha, "percentile", est = est)
		},
			infer_original_se = function(){
				se_current = as.numeric(private$cached_values$s_beta_hat_T %||% NA_real_)[1L]
				if (is.finite(se_current) && se_current > 0) {
					return(se_current)
				}
				fresh = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
				fresh$.__enclos_env__$private$cached_values = list()
				tryCatch({
				fresh$compute_estimate(estimate_only = FALSE)
				as.numeric(fresh$.__enclos_env__$private$cached_values$s_beta_hat_T)[1]
			}, error = function(e) NA_real_)
		},
		# BCa p-value via closed-form CI inversion (Efron 1987; Efron & Tibshirani 1993).
		# Derives the bias-correction (z0) and acceleration (a) constants from the bootstrap
		# distribution and jackknife estimates, then computes the BCa-adjusted CDF value at
		# delta_0.  The two-sided p-value is  2 * min(Phi(adj_z), 1 - Phi(adj_z))  where
		#   z_delta = Phi^{-1}(F_boot(delta_0)),   s = z_delta - z0,
		#   adj_z   = s / (1 + a*s) - z0.
		pval_bca = function(boot_distr, est, delta){
			jack = private$approximate_jackknife_distribution_beta_hat_T_private(unit = "auto")
			jack = jack[is.finite(jack)]
			if (length(jack) < 2L) {
				reason = if (private$jackknife_block_size_gt_one_unsupported(unit = "auto")) {
					"jackknife_block_size_gt_one_not_supported"
				} else {
					"bootstrap_bca_jackknife_unavailable"
				}
				private$cache_nonestimable_se(reason)
				return(NA_real_)
			}
			if (should_run_asserts()) {
				if (length(jack) < 2L) stop("BCa p-value requires jackknife estimates.")
			}
			p_less = mean(boot_distr < est)
			p_less = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_less))
			z0 = stats::qnorm(p_less)
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			if (!is.finite(z0) || !is.finite(a) || abs(z0) > 2.5 || abs(a) > 1) {
				private$cache_nonestimable_se("bootstrap_bca_unstable_bias_or_acceleration")
				return(NA_real_)
			}
			z_alpha = stats::qnorm(c(0.025, 0.975))
			denom_ci = 1 - a * (z0 + z_alpha)
			if (any(!is.finite(denom_ci)) || any(abs(denom_ci) < sqrt(.Machine$double.eps))) {
				private$cache_nonestimable_se("bootstrap_bca_unstable_bias_or_acceleration")
				return(NA_real_)
			}
			prob_eps = 1 / (length(boot_distr) + 1)
			adj_ci = sort(stats::pnorm(z0 + (z0 + z_alpha) / denom_ci))
			if (any(!is.finite(adj_ci)) || any(adj_ci <= 2 * prob_eps) || any(adj_ci >= 1 - 2 * prob_eps)) {
				private$cache_nonestimable_se("bootstrap_bca_adjustment_on_boundary")
				return(NA_real_)
			}
			p_delta = mean(boot_distr < delta)
			p_delta = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_delta))
			z_delta = stats::qnorm(p_delta)
			s = z_delta - z0
			denom = 1 + a * s
			if (!is.finite(denom) || abs(denom) < sqrt(.Machine$double.eps)) {
				private$cache_nonestimable_se("bootstrap_bca_unstable_bias_or_acceleration")
				return(NA_real_)
			}
			adj_z = s / denom - z0
			if (!is.finite(adj_z) || abs(adj_z) > 8) {
				private$cache_nonestimable_se("bootstrap_bca_adjustment_on_boundary")
				return(NA_real_)
			}
			p_raw = min(1, 2 * min(stats::pnorm(adj_z), 1 - stats::pnorm(adj_z)))
			min(1, max(2 / length(boot_distr), p_raw))
		}
	)
)
