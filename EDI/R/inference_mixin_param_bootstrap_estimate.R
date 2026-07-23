#' Mixin for Parametric-Bootstrap Point-Estimate Bias Correction
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' providing the private implementation of the parametric-bootstrap
#' bias-corrected point estimate, its "basic"/reflected-quantile confidence
#' interval, and its reflected-empirical-CDF p-value -- the Monte-Carlo analog
#' of an analytic Cox-Snell (1968) first-order bias correction. Shared by every
#' \code{InferenceParamBootstrap} family that already implements
#' \code{simulate_under_lik_null()} for parametric-bootstrap LR calibration,
#' since all three public methods (defined directly on \code{InferenceParamBootstrap},
#' not this mixin -- see note below) reuse that same simulate-and-refit contract,
#' just anchored at the unrestricted fit instead of a null-restricted one.
#'
#' Public methods are deliberately \emph{not} spliced from this mixin:
#' roxygen2's R6 method-doc collector only associates description and argument
#' documentation comments with methods physically defined in the same file as
#' the \code{R6Class()} call, so a mixin-supplied public method silently loses
#' its documentation. This mirrors every other Pattern-1 mixin in this package
#' (e.g. \code{InferenceMixinBartlettApprox}), which are all private-only for
#' the same reason. Splice this mixin's private list into a daughter class (in
#' practice, just \code{InferenceParamBootstrap} itself, once) via
#' \code{private = c(InferenceMixinParamBootstrapEstimate$private, list(...))}.
#'
#' Depends on host-private \code{get_likelihood_test_spec()},
#' \code{effective_parallel_cores()}, \code{par_lapply()},
#' \code{use_deterministic_param_bootstrap()}, \code{with_param_bootstrap_seed()},
#' \code{simulate_under_lik_null()}, and \code{cache_nonestimable_se()}, all
#' defined on \code{InferenceParamBootstrap} itself or an ancestor class (not
#' part of this mixin), so this mixin is only meaningful spliced into that class
#' or a class providing the same contract.
#'
#' Capability flag: \code{private$supports_param_bootstrap_estimate()}.
#'
#' @keywords internal
#' @noRd
InferenceMixinParamBootstrapEstimate = list(
	public = list(),
	private = list(
		# Shared batch runner used by both compute_param_bootstrap_estimate()
		# and compute_param_bootstrap_confidence_interval(): extracts the raw
		# (unrestricted) estimate from the likelihood-test spec, runs B
		# simulate-and-refit replicates anchored at that same fit, and stores
		# diagnostics. Returns NULL (having already recorded a nonestimable
		# reason) if the spec or raw estimate is unavailable; otherwise a list
		# with raw_estimate, finite_reps, and n_success -- callers apply their
		# own min_number_usable_samples gate and final formula (mean-reflection
		# for the point estimate, quantile-reflection for the interval).
		run_param_bootstrap_estimate_batch = function(B, max_attempts_per_replicate, show_progress = FALSE){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec) || is.null(spec$full_fit) || is.null(spec$j)) {
				if (!isTRUE(self$is_nonestimable("estimate")))
					private$cache_nonestimable_se("param_bootstrap_estimate_spec_unavailable")
				return(NULL)
			}
			j = as.integer(spec$j)
			raw_estimate = as.numeric(spec$full_fit$b[j])
			if (!is.finite(raw_estimate)) {
				if (!isTRUE(self$is_nonestimable("estimate")))
					private$cache_nonestimable_se("param_bootstrap_estimate_raw_estimate_nonfinite")
				return(NULL)
			}

			run = private$run_param_bootstrap_estimate_replicates(
				spec = spec,
				full_fit = spec$full_fit,
				B = as.integer(B),
				max_attempts_per_replicate = max_attempts_per_replicate,
				show_progress = show_progress
			)
			replicate_estimates = vapply(run$results, function(res) as.numeric(res$b %||% NA_real_)[1L], numeric(1))
			finite_reps = replicate_estimates[is.finite(replicate_estimates)]
			n_success = length(finite_reps)
			n_failure = length(replicate_estimates) - n_success

			private$cached_values$last_param_bootstrap_estimate_diagnostics = list(
				B = as.integer(B),
				raw_estimate = raw_estimate,
				replicate_estimates = replicate_estimates,
				n_success = as.integer(n_success),
				n_failure = as.integer(n_failure),
				used_deterministic_mode = isTRUE(run$used_deterministic_mode)
			)

			list(raw_estimate = raw_estimate, finite_reps = finite_reps, n_success = n_success)
		},
		# Shared replicate-running core for compute_param_bootstrap_estimate():
		# simulates B datasets anchored at the unrestricted full_fit (not a
		# null-restricted fit), refits each unrestricted via the same
		# simulate_under_lik_null() contract used by the LR-bootstrap path, and
		# collects each replicate's fitted coefficient vector. Deliberately
		# simpler than run_param_bootstrap_replicates() -- no null refit is
		# needed here, and no reusable-worker-state duplication path, since a
		# point estimate's bias correction is not on the same CI-inversion hot
		# path (repeated at many delta values) that motivated that optimization
		# for the LR-bootstrap methods.
		run_param_bootstrap_estimate_replicates = function(spec, full_fit, B, max_attempts_per_replicate, show_progress = FALSE){
			B = as.integer(B)
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
			replicate_seeds = sample.int(.Machine$integer.max, B, replace = TRUE)
			deterministic_mode = isTRUE(private$use_deterministic_param_bootstrap())

			run_one = function(b){
				private$compute_param_bootstrap_estimate_impl(
					spec = spec,
					full_fit = full_fit,
					seed = replicate_seeds[[b]],
					max_attempts_per_replicate = max_attempts_per_replicate
				)
			}

			results = if (actual_cores <= 1L) {
				lapply(seq_len(B), run_one)
			} else {
				unlist(private$par_lapply(
					as.list(seq_len(B)),
					function(idx) list(run_one(as.integer(idx)[1L])),
					n_cores = actual_cores,
					budget = 1L,
					show_progress = show_progress
				), recursive = FALSE, use.names = FALSE)
			}

			list(results = results, used_deterministic_mode = deterministic_mode)
		},
		# One parametric-bootstrap estimate replicate: simulate a dataset
		# anchored at full_fit's own coefficients, refit unrestricted, and
		# return the refit's coefficient vector (not an LR statistic).
		compute_param_bootstrap_estimate_impl = function(spec, full_fit, seed = NULL, max_attempts_per_replicate = 1L){
			j = as.integer(spec$j)
			private$with_param_bootstrap_seed(seed, {
				last_result = list(success = FALSE, b = NA_real_, reason = "simulated_data_failure", attempts = 0L)
				for (attempt in seq_len(max(1L, as.integer(max_attempts_per_replicate)))) {
					boot_spec = tryCatch(
						private$simulate_under_lik_null(spec, delta = full_fit$b[j], null_fit = full_fit),
						error = function(e) NULL
					)
					if (is.null(boot_spec) || is.null(boot_spec$full_fit) ||
						length(boot_spec$full_fit$b) < j || !is.finite(boot_spec$full_fit$b[j])) {
						last_result = list(success = FALSE, b = NA_real_, reason = "simulated_refit_failure", attempts = as.integer(attempt))
						next
					}
					return(list(success = TRUE, b = as.numeric(boot_spec$full_fit$b[j]), reason = "success", attempts = as.integer(attempt)))
				}
				last_result
			})
		},
		# Parametric-bootstrap point-estimate bias correction automatically follows
		# parametric-bootstrap LR support: any family that already implements
		# simulate_under_lik_null() gets compute_param_bootstrap_estimate() for
		# free, since both only depend on that same simulate-and-refit contract.
		# Families wanting to withhold estimate bias correction for a reason that
		# doesn't affect LR bootstrap itself should override this directly.
		supports_param_bootstrap_estimate = function(){
			isTRUE(private$supports_lik_ratio_param_bootstrap())
		}
	)
)
