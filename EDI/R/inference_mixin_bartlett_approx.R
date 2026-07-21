#' Mixin for Approximate (Monte-Carlo) Bartlett-Corrected Likelihood-Ratio Inference
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' providing the generic Monte-Carlo Bartlett correction factor for the LR test,
#' shared by every \code{InferenceParamBootstrap} family that already implements
#' \code{simulate_under_lik_null()} for parametric-bootstrap LR calibration.
#' Splice into a daughter class (in practice, just \code{InferenceParamBootstrap}
#' itself, once) via
#' \code{public = c(InferenceMixinBartlettApprox$public, list(...))} and
#' \code{private = c(InferenceMixinBartlettApprox$private, list(...))}.
#'
#' Depends on \code{private$run_param_bootstrap_replicates(...)} and
#' \code{private$param_bootstrap_lr_extreme(...)}, both defined on
#' \code{InferenceParamBootstrap} itself (not part of this mixin), so this mixin
#' is only meaningful spliced into that class or a class providing the same
#' contract.
#'
#' Capability flag: \code{private$supports_bartlett_likelihood_ratio_approx()}.
#'
#' @keywords internal
#' @noRd
InferenceMixinBartlettApprox = list(
	public = list(),
	private = list(
		bartlett_factor_mc_min_usable_fraction = 0.2,
		bartlett_factor_mc_max_attempts_per_replicate = 2L,
		# Approximate (Monte-Carlo) Bartlett support automatically follows
		# parametric-bootstrap LR support: any family that already implements
		# simulate_under_lik_null() gets the shared Monte-Carlo factor below for
		# free. Families that need to withhold Bartlett support for a reason that
		# doesn't affect parametric-bootstrap LR itself (e.g. a known raw-LR
		# miscalibration) should override this method directly.
		supports_bartlett_likelihood_ratio_approx = function(){
			isTRUE(private$supports_lik_ratio_param_bootstrap())
		},
		# Generic Monte-Carlo Bartlett correction factor, shared by every
		# InferenceParamBootstrap family that already implements
		# simulate_under_lik_null() for parametric-bootstrap LR calibration.
		#
		# Simulates B datasets under H0: (tested parameter) = delta using the same
		# simulate_under_lik_null()/refit machinery as
		# compute_lik_ratio_bootstrap_two_sided_pval() -- including the same
		# multi-core parallelism and reusable-worker-state optimizations, via the
		# shared private$run_param_bootstrap_replicates() helper -- and sets
		#   c(delta) = mean(simulated LR statistics)
		# which is a Monte-Carlo estimate of E[LR | H0], the quantity a classical
		# analytic Bartlett correction targets exactly (E[LR]/df = 1 under a
		# chi-square(df) reference). This is deliberately model-agnostic: any
		# family wanting a faster closed-form factor instead can override
		# get_bartlett_factor_exact()/supports_bartlett_likelihood_ratio_exact()
		# (a subclass override of this method wins too, if preferred).
		get_bartlett_factor_approx = function(spec, delta, full_fit, null_fit, B = 99){
			if (!isTRUE(private$supports_lik_ratio_param_bootstrap())) return(NULL)
			if (is.null(null_fit)) return(NULL)

			B = as.integer(B)
			# Scales with B rather than a fixed count, so small-B calls (e.g. quick
			# smoke tests) aren't structurally impossible to satisfy.
			min_usable = max(2L, as.integer(ceiling(private$bartlett_factor_mc_min_usable_fraction * B)))
			max_attempts = as.integer(private$bartlett_factor_mc_max_attempts_per_replicate)

			private$active_resampling_operation = "param_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)

			run = tryCatch(
				private$run_param_bootstrap_replicates(
					spec = spec,
					delta = delta,
					null_fit = null_fit,
					B = B,
					max_attempts_per_replicate = max_attempts,
					show_progress = FALSE,
					# Empirically slower here, not faster: this is called at many
					# different delta values (CI root-finding), so the worker-state's
					# one-time duplication cost never amortizes across delta the way it
					# does for a single-delta call, and its per-replicate cost was
					# observed to be ~5x the plain path's for InferenceIncidLogRegr.
					allow_worker_reuse = FALSE
				),
				error = function(e) NULL
			)
			if (is.null(run)) return(NULL)

			lr_boots = vapply(run$results, function(res) as.numeric(res$lr %||% NA_real_)[1L], numeric(1))
			extreme = is.finite(lr_boots) & private$param_bootstrap_lr_extreme(lr_boots)
			lr_boots[extreme] = NA_real_
			finite_lr = lr_boots[is.finite(lr_boots)]
			if (length(finite_lr) < min_usable) return(NULL)

			factor = mean(finite_lr)
			if (!is.finite(factor) || factor <= 0) return(NULL)
			factor
		}
	)
)
