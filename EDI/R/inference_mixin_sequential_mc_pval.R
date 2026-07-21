#' Mixin for Sequential Monte Carlo Early-Stopping in Randomization P-values
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' providing the two pieces of Sequential Monte Carlo (SMC) early-stopping logic
#' that are byte-identical between \code{InferenceRand}'s
#' \code{compute_two_sided_pval_with_sequential_mc()} and
#' \code{InferenceRandBootstrap}'s \code{compute_two_sided_brt_pval_with_sequential_mc()}:
#' the "is SMC configured and enabled" gate, and the "does the current
#' Clopper-Pearson-style band for the running p-value estimate exclude the
#' stopping threshold" decision.
#'
#' Deliberately does \emph{not} unify the two callers' outer batching loops --
#' those differ for real (the plain-randomization side grows its prefix via
#' \code{get_randomization_distribution_prefix()}'s own incremental cache, while
#' the BRT side tracks an explicit \code{target} and has an extra
#' full-distribution cache short-circuit before it ever enters the loop), so
#' forcing them into one shape would risk changing behavior rather than just
#' removing duplication.
#'
#' Splice into \code{InferenceRand} via
#' \code{private = c(InferenceMixinSequentialMCPval$private, list(...))}.
#' \code{InferenceRandBootstrap} does not need to splice this in separately --
#' it inherits \code{InferenceRand} (via \code{InferenceRandCI} /
#' \code{InferenceNonParamBootstrap}), so \code{private$sequential_mc_control_enabled()}
#' / \code{private$sequential_mc_band_excludes_threshold()} are already visible
#' there through ordinary R6 private-environment inheritance.
#'
#' @keywords internal
#' @noRd
InferenceMixinSequentialMCPval = list(
	public = list(),
	private = list(
		# TRUE iff `mc_ctrl` is a usable, enabled SMC control object (i.e. has
		# mc_enable = TRUE and a finite mc_stop_threshold). Both SMC call sites
		# return NULL immediately when this is FALSE.
		sequential_mc_control_enabled = function(mc_ctrl){
			!is.null(mc_ctrl) && isTRUE(mc_ctrl$mc_enable) && is.finite(mc_ctrl$mc_stop_threshold)
		},
		# TRUE iff the two-sided randomization-pval confidence band for the current
		# (possibly partial) null distribution `t0s` excludes `threshold`, i.e. we
		# are confident enough in the current running p-value estimate to stop
		# drawing more null replicates early.
		sequential_mc_band_excludes_threshold = function(t0s, t, threshold, conf_level){
			band = private$compute_two_sided_randomization_pval_band(t0s, t, conf_level)
			is.finite(band[1]) && is.finite(band[2]) && (band[2] < threshold || band[1] > threshold)
		}
	)
)
