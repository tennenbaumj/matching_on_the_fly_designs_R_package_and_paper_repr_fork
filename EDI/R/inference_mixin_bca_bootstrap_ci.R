#' Mixin for BCa (Bias-Corrected and Accelerated) Bootstrap Confidence Intervals and P-values
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' providing the generic Efron BCa correction math shared by
#' \code{InferenceNonParamBootstrap}'s plain-bootstrap \code{ci_bca()}/\code{pval_bca()}
#' and \code{InferenceBayesianBootstrap}'s \code{ci_bayesian_bca()}/\code{pval_bayesian_bca()}.
#' The two call sites differ only in where the jackknife replicate distribution comes
#' from and how a "can't compute this" failure is reported back to the caller; both
#' differences are captured by the \code{jack}/\code{reason_prefix}/\code{on_failure}
#' arguments below, so the bias-correction (z0) / acceleration (a) / adjusted-quantile
#' math itself lives in exactly one place.
#'
#' Splice into \code{InferenceNonParamBootstrap} via
#' \code{private = c(InferenceMixinBcaBootstrapCI$private, list(...))}.
#' \code{InferenceBayesianBootstrap} does not need to splice this in separately --
#' it inherits \code{InferenceNonParamBootstrap} (via \code{InferenceRandBootstrap} /
#' \code{InferenceRandBootstrapCI}), so \code{private$bca_ci_core()} /
#' \code{private$bca_pval_core()} are already visible there through ordinary R6
#' private-environment inheritance.
#'
#' @keywords internal
#' @noRd
InferenceMixinBcaBootstrapCI = list(
	public = list(),
	private = list(
		# Shared BCa confidence-interval math (Efron 1987; Efron & Tibshirani 1993).
		# `jack` must already be the finite jackknife replicate distribution (length >= 2);
		# callers are responsible for sourcing it (plain vs. Bayesian jackknife) and for
		# handling the "too few jackknife replicates" case themselves, since that failure
		# is reported differently at each call site.
		# `on_failure(reason)` is called with a bare reason suffix (e.g.
		# "unstable_bias_or_acceleration") -- the caller supplies `reason_prefix` so the
		# cached nonestimable-reason string matches its own family's existing vocabulary
		# ("bootstrap_bca_..." vs "bayesian_bootstrap_bca_...").
		bca_ci_core = function(boot_distr, alpha, est, jack, reason_prefix, on_failure){
			p_less = mean(boot_distr < est)
			p_less = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_less))
			z0 = stats::qnorm(p_less)
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			if (!is.finite(z0) || !is.finite(a) || abs(z0) > 2.5 || abs(a) > 1) {
				return(on_failure(paste0(reason_prefix, "unstable_bias_or_acceleration")))
			}
			alpha_vec = c(alpha / 2, 1 - alpha / 2)
			z_alpha = stats::qnorm(alpha_vec)
			denom = 1 - a * (z0 + z_alpha)
			if (any(!is.finite(denom)) || any(abs(denom) < sqrt(.Machine$double.eps))) {
				return(on_failure(paste0(reason_prefix, "unstable_bias_or_acceleration")))
			}
			adj = stats::pnorm(z0 + (z0 + z_alpha) / denom)
			prob_eps = 1 / (length(boot_distr) + 1)
			adj = pmin(1 - prob_eps, pmax(prob_eps, adj))
			adj = sort(adj)
			if (any(adj <= 2 * prob_eps) || any(adj >= 1 - 2 * prob_eps)) {
				return(on_failure(paste0(reason_prefix, "adjustment_on_boundary")))
			}
			if (diff(adj) < prob_eps) {
				return(stats::quantile(boot_distr, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8))
			}
			stats::quantile(boot_distr, probs = adj, names = FALSE, type = 8)
		},
		# Shared BCa two-sided p-value math -- see bca_ci_core() above for the shared
		# bias-correction (z0) / acceleration (a) derivation; this additionally maps a
		# tested null value `delta` through the same adjustment to get a two-sided p-value.
		bca_pval_core = function(boot_distr, est, delta, jack, reason_prefix, on_failure){
			p_less = mean(boot_distr < est)
			p_less = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_less))
			z0 = stats::qnorm(p_less)
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			if (!is.finite(z0) || !is.finite(a) || abs(z0) > 2.5 || abs(a) > 1) {
				return(on_failure(paste0(reason_prefix, "unstable_bias_or_acceleration")))
			}
			z_alpha = stats::qnorm(c(0.025, 0.975))
			denom_ci = 1 - a * (z0 + z_alpha)
			if (any(!is.finite(denom_ci)) || any(abs(denom_ci) < sqrt(.Machine$double.eps))) {
				return(on_failure(paste0(reason_prefix, "unstable_bias_or_acceleration")))
			}
			prob_eps = 1 / (length(boot_distr) + 1)
			adj_ci = sort(stats::pnorm(z0 + (z0 + z_alpha) / denom_ci))
			if (any(!is.finite(adj_ci)) || any(adj_ci <= 2 * prob_eps) || any(adj_ci >= 1 - 2 * prob_eps)) {
				return(on_failure(paste0(reason_prefix, "adjustment_on_boundary")))
			}
			p_delta = mean(boot_distr < delta)
			p_delta = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_delta))
			z_delta = stats::qnorm(p_delta)
			s = z_delta - z0
			denom = 1 + a * s
			if (!is.finite(denom) || abs(denom) < sqrt(.Machine$double.eps)) {
				return(on_failure(paste0(reason_prefix, "unstable_bias_or_acceleration")))
			}
			adj_z = s / denom - z0
			if (!is.finite(adj_z) || abs(adj_z) > 8) {
				return(on_failure(paste0(reason_prefix, "adjustment_on_boundary")))
			}
			p_raw = min(1, 2 * min(stats::pnorm(adj_z), 1 - stats::pnorm(adj_z)))
			min(1, max(2 / length(boot_distr), p_raw))
		}
	)
)
