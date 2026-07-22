# Internal pure helpers for bootstrap confidence-interval calculations.
# These functions intentionally do not access inference objects, caches, or
# resampling state so callers can test and reuse the numerical operations directly.

bootstrap_ci_from_distribution = function(boot_distr, alpha, type, est = NULL){
	type = tolower(type)
	if (identical(type, "percentile")) {
		return(stats::quantile(
			boot_distr,
			probs = c(alpha / 2, 1 - alpha / 2),
			names = FALSE,
			type = 8
		))
	}
	if (is.null(est)) stop("Basic bootstrap confidence intervals require an estimate.")
	2 * est - stats::quantile(
		boot_distr,
		probs = c(1 - alpha / 2, alpha / 2),
		names = FALSE,
		type = 8
	)
}

bootstrap_studentized_interval_scale_unstable = function(theta, ci = NULL, se_hat = NULL, pivots = NULL, est = 0, alpha = 0.05, max_width_ratio = 5){
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
}

bootstrap_studentized_pivots = function(theta, se, est, se_hat, min_number_usable_samples = 10L, symmetric = FALSE){
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
}

bootstrap_ci_studentized = function(boot_stats, alpha, est, se_hat, min_number_usable_samples = 10L){
	pivots = bootstrap_studentized_pivots(
		theta = boot_stats$theta,
		se = boot_stats$se,
		est = est,
		se_hat = se_hat,
		min_number_usable_samples = min_number_usable_samples
	)
	q = stats::quantile(pivots, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
	c(est - q[1L] * se_hat, est - q[2L] * se_hat)
}

bootstrap_ci_symmetric_studentized = function(boot_stats, alpha, est, se_hat, min_number_usable_samples = 10L){
	pivots = bootstrap_studentized_pivots(
		theta = boot_stats$theta,
		se = boot_stats$se,
		est = est,
		se_hat = se_hat,
		min_number_usable_samples = min_number_usable_samples,
		symmetric = TRUE
	)
	q = stats::quantile(pivots, probs = 1 - alpha, names = FALSE, type = 8)
	c(est - q[1L] * se_hat, est + q[1L] * se_hat)
}
