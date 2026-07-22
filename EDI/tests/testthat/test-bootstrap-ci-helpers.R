test_that("bootstrap CI helpers compute percentile and basic intervals", {
	boot_distr = c(-2, -1, 0, 2, 5, 9)
	alpha = 0.2
	est = 1.5

	percentile = EDI:::bootstrap_ci_from_distribution(boot_distr, alpha, "percentile")
	expect_equal(
		percentile,
		stats::quantile(boot_distr, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8)
	)
	expect_equal(
		EDI:::bootstrap_ci_from_distribution(boot_distr, alpha, "basic", est = est),
		2 * est - rev(percentile)
	)
	expect_error(
		EDI:::bootstrap_ci_from_distribution(boot_distr, alpha, "basic"),
		"require an estimate"
	)
})

test_that("bootstrap CI helpers compute stable studentized intervals", {
	boot_stats = list(
		theta = c(2.1, 1.8, 2.4, 2.0, 2.2),
		se = c(0.5, 0.5, 0.5, 0.5, 1e-10)
	)
	est = 2
	se_hat = 0.4
	alpha = 0.1
	pivots = EDI:::bootstrap_studentized_pivots(
		theta = boot_stats$theta,
		se = boot_stats$se,
		est = est,
		se_hat = se_hat,
		min_number_usable_samples = 4L
	)
	expect_equal(pivots, c(0.2, -0.4, 0.8, 0))

	studentized = EDI:::bootstrap_ci_studentized(
		boot_stats = boot_stats,
		alpha = alpha,
		est = est,
		se_hat = se_hat,
		min_number_usable_samples = 4L
	)
	q = stats::quantile(pivots, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
	expect_equal(studentized, c(est - q[1L] * se_hat, est - q[2L] * se_hat))

	symmetric = EDI:::bootstrap_ci_symmetric_studentized(
		boot_stats = boot_stats,
		alpha = alpha,
		est = est,
		se_hat = se_hat,
		min_number_usable_samples = 4L
	)
	expect_equal(sum(symmetric), 2 * est)
	expect_true(EDI:::bootstrap_studentized_interval_scale_unstable(
		theta = boot_stats$theta,
		ci = c(-100, 100),
		est = est
	))
})

test_that("bootstrap CI helpers reject unusable studentized standard errors", {
	expect_error(
		EDI:::bootstrap_studentized_pivots(
			theta = c(1, 2, 3),
			se = c(0, NA_real_, -1),
			est = 2,
			se_hat = 1
		),
		"finite positive standard errors"
	)
})
