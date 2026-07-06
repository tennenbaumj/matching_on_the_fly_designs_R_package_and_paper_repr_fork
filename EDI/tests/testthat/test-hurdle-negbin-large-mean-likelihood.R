test_that("truncated NB large-mean likelihood matches exact normalization", {
	theta <- 8
	p_zero <- c(1e-10, 1e-8, 5e-8, 2e-7, 1e-4, 0.1)
	mu <- theta * (p_zero^(-1 / theta) - 1)
	eta <- log(mu)
	y <- c(120L, 100L, 80L, 60L, 20L, 2L)

	# Leave a zero-valued dummy coefficient free while fixing every parameter
	# that determines the likelihood.
	X <- cbind(eta, dummy = 0)
	params <- c(1, 7, log(theta))
	fit <- EDI:::fast_truncated_negbin_count_cpp(
		X, y,
		warm_start_params = params,
		fixed_idx = c(1L, 3L),
		fixed_values = params[c(1L, 3L)],
		estimate_only = TRUE
	)

	exact_loglik <- sum(
		dnbinom(y, size = theta, mu = mu, log = TRUE) - log1p(-p_zero)
	)
	expect_true(fit$converged)
	expect_equal(as.numeric(fit$params), params, tolerance = 0)
	expect_equal(fit$neg_ll, -exact_loglik, tolerance = 1e-12)
})
