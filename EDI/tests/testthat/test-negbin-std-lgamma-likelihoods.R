test_that("standard and truncated NB std::lgamma likelihoods match R", {
	mu <- c(0.1, 0.4, 1, 2.5, 8, 25)
	eta <- log(mu)
	theta <- 2.3

	# A zero predictor column leaves one harmless free parameter while all
	# likelihood-determining parameters are fixed at known values.
	X <- cbind(eta, dummy = 0)
	params <- c(1, 7, log(theta))
	fixed_idx <- c(1L, 3L)

	y_nb <- c(0L, 1L, 2L, 5L, 15L, 40L)
	fit_nb <- EDI:::fast_neg_bin_cpp(
		X, y_nb,
		warm_start_params = params,
		fixed_idx = fixed_idx,
		fixed_values = params[fixed_idx],
		estimate_only = TRUE
	)
	loglik_nb <- sum(dnbinom(y_nb, size = theta, mu = mu, log = TRUE))

	y_truncated <- c(1L, 2L, 3L, 6L, 15L, 40L)
	fit_truncated <- EDI:::fast_truncated_negbin_count_cpp(
		X, y_truncated,
		warm_start_params = params,
		fixed_idx = fixed_idx,
		fixed_values = params[fixed_idx],
		estimate_only = TRUE
	)
	log_p_zero <- dnbinom(0, size = theta, mu = mu, log = TRUE)
	loglik_truncated <- sum(
		dnbinom(y_truncated, size = theta, mu = mu, log = TRUE) -
			log1p(-exp(log_p_zero))
	)

	expect_true(fit_nb$converged)
	expect_true(fit_truncated$converged)
	expect_equal(as.numeric(fit_nb$b), params[1:2], tolerance = 0)
	expect_equal(as.numeric(fit_truncated$params), params, tolerance = 0)
	expect_equal(fit_nb$logLik, loglik_nb, tolerance = 1e-12)
	expect_equal(fit_truncated$neg_ll, -loglik_truncated, tolerance = 1e-12)
})
