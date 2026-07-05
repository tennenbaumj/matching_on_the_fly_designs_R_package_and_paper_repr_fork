test_that("hurdle Poisson likelihood is stable on both log1mexp branches", {
	log_two <- log(2)
	lambda <- c(1e-8, 0.25, log_two * (1 - 1e-8),
		log_two * (1 + 1e-8), 2, 30)
	eta_cond <- log(lambda)
	eta_zi <- c(-2, -0.5, 0.25, 1, -1.5, 2)
	y <- c(1, 2, 1, 3, 0, 4)

	# The zero column leaves one harmless free parameter, allowing every
	# parameter that determines the likelihood to be held at a known value.
	X <- cbind(eta_cond, dummy = 0)
	Xzi <- cbind(intercept = 1, eta_zi)
	params <- c(1, 7, 0, 1)
	fit <- EDI:::fast_zero_augmented_poisson_cpp(
		X, y, Xzi, is_hurdle = TRUE,
		warm_start_params = params,
		fixed_idx = c(1L, 3L, 4L),
		fixed_values = params[c(1L, 3L, 4L)],
		estimate_only = TRUE
	)

	log_likelihood <- ifelse(
		y == 0,
		plogis(eta_zi, log.p = TRUE),
		plogis(eta_zi, log.p = TRUE, lower.tail = FALSE) +
			y * eta_cond - lambda - log(-expm1(-lambda))
	)

	expect_true(fit$converged)
	expect_equal(as.numeric(fit$params), params, tolerance = 0)
	expect_equal(fit$neg_loglik, -sum(log_likelihood), tolerance = 1e-12)
})
