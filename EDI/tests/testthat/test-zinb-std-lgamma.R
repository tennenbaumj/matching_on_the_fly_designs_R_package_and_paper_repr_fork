test_that("ZINB std::lgamma likelihood matches an independent R calculation", {
	mu <- c(0.05, 0.5, 1, 3, 10, 40)
	eta_cond <- log(mu)
	eta_zi <- c(-2, -0.5, 0.25, 1, -1.5, 2)
	y <- c(0, 1, 2, 5, 20, 50)
	theta <- 2.3

	# The zero column leaves one harmless free parameter, allowing every
	# parameter that determines the likelihood to be held at a known value.
	X <- cbind(eta_cond, dummy = 0)
	Xzi <- cbind(intercept = 1, eta_zi)
	params <- c(1, 7, 0, 1, log(theta))
	fixed_idx <- c(1L, 3L, 4L, 5L)
	fit <- EDI:::fast_zinb_cpp(
		X, Xzi, y,
		warm_start_params = params,
		fixed_idx = fixed_idx,
		fixed_values = params[fixed_idx],
		estimate_only = FALSE
	)

	p_zero <- plogis(eta_zi)
	log_nb <- lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
		theta * (log(theta) - log(theta + mu)) +
		y * (eta_cond - log(theta + mu))
	log_likelihood <- ifelse(
		y == 0,
		log(p_zero + (1 - p_zero) * exp(log_nb)),
		plogis(eta_zi, lower.tail = FALSE, log.p = TRUE) + log_nb
	)

	expect_true(fit$converged)
	expect_equal(as.numeric(fit$params), params, tolerance = 0)
	expect_equal(fit$neg_loglik, -sum(log_likelihood), tolerance = 1e-12)
})
