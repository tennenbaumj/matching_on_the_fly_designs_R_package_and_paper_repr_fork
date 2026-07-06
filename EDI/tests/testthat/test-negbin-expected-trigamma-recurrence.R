test_that("negbin expected Hessian trigamma recurrence matches direct summation", {
	set.seed(6401)
	n <- 18L
	X <- cbind(1, matrix(rnorm(n * 2L), ncol = 2L))
	y <- rep(c(0L, 1L, 3L), length.out = n)
	params <- c(0.15, -0.25, 0.2, log(2.4))
	theta <- exp(params[4])
	mu <- exp(drop(X %*% params[1:3]))

	expected_trigamma_direct <- function(mu_i) {
		prob_success <- theta / (theta + mu_i)
		pk <- exp(theta * log(prob_success))
		value <- pk * trigamma(theta)
		cdf <- pk
		ratio_base <- mu_i / (theta + mu_i)
		min_iter <- ceiling(mu_i + 10 * sqrt(mu_i + mu_i^2 / theta))
		for (k in 0:99999) {
			pk <- pk * (k + theta) / (k + 1) * ratio_base
			y_k <- k + 1
			value <- value + pk * trigamma(y_k + theta)
			cdf <- cdf + pk
			if (y_k > min_iter && pk < 1e-14 && 1 - cdf < 1e-12) break
		}
		value
	}

	weights <- mu * theta / (theta + mu)
	reference <- matrix(0, ncol(X) + 1L, ncol(X) + 1L)
	reference[seq_len(ncol(X)), seq_len(ncol(X))] <- crossprod(X, weights * X)
	reference[ncol(reference), ncol(reference)] <- sum(
		-theta^2 * (vapply(mu, expected_trigamma_direct, numeric(1)) -
			trigamma(theta) + 1 / theta - 1 / (theta + mu))
	)

	actual <- EDI:::get_negbin_regression_expected_hessian_cpp(X, y, params)
	expect_equal(unname(actual), unname(reference), tolerance = 1e-11)
})
