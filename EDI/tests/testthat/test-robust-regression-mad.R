library(testthat)
library(EDI)

test_that("robust regression MAD scale selects odd and even medians", {
	beta <- c(0.5, -0.25)

	check_scale <- function(deviations) {
		x <- seq_along(deviations) - 1L
		X <- cbind(1, x)
		y <- drop(X %*% beta) + deviations
		fit <- EDI:::fast_robust_regression_cpp(
			X,
			y,
			warm_start_beta = beta,
			maxit = 0L,
			estimate_only = TRUE
		)
		expect_equal(fit$scale, median(abs(deviations)) / 0.6745, tolerance = 1e-14)
	}

	check_scale(c(-9, -5, -3, -1, 2, 4, 8))
	check_scale(c(-9, -5, -3, -1, 2, 4, 6, 8))
	check_scale(rep(c(-9, -5, -3, -1, 2, 4, 8), length.out = 513L))
	check_scale(rep(c(-9, -5, -3, -1, 2, 4, 6, 8), length.out = 514L))
})
