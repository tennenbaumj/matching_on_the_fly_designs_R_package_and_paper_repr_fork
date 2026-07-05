library(testthat)
library(EDI)

test_that("Poisson IRLS information matches direct weighted crossproducts", {
	set.seed(37)
	n <- 300L
	X <- cbind(1, matrix(rnorm(n * 3L), ncol = 3L))
	beta <- c(0.25, -0.2, 0.15, 0.1)
	y <- stats::rpois(n, exp(drop(X %*% beta)))

	fit <- fast_poisson_regression_cpp(X, y, optimization_alg = "irls")
	information_reference <- crossprod(X, X * as.numeric(fit$mu))
	expect_equal(unname(fit$XtWX), unname(information_reference), tolerance = 1e-10)
	expect_equal(
		as.numeric(fit$b),
		as.numeric(stats::glm.fit(X, y, family = stats::poisson())$coefficients),
		tolerance = 1e-7
	)

	weights <- seq(0.5, 1.5, length.out = n)
	weighted_fit <- fast_poisson_regression_weighted_cpp(
		X,
		y,
		weights,
		optimization_alg = "irls"
	)
	weighted_information_reference <- crossprod(
		X,
		X * (weights * as.numeric(weighted_fit$mu))
	)
	expect_equal(
		unname(weighted_fit$XtWX),
		unname(weighted_information_reference),
		tolerance = 1e-10
	)
	expect_equal(
		as.numeric(weighted_fit$b),
		as.numeric(stats::glm.fit(X, y, weights = weights, family = stats::poisson())$coefficients),
		tolerance = 1e-7
	)
})
