library(testthat)
library(EDI)

test_that("ordinal information outputs and variance use the same Hessian", {
	set.seed(31)
	n <- 240L
	X <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L)
	latent <- drop(X %*% c(0.5, -0.3, 0.2)) + stats::rlogis(n)
	y <- as.numeric(cut(latent, breaks = c(-Inf, -1, 0, 1, Inf)))

	fit <- fast_ordinal_regression_cpp(X, y)
	expect_identical(fit$observed_information, fit$fisher_information)
	expect_identical(fit$observed_information, fit$information)

	fit_var <- fast_ordinal_regression_with_var_cpp(X, y)
	expect_identical(fit_var$fisher_information, fit$observed_information)
	expect_equal(
		unname(fit_var$vcov),
		unname(solve(fit_var$fisher_information)),
		tolerance = 1e-12
	)

	weighted_fit <- fast_ordinal_regression_weighted_cpp(
		X,
		y,
		weights = seq(0.5, 1.5, length.out = n)
	)
	expect_identical(
		weighted_fit$observed_information,
		weighted_fit$fisher_information
	)
	expect_identical(weighted_fit$observed_information, weighted_fit$information)
})
