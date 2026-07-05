library(testthat)
library(EDI)

test_that("logistic score and Hessian handle extreme linear predictors", {
	x <- c(-1000, -50, -5, 0, 5, 50, 1000)
	X <- cbind(1, x)
	beta <- c(0, 1)
	y <- c(0, 1, 0, 1, 0, 1, 1)
	weights <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5)
	mu <- plogis(x)
	working_weights <- mu * (1 - mu)

	expect_equal(
		as.numeric(EDI:::get_logistic_regression_score_cpp(X, y, beta)),
		as.numeric(crossprod(X, y - mu)),
		tolerance = 1e-12
	)
	expect_equal(
		unname(EDI:::get_logistic_regression_hessian_cpp(X, beta)),
		unname(-crossprod(X, X * working_weights)),
		tolerance = 1e-12
	)
	expect_equal(
		as.numeric(EDI:::get_logistic_regression_weighted_score_cpp(X, y, weights, beta)),
		as.numeric(crossprod(X, weights * (y - mu))),
		tolerance = 1e-12
	)
	expect_equal(
		unname(EDI:::get_logistic_regression_weighted_hessian_cpp(X, weights, beta)),
		unname(-crossprod(X, X * (weights * working_weights))),
		tolerance = 1e-12
	)
})
