library(testthat)
library(EDI)

test_that("Weibull AFT reusable buffers preserve score and Hessian", {
	skip_if_not_installed("numDeriv")
	set.seed(48)
	n <- 17L
	X <- cbind(1, rnorm(n), runif(n, -1, 1))
	y <- exp(0.3 + 0.2 * X[, 2L] - 0.1 * X[, 3L] + rnorm(n, sd = 0.4))
	dead <- as.numeric(c(rep(1, 12L), rep(0, 5L)))
	params <- c(0.25, 0.15, -0.08, log(0.7))

	nll <- function(par) {
		eta <- drop(X %*% par[seq_len(ncol(X))])
		log_sigma <- par[ncol(X) + 1L]
		sigma <- exp(log_sigma)
		w <- pmin((log(y) - eta) / sigma, 700)
		-sum(dead * (w - log_sigma - log(y)) - exp(w))
	}

	score <- EDI:::get_weibull_regression_score_cpp(X, y, dead, params)
	hessian <- EDI:::get_weibull_regression_hessian_cpp(X, y, dead, params)
	expect_equal(
		as.numeric(score),
		as.numeric(-numDeriv::grad(nll, params)),
		tolerance = 1e-7
	)
	expect_equal(
		unname(hessian),
		unname(-numDeriv::hessian(nll, params)),
		tolerance = 1e-5
	)
})
