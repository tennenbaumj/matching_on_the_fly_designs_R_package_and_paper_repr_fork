test_that("negbin distinct-y Hessian tables match numerical derivatives", {
	skip_if_not_installed("numDeriv")
	set.seed(6301)
	n <- 75L
	X <- cbind(1, matrix(rnorm(n * 2L), ncol = 2L))
	y <- rep(c(0L, 1L, 2L, 5L, 10L), length.out = n)
	params <- c(0.2, -0.15, 0.3, log(2.4))

	nll <- function(par) {
		mu <- exp(drop(X %*% par[1:3]))
		-sum(dnbinom(y, mu = mu, size = exp(par[4]), log = TRUE))
	}

	hessian <- EDI:::get_negbin_regression_hessian_cpp(X, y, params)
	expect_equal(
		unname(hessian),
		unname(-numDeriv::hessian(nll, params)),
		tolerance = 1e-5
	)
})
