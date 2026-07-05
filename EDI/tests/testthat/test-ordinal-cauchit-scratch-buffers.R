test_that("ordinal cauchit scratch-buffer score and Hessian match numerical derivatives", {
	skip_if_not_installed("numDeriv")
	set.seed(5101)
	X <- matrix(rnorm(90), ncol = 3)
	y <- rep(1:3, length.out = nrow(X))
	params <- c(-0.8, 1.1, 0.25, -0.4, 0.15)

	nll <- function(par) {
		alpha <- par[1:2]
		eta <- drop(X %*% par[3:5])
		upper <- ifelse(y == 3, 1, pcauchy(alpha[pmin(y, 2)] - eta))
		lower <- ifelse(y == 1, 0, pcauchy(alpha[pmax(y - 1, 1)] - eta))
		-sum(log(pmax(upper - lower, 1e-12)))
	}

	score <- EDI:::get_ordinal_cauchit_regression_score_cpp(X, y, params)
	hessian <- EDI:::get_ordinal_cauchit_regression_hessian_cpp(X, y, params)

	expect_equal(as.numeric(score), as.numeric(-numDeriv::grad(nll, params)), tolerance = 1e-7)
	expect_equal(unname(hessian), unname(-numDeriv::hessian(nll, params)), tolerance = 1e-5)
})
