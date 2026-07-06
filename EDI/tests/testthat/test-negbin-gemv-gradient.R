library(testthat)
library(EDI)

test_that("negbin GEMV gradient matches numerical derivative of log-likelihood", {
	skip_if_not_installed("numDeriv")
	set.seed(1701)
	n <- 120L; p <- 4L
	X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
	y <- as.integer(rnbinom(n, size = 3, mu = exp(X %*% c(0.5, -0.3, 0.2, 0.1))))
	params <- c(0.4, -0.2, 0.15, 0.08, log(2.5))  # [beta, log_theta]

	neg_ll_r <- function(par) {
		beta  <- par[seq_len(p)]
		theta <- exp(par[p + 1L])
		mu    <- exp(as.numeric(X %*% beta))
		-sum(dnbinom(y, size = theta, mu = mu, log = TRUE))
	}

	score_cpp <- as.numeric(EDI:::get_negbin_regression_score_cpp(X, y, params))
	score_num <- as.numeric(numDeriv::grad(neg_ll_r, params)) * -1  # gradient of log-lik

	expect_equal(score_cpp, score_num, tolerance = 1e-6)
})

test_that("negbin GEMV fit agrees with MASS::glm.nb coefficients", {
	skip_if_not_installed("MASS")
	set.seed(7777)
	n <- 200L; p <- 3L
	X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
	y <- as.integer(rnbinom(n, size = 2, mu = exp(X %*% c(0.5, 0.3, -0.2))))
	fit_cpp <- EDI:::fast_neg_bin_with_var_cpp(X, y)
	df <- data.frame(y = y, x1 = X[, 2L], x2 = X[, 3L])
	fit_r <- MASS::glm.nb(y ~ x1 + x2, data = df)
	expect_equal(as.numeric(fit_cpp$b), as.numeric(coef(fit_r)), tolerance = 1e-3)
	expect_equal(fit_cpp$theta, fit_r$theta, tolerance = 1e-2)
})
