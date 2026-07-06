library(testthat)
library(EDI)

kk21_negbin_weight_reference <- function(x, y) {
	design <- cbind(1, x)
	beta <- numeric(2L)
	y_mean <- mean(y)
	y_var <- sum((y - y_mean)^2) / (length(y) - 1)
	theta <- if (y_var > y_mean) y_mean^2 / (y_var - y_mean) else 10
	theta <- pmax(0.01, pmin(1000, theta))

	for (outer in seq_len(50L)) {
		for (iteration in seq_len(25L)) {
			eta <- drop(design %*% beta)
			mu <- exp(pmin(eta, 20))
			w <- pmax(1e-10, mu / (1 + mu / theta))
			z <- eta + (y - mu) / mu
			beta_new <- drop(solve(crossprod(design, design * w), crossprod(design, w * z)))
			if (sqrt(sum((beta - beta_new)^2)) < 1e-6) {
				beta <- beta_new
				break
			}
			beta <- beta_new
		}

		eta <- drop(design %*% beta)
		mu <- exp(pmin(eta, 20))
		score <- sum(
			digamma(y + theta) - digamma(theta) + log(theta) - log(theta + mu) + 1 -
				(y + theta) / (theta + mu)
		)
		information <- sum(
			-trigamma(y + theta) + trigamma(theta) - 1 / theta + 2 / (theta + mu) -
				(y + theta) / (theta + mu)^2
		)
		if (abs(information) < 1e-10) break
		theta_new <- pmax(0.01, pmin(1000, theta - score / information))
		if (abs(theta_new - theta) < 1e-6) {
			theta <- theta_new
			break
		}
		theta <- theta_new
	}

	eta <- drop(design %*% beta)
	mu <- exp(pmin(eta, 20))
	w <- pmax(1e-10, mu / (1 + mu / theta))
	variance <- solve(crossprod(design, design * w))[2L, 2L]
	abs(beta[2L] / sqrt(variance))
}

test_that("KK21 negative-binomial weights match an independent R implementation", {
	set.seed(78)
	n <- 80L
	X <- cbind(x1 = rnorm(n), x2 = runif(n, -1, 1), x3 = rep(c(-1, 1), n / 2L))
	mu <- exp(0.3 + 0.45 * X[, 1L] - 0.25 * X[, 2L])
	y <- as.numeric(rnbinom(n, size = 2.5, mu = mu))

	expected <- apply(X, 2L, kk21_negbin_weight_reference, y = y)
	actual <- EDI:::kk21_negbin_weights_cpp(X, y)
	expect_equal(as.numeric(actual), unname(expected), tolerance = 1e-7)
})

test_that("KK21 negative-binomial weights preserve edge and stepwise invariants", {
	expect_identical(
		as.numeric(EDI:::kk21_negbin_weights_cpp(matrix(1:4, 2L), c(0, 1))),
		rep(.Machine$double.eps, 2L)
	)

	set.seed(781)
	X <- matrix(rnorm(120L), 40L, 3L)
	y <- as.numeric(rnbinom(40L, size = 2, mu = exp(0.2 + 0.3 * X[, 1L])))
	weights <- EDI:::kk21_stepwise_negbin_weights_cpp(X, y, rep(1, 40L))
	expect_length(weights, ncol(X))
	expect_true(all(is.finite(weights)))
	expect_true(all(weights >= 0))
})
