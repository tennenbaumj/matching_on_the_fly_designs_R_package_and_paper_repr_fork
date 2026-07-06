library(testthat)
library(EDI)

kk21_beta_weight_reference <- function(x, y) {
	design <- cbind(1, x)
	beta <- numeric(2L)
	for (iteration in seq_len(100L)) {
		eta <- drop(design %*% beta)
		prob <- pmax(1e-10, pmin(1 - 1e-10, plogis(eta)))
		w <- pmax(1e-10, prob * (1 - prob))
		z <- eta + (y - prob) / w
		beta_new <- drop(solve(crossprod(design, design * w), crossprod(design, w * z)))
		if (sqrt(sum((beta - beta_new)^2)) < 1e-8) {
			beta <- beta_new
			break
		}
		beta <- beta_new
	}

	eta <- drop(design %*% beta)
	prob <- pmax(1e-12, pmin(1 - 1e-12, plogis(eta)))
	w <- pmax(1e-10, prob * (1 - prob))
	log_y <- log(y)
	log_one_minus_y <- log1p(-y)
	grid <- seq(-2, 7, by = 0.5)
	loglik <- vapply(grid, function(log_phi) {
		phi <- exp(log_phi)
		mu_phi <- pmax(1e-12, prob * phi)
		one_minus_mu_phi <- pmax(1e-12, (1 - prob) * phi)
		sum(
			lgamma(phi) - lgamma(mu_phi) - lgamma(one_minus_mu_phi) +
				(mu_phi - 1) * log_y +
				(one_minus_mu_phi - 1) * log_one_minus_y
		)
	}, numeric(1L))
	best_phi <- exp(grid[which.max(loglik)])
	information <- crossprod(design, design * ((1 + best_phi) * w))
	variance <- solve(information)[2L, 2L]
	abs(beta[2L] / sqrt(variance))
}

test_that("KK21 beta weights match an independent R implementation", {
	set.seed(77)
	n <- 80L
	X <- cbind(x1 = rnorm(n), x2 = runif(n, -1, 1), x3 = rep(c(-1, 1), n / 2L))
	mu <- plogis(0.2 + 0.7 * X[, 1L] - 0.4 * X[, 2L])
	phi <- 12
	y <- pmax(1e-6, pmin(1 - 1e-6, rbeta(n, mu * phi, (1 - mu) * phi)))

	expected <- apply(X, 2L, kk21_beta_weight_reference, y = y)
	actual <- EDI:::kk21_beta_weights_cpp(X, y)
	expect_equal(as.numeric(actual), unname(expected), tolerance = 1e-8)
})

test_that("KK21 beta weights preserve edge-case and stepwise invariants", {
	expect_identical(
		as.numeric(EDI:::kk21_beta_weights_cpp(matrix(1:4, 2L), c(0.2, 0.8))),
		rep(.Machine$double.eps, 2L)
	)

	set.seed(771)
	X <- matrix(rnorm(120L), 40L, 3L)
	y <- stats::rbeta(40L, 2, 3)
	weights <- EDI:::kk21_stepwise_beta_weights_cpp(X, y, rep(1, 40L))
	expect_length(weights, ncol(X))
	expect_true(all(is.finite(weights)))
	expect_true(all(weights >= 0))
})
