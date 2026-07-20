library(testthat)
library(EDI)

kk21_logistic_reduced_fit_reference <- function(X, y, beta) {
	for (iteration in seq_len(100L)) {
		eta <- drop(X %*% beta)
		prob <- plogis(eta)
		weights <- pmax(1e-10, prob * (1 - prob))
		z <- eta + (y - prob) / weights
		information <- crossprod(X, X * weights)
		beta_new <- drop(solve(information, crossprod(X, weights * z)))
		if (sqrt(sum((beta - beta_new)^2)) < 1e-8) {
			beta <- beta_new
			break
		}
		beta <- beta_new
	}
	eta <- drop(X %*% beta)
	prob <- plogis(eta)
	weights <- pmax(1e-10, prob * (1 - prob))
	list(
		beta = beta,
		weights = weights,
		residuals = y - prob,
		information = crossprod(X, X * weights)
	)
}

kk21_stepwise_logistic_reference <- function(X, y, treatment) {
	p <- ncol(X)
	selected <- integer()
	unused <- seq_len(p)
	weights_out <- rep(NA_real_, p)
	model <- cbind(1, treatment)
	beta <- numeric(2L)

	for (step in seq_len(p)) {
		fit <- kk21_logistic_reduced_fit_reference(model, y, beta)
		stats <- vapply(unused, function(j) {
			x <- X[, j]
			score <- sum(x * fit$residuals)
			v <- drop(crossprod(model, fit$weights * x))
			adjusted_information <- sum(fit$weights * x^2) -
				drop(crossprod(v, solve(fit$information, v)))
			if (adjusted_information <= 0) return(0)
			abs(score) / sqrt(adjusted_information)
		}, numeric(1L))
		best_position <- which.max(stats)
		best <- unused[best_position]
		weights_out[best] <- stats[best_position]
		selected <- c(selected, best)
		unused <- unused[-best_position]
		model <- cbind(model, X[, best])
		beta <- c(fit$beta, 0)
	}
	weights_out
}

kk21_weibull_weight_reference <- function(x, time, event) {
	log_time <- log(pmax(time, 1e-10))
	design <- cbind(1, x)
	beta <- drop(solve(crossprod(design), crossprod(design, log_time)))
	residuals <- log_time - drop(design %*% beta)
	scale <- sqrt(sum(residuals[event > 0.5]^2) / max(1, sum(event) - 2))
	scale <- pmax(0.01, pmin(10, scale))

	for (iteration in seq_len(30L)) {
		z <- pmax(-20, pmin(20, residuals / scale))
		weights <- exp(z)
		adjustment <- ifelse(event > 0.5, z - 1 + weights, weights)
		pseudo_response <- log_time - scale * adjustment / weights
		beta_new <- drop(solve(
			crossprod(design, design * weights),
			crossprod(design, weights * pseudo_response)
		))
		residuals_new <- log_time - drop(design %*% beta_new)
		z_new <- pmax(-20, pmin(20, residuals_new / scale))
		score <- sum(event * (z_new - 1) - z_new * exp(z_new)) / scale
		information <- sum(event) / scale^2
		scale_new <- if (abs(information) > 1e-10) scale - score / information else scale
		scale_new <- pmax(0.01, pmin(10, scale_new))
		difference <- sqrt(sum((beta_new - beta)^2)) + abs(scale_new - scale)
		beta <- beta_new
		scale <- scale_new
		residuals <- residuals_new
		if (difference < 1e-5) break
	}

	weights <- exp(pmax(-20, pmin(20, residuals / scale))) / scale^2
	variance <- solve(crossprod(design, design * weights))[2L, 2L]
	abs(beta[2L] / sqrt(variance))
}

test_that("KK21 stepwise logistic weights match an R score-test reference", {
	set.seed(79)
	n <- 100L
	X <- cbind(x1 = rnorm(n), x2 = runif(n, -1, 1), x3 = rep(c(-1, 1, 1, -1), length.out = n))
	treatment <- rep(c(0, 1), n / 2L)
	y <- rbinom(n, 1, plogis(-0.2 + 0.5 * treatment + 0.7 * X[, 1L]))

	expected <- kk21_stepwise_logistic_reference(X, y, treatment)
	actual <- EDI:::kk21_stepwise_logistic_weights_cpp(X, y, treatment)
	expect_equal(as.numeric(actual), expected, tolerance = 1e-8)
})

test_that("KK21 survival weights match an R Weibull reference", {
	set.seed(791)
	n <- 100L
	X <- cbind(x1 = rnorm(n), x2 = runif(n, -1, 1), x3 = rep(c(-1, 1), n / 2L))
	time <- rexp(n, rate = exp(0.3 * X[, 1L] - 0.2 * X[, 2L]))
	event <- rbinom(n, 1, 0.8)

	expected <- apply(X, 2L, kk21_weibull_weight_reference, time = time, event = event)
	actual <- EDI:::kk21_survival_weights_cpp(X, time, event)
	expect_equal(as.numeric(actual), unname(expected), tolerance = 1e-7)
})
