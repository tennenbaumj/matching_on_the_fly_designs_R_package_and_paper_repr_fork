test_that("Poisson IRLS weighted_crossprod matches glm() on standard data", {
	set.seed(1)
	n <- 200
	X <- cbind(1, matrix(rnorm(n * 3), n, 3))
	beta_true <- c(0.5, -0.3, 0.2, 0.1)
	mu <- exp(X %*% beta_true)
	y <- rpois(n, mu)
	fit <- EDI:::fast_poisson_regression_cpp(X, y)
	ref <- glm.fit(X, y, family = poisson())
	expect_equal(as.numeric(fit$b), unname(ref$coefficients), tolerance = 1e-6)
})

test_that("Poisson IRLS weighted_crossprod matches glm() with weights", {
	set.seed(2)
	n <- 150
	X <- cbind(1, rnorm(n), rnorm(n))
	beta_true <- c(1, 0.4, -0.5)
	y <- rpois(n, exp(X %*% beta_true))
	w <- runif(n, 0.5, 2)
	fit <- EDI:::fast_poisson_regression_weighted_cpp(X, y, w)
	ref <- glm.fit(X, y, weights = w, family = poisson())
	expect_equal(as.numeric(fit$b), unname(ref$coefficients), tolerance = 1e-6)
})

test_that("Poisson information matrix is symmetric positive definite", {
	set.seed(3)
	n <- 300
	X <- cbind(1, matrix(rnorm(n * 4), n, 4))
	y <- rpois(n, exp(X %*% c(0.3, 0.1, -0.2, 0.4, -0.1)))
	fit <- EDI:::fast_poisson_regression_cpp(X, y)
	info <- fit$fisher_information
	expect_equal(info, t(info))
	expect_true(all(eigen(info, only.values = TRUE)$values > 0))
})

test_that("Poisson IRLS score is near zero at convergence", {
	set.seed(4)
	n <- 400
	X <- cbind(1, rnorm(n), rnorm(n)^2)
	y <- rpois(n, exp(X %*% c(1, -0.3, 0.1)))
	fit <- EDI:::fast_poisson_regression_cpp(X, y)
	expect_true(fit$converged)
	expect_true(max(abs(fit$score)) < 1e-4)
})

test_that("Poisson IRLS weighted_crossprod matches on intercept-only model", {
	set.seed(5)
	n <- 100
	X <- matrix(1, n, 1)
	y <- rpois(n, 3.7)
	fit <- EDI:::fast_poisson_regression_cpp(X, y)
	# MLE for intercept-only Poisson is log(mean(y))
	expect_equal(as.numeric(fit$b), log(mean(y)), tolerance = 1e-8)
})
