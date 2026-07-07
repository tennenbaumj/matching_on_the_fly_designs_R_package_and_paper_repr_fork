test_that("OLS symmetric_crossprod matches lm() on standard data", {
	set.seed(1)
	n <- 300
	X <- cbind(1, matrix(rnorm(n * 4), n, 4))
	beta_true <- c(2, 0.5, -0.3, 0.7, -0.1)
	y <- X %*% beta_true + rnorm(n)
	fit <- EDI:::fast_ols_with_var_cpp(X, y)
	ref <- lm.fit(X, y)
	expect_equal(as.numeric(fit$b), unname(ref$coefficients), tolerance = 1e-10)
})

test_that("OLS XtX matrix is symmetric and positive definite", {
	set.seed(2)
	n <- 200
	X <- cbind(1, rnorm(n), rnorm(n), rnorm(n))
	y <- X %*% c(1, -0.5, 0.3, 0.2) + rnorm(n)
	fit <- EDI:::fast_ols_with_var_cpp(X, y)
	XtX <- fit$XtX
	expect_equal(XtX, t(XtX))
	expect_true(all(eigen(XtX, only.values = TRUE)$values > 0))
})

test_that("OLS with fixed params uses symmetric_crossprod correctly", {
	set.seed(3)
	n <- 250
	X <- cbind(1, matrix(rnorm(n * 3), n, 3))
	beta_true <- c(1, 0.5, -0.3, 0.2)
	y <- X %*% beta_true + rnorm(n)
	# Fix first param (intercept) to 1
	fit <- EDI:::fast_ols_with_var_cpp(X, y, fixed_idx = 1L, fixed_values = 1)
	y_adj <- y - X[, 1]
	ref <- lm.fit(X[, -1], y_adj)
	expect_equal(as.numeric(fit$b[-1]), unname(ref$coefficients), tolerance = 1e-10)
})

test_that("OLS sigma2_hat is positive and matches manual calculation", {
	set.seed(4)
	n <- 400
	p <- 5
	X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
	beta_true <- rnorm(p, 0, 0.5)
	sigma_true <- 1.5
	y <- X %*% beta_true + rnorm(n, 0, sigma_true)
	fit <- EDI:::fast_ols_with_var_cpp(X, y)
	expect_true(fit$sigma2_hat > 0)
	# Should be within 2 SEs of truth
	expect_true(abs(sqrt(fit$sigma2_hat) - sigma_true) < 2 * sigma_true / sqrt(2 * (n - p)))
})

test_that("OLS intercept-only model matches sample mean", {
	set.seed(5)
	n <- 100
	y <- rnorm(n, 3.5)
	X <- matrix(1, n, 1)
	fit <- EDI:::fast_ols_with_var_cpp(X, y)
	expect_equal(as.numeric(fit$b), mean(y), tolerance = 1e-12)
})

test_that("Robust regression XtX path gives finite result", {
	set.seed(6)
	n <- 150
	X <- cbind(1, rnorm(n), rnorm(n))
	y <- X %*% c(1, 0.5, -0.3) + rnorm(n)
	# Add outlier to trigger robust path
	y[1] <- y[1] + 20
	fit <- EDI:::fast_robust_regression_cpp(X, y)
	expect_true(is.finite(fit$coefficients[1]))
	expect_true(is.finite(fit$ssq_b_j))
})
