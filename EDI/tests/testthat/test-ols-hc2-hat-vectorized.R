test_that("ols HC2 setup computes hat diagonal from vectorized row products", {
  set.seed(5501)
  n <- 160L
  x1 <- rnorm(n)
  x2 <- 0.35 * x1 + rnorm(n)
  x3 <- rbinom(n, 1L, 0.4)
  X <- cbind(1, x1, x2, x3, x1 * x3)

  setup <- EDI:::ols_hc2_setup_cpp(X)
  bread_ref <- solve(crossprod(X))
  hat_ref <- rowSums((X %*% bread_ref) * X)

  expect_equal(unname(setup$bread), unname(bread_ref), tolerance = 1e-11)
  expect_equal(unname(setup$hat), unname(hat_ref), tolerance = 1e-11)
})

test_that("ols HC2 setup hat diagonal feeds the same post-fit covariance", {
  set.seed(5502)
  n <- 120L
  X <- cbind(1, matrix(rnorm(n * 4L), n, 4L))
  beta <- c(0.3, -0.15, 0.2, 0.05, -0.1)
  y <- as.numeric(X %*% beta + rnorm(n, sd = 1 + 0.2 * abs(X[, 2L])))
  coef_hat <- as.numeric(lm.fit(X, y)$coefficients)

  setup <- EDI:::ols_hc2_setup_cpp(X)
  fit <- EDI:::ols_hc2_post_fit_precomputed_cpp(
    X,
    y,
    coef_hat,
    setup$bread,
    setup$hat,
    j_treat = 2L
  )

  resid <- as.numeric(y - X %*% coef_hat)
  omega <- resid^2 / pmax(1 - setup$hat, .Machine$double.eps)
  meat_ref <- crossprod(X * sqrt(omega))
  vcov_ref <- setup$bread %*% meat_ref %*% setup$bread
  vcov_ref <- 0.5 * (vcov_ref + t(vcov_ref))

  expect_equal(unname(fit$vcov), unname(vcov_ref), tolerance = 1e-11)
  expect_equal(fit$ssq_hat, vcov_ref[2L, 2L], tolerance = 1e-11)
})
