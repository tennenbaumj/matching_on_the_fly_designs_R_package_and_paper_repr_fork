test_that("fractional-logit post-fit statistics match an independent R reference", {
  set.seed(8701)
  n <- 320L
  X <- cbind(
    intercept = 1,
    treatment = rep(0:1, length.out = n),
    x1 = rnorm(n),
    x2 = runif(n, -1, 1),
    x3 = rnorm(n)
  )
  coef_hat <- c(-0.3, 0.55, -0.2, 0.35, 0.1)
  mu_hat <- plogis(drop(X %*% coef_hat))
  y <- rbeta(n, 12 * mu_hat, 12 * (1 - mu_hat))

  actual <- EDI:::gcomp_fractional_logit_post_fit_cpp(
    X, y, coef_hat, mu_hat, 2L
  )

  working_weights <- mu_hat * (1 - mu_hat)
  bread <- solve(crossprod(X, X * working_weights))
  score_residual <- y - mu_hat
  meat <- crossprod(X, X * score_residual^2)
  vcov <- bread %*% meat %*% bread
  vcov <- 0.5 * (vcov + t(vcov))

  eta <- drop(X %*% coef_hat)
  eta_base <- eta - coef_hat[[2L]] * X[, 2L]
  mean1_i <- plogis(eta_base + coef_hat[[2L]])
  mean0_i <- plogis(eta_base)
  mean1 <- mean(mean1_i)
  mean0 <- mean(mean0_i)
  wt1 <- mean1_i * (1 - mean1_i) / n
  wt0 <- mean0_i * (1 - mean0_i) / n
  grad1 <- drop(crossprod(X, wt1))
  grad0 <- drop(crossprod(X, wt0))
  grad1[[2L]] <- sum(wt1)
  grad0[[2L]] <- 0
  grad_md <- grad1 - grad0

  expect_named(
    actual,
    c("vcov", "std_err", "z_vals", "mean1", "mean0", "md", "se_md"),
    ignore.order = FALSE
  )
  expect_equal(unname(actual$vcov), unname(vcov), tolerance = 1e-12)
  expect_equal(unname(actual$std_err), unname(sqrt(diag(vcov))), tolerance = 1e-12)
  expect_equal(unname(actual$z_vals), unname(coef_hat / sqrt(diag(vcov))), tolerance = 1e-12)
  expect_equal(actual$mean1, mean1, tolerance = 1e-14)
  expect_equal(actual$mean0, mean0, tolerance = 1e-14)
  expect_equal(actual$md, mean1 - mean0, tolerance = 1e-14)
  expect_equal(
    actual$se_md,
    sqrt(drop(t(grad_md) %*% vcov %*% grad_md)),
    tolerance = 1e-12
  )

  logistic <- EDI:::gcomp_logistic_post_fit_cpp(
    X, y, coef_hat, mu_hat, 2L
  )
  expect_equal(logistic$vcov, actual$vcov, tolerance = 0)
  expect_equal(logistic$risk1, actual$mean1, tolerance = 0)
  expect_equal(logistic$risk0, actual$mean0, tolerance = 0)
  expect_equal(logistic$rd, actual$md, tolerance = 0)
  expect_equal(logistic$se_rd, actual$se_md, tolerance = 0)

  grad_log_rr <- grad1 / mean1 - grad0 / mean0
  expect_equal(logistic$log_rr, log(mean1 / mean0), tolerance = 1e-14)
  expect_equal(logistic$rr, mean1 / mean0, tolerance = 1e-14)
  expect_equal(
    logistic$se_log_rr,
    sqrt(drop(t(grad_log_rr) %*% vcov %*% grad_log_rr)),
    tolerance = 1e-12
  )
})

test_that("fractional-logit post-fit validation is retained", {
  X <- cbind(1, c(0, 1, 0, 1), c(-1, 0, 1, 2))
  coef_hat <- c(0.1, 0.3, -0.2)
  mu_hat <- plogis(drop(X %*% coef_hat))
  y <- c(0.2, 0.7, 0.4, 0.8)

  expect_error(
    EDI:::gcomp_fractional_logit_post_fit_cpp(X, y, coef_hat, mu_hat, 0L),
    "treatment column index is out of bounds"
  )
  expect_error(
    EDI:::gcomp_fractional_logit_post_fit_cpp(X, y[-1L], coef_hat, mu_hat, 2L),
    "dimension mismatch"
  )
  expect_error(
    EDI:::gcomp_fractional_logit_post_fit_cpp(X, y, coef_hat, replace(mu_hat, 1L, 1), 2L),
    "boundary fitted values"
  )
})
