make_zinb_workspace_data <- function(seed = 11001L, n = 180L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  z1 <- rnorm(n)
  X <- cbind(1, x1, x2)
  Xzi <- cbind(1, z1, x1)
  beta <- c(0.35, 0.25, -0.15)
  gamma <- c(-0.75, 0.30, -0.20)
  theta <- 1.7
  mu <- exp(drop(X %*% beta))
  p_zero <- plogis(drop(Xzi %*% gamma))
  y <- ifelse(runif(n) < p_zero, 0, rnbinom(n, size = theta, mu = mu))
  list(X = X, Xzi = Xzi, y = as.numeric(y), start = c(beta, gamma, log(theta)))
}

test_that("ZINB operator workspace is repeatable across estimate-only fits", {
  d <- make_zinb_workspace_data()

  fit1 <- EDI:::fast_zinb_cpp(
    d$X,
    d$Xzi,
    d$y,
    warm_start_params = d$start,
    maxit = 300L,
    tol = 1e-8,
    smart_cold_start = FALSE,
    estimate_only = TRUE
  )
  fit2 <- EDI:::fast_zinb_cpp(
    d$X,
    d$Xzi,
    d$y,
    warm_start_params = d$start,
    maxit = 300L,
    tol = 1e-8,
    smart_cold_start = FALSE,
    estimate_only = TRUE
  )

  expect_true(fit1$converged)
  expect_true(fit2$converged)
  expect_equal(fit2$params, fit1$params, tolerance = 0)
})

test_that("ZINB variance path remains finite after operator workspace reuse", {
  d <- make_zinb_workspace_data()
  fit0 <- EDI:::fast_zinb_cpp(
    d$X,
    d$Xzi,
    d$y,
    warm_start_params = d$start,
    maxit = 300L,
    tol = 1e-8,
    smart_cold_start = FALSE,
    estimate_only = TRUE
  )
  fit <- EDI:::fast_zinb_cpp(
    d$X,
    d$Xzi,
    d$y,
    warm_start_params = fit0$params,
    maxit = 150L,
    tol = 1e-8,
    smart_cold_start = FALSE,
    estimate_only = FALSE
  )

  expect_true(fit$converged)
  expect_true(all(is.finite(fit$score)))
  expect_equal(unname(fit$hessian), unname(t(fit$hessian)), tolerance = 1e-7)
  expect_true(all(is.finite(diag(fit$vcov))))
})
