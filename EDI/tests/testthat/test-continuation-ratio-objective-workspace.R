make_continuation_ratio_workspace_data <- function() {
  set.seed(6101)
  n <- 180L
  p <- 4L
  X <- matrix(rnorm(n * p, sd = 0.65), n, p)
  eta <- drop(0.3 * X[, 1] - 0.2 * X[, 2] + 0.1 * X[, 3])
  logits <- cbind(
    -0.5 + 0.4 * eta,
    -0.1 + 0.7 * eta,
    0.3 + 0.9 * eta
  )
  hazards <- plogis(logits)
  y <- integer(n)
  for (i in seq_len(n)) {
    y[i] <- 4L
    for (k in seq_len(3L)) {
      if (runif(1) < hazards[i, k]) {
        y[i] <- k
        break
      }
    }
  }
  list(X = X, y = as.numeric(y))
}

test_that("continuation-ratio objective workspace is repeatable and resets across optimizer calls", {
  d <- make_continuation_ratio_workspace_data()

  fit1 <- EDI:::fast_continuation_ratio_regression_cpp(d$X, d$y, maxit = 200L, tol = 1e-8)
  fit2 <- EDI:::fast_continuation_ratio_regression_cpp(d$X, d$y, maxit = 200L, tol = 1e-8)

  expect_true(fit1$converged)
  expect_true(fit2$converged)
  expect_equal(fit2$params, fit1$params, tolerance = 0)

  score <- EDI:::get_continuation_ratio_regression_score_cpp(d$X, d$y, fit1$params)
  observed_information <- -EDI:::get_continuation_ratio_regression_hessian_cpp(d$X, d$y, fit1$params)

  expect_lt(sqrt(sum(score^2)), 1e-2)
  expect_equal(unname(fit1$fisher_information), unname(observed_information), tolerance = 1e-10)
})

test_that("continuation-ratio variance path reuses workspace without corrupting information", {
  d <- make_continuation_ratio_workspace_data()
  fit <- EDI:::fast_continuation_ratio_regression_with_var_cpp(
    d$X,
    d$y,
    maxit = 200L,
    tol = 1e-8
  )

  observed_information <- -EDI:::get_continuation_ratio_regression_hessian_cpp(
    d$X,
    d$y,
    fit$params
  )

  expect_true(fit$converged)
  expect_equal(unname(fit$fisher_information), unname(observed_information), tolerance = 1e-10)
  expect_equal(unname(fit$vcov), unname(solve(observed_information)), tolerance = 1e-8)
})
