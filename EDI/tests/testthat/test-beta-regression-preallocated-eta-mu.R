make_beta_prealloc_data <- function(n = 180L, p = 4L, seed = 5601L, phi = 11) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  beta <- seq(-0.2, 0.25, length.out = p)
  beta[1L] <- -0.1
  mu <- plogis(drop(X %*% beta))
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  list(X = X, y = y)
}

test_that("beta regression preallocated eta/mu path is repeatable and stationary", {
  dat <- make_beta_prealloc_data()

  fit1 <- EDI:::fast_beta_regression_cpp(dat$X, dat$y, estimate_only = TRUE)
  fit2 <- EDI:::fast_beta_regression_cpp(dat$X, dat$y, estimate_only = TRUE)

  expect_true(isTRUE(fit1$converged))
  expect_true(isTRUE(fit2$converged))
  expect_equal(fit2$coefficients, fit1$coefficients, tolerance = 1e-12)
  expect_equal(fit2$phi, fit1$phi, tolerance = 1e-12)
  expect_equal(fit2$neg_loglik, fit1$neg_loglik, tolerance = 1e-10)

  params <- c(fit1$coefficients, log(fit1$phi))
  score <- EDI:::get_beta_regression_score_cpp(dat$X, dat$y, params)
  expect_lt(max(abs(score)) / nrow(dat$X), 5e-4)
})

test_that("beta regression preallocated eta/mu path handles weights and fixed parameters", {
  dat <- make_beta_prealloc_data(n = 150L, p = 5L, seed = 5602L, phi = 9)
  set.seed(5603)
  weights <- runif(nrow(dat$X), 0.25, 1.75)

  fit_w1 <- EDI:::fast_beta_regression_weighted_cpp(
    dat$X,
    dat$y,
    weights,
    estimate_only = TRUE
  )
  fit_w2 <- EDI:::fast_beta_regression_weighted_cpp(
    dat$X,
    dat$y,
    weights,
    estimate_only = TRUE
  )

  expect_true(isTRUE(fit_w1$converged))
  expect_true(isTRUE(fit_w2$converged))
  expect_equal(fit_w2$coefficients, fit_w1$coefficients, tolerance = 1e-12)
  expect_equal(fit_w2$phi, fit_w1$phi, tolerance = 1e-12)

  fixed_idx <- 2L
  fixed_value <- 0.05
  fit_fixed <- EDI:::fast_beta_regression_cpp(
    dat$X,
    dat$y,
    fixed_idx = fixed_idx,
    fixed_values = fixed_value,
    estimate_only = TRUE
  )

  expect_true(isTRUE(fit_fixed$converged))
  expect_equal(fit_fixed$coefficients[fixed_idx], fixed_value, tolerance = 0)
})
