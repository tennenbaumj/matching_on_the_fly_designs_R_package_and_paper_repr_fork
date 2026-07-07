make_logbin_var_data <- function(n = 220L, p = 5L, seed = 5701L) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L), sd = 0.45), n, p - 1L))
  beta <- c(-2.6, seq(-0.12, 0.12, length.out = p - 1L))
  eta <- pmin(drop(X %*% beta), -0.05)
  y <- rbinom(n, 1L, exp(eta))
  list(X = X, y = y)
}

make_identity_var_data <- function(n = 180L, p = 4L, seed = 5702L) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L), sd = 0.25), n, p - 1L))
  beta <- c(0.28, seq(-0.035, 0.035, length.out = p - 1L))
  pr <- pmin(pmax(drop(X %*% beta), 0.05), 0.85)
  y <- rbinom(n, 1L, pr)
  list(X = X, y = y)
}

expect_logbin_var_matches_fisher <- function(fit, X, link = c("log", "identity"), j,
                                             fixed_idx = integer()) {
  link <- match.arg(link)
  beta <- as.numeric(fit$b)
  eta <- drop(X %*% beta)
  if (link == "log") {
    eta <- pmin(eta, -1e-8)
    mu <- exp(eta)
    w <- pmax(mu / pmax(1 - mu, 1e-8), 1e-8)
  } else {
    mu <- pmin(pmax(eta, 1e-8), 1 - 1e-8)
    w <- pmax(1 / pmax(mu * (1 - mu), 1e-8), 1e-8)
  }
  fisher_ref <- crossprod(X * sqrt(w))
  free_idx <- setdiff(seq_len(ncol(X)), fixed_idx)
  free_j <- match(j, free_idx)

  expect_equal(unname(fit$fisher_information), unname(fisher_ref), tolerance = 1e-10)
  if (is.na(free_j)) {
    expect_true(is.na(fit$ssq_b_j))
  } else {
    ssq_ref <- solve(fisher_ref[free_idx, free_idx, drop = FALSE])[free_j, free_j]
    expect_equal(fit$ssq_b_j, ssq_ref, tolerance = 1e-10)
  }
}

test_that("log-binomial variance path reuses returned Fisher information", {
  dat <- make_logbin_var_data()

  fit <- EDI:::fast_log_binomial_regression_with_var_cpp(dat$X, dat$y, j = 3L)

  expect_true(isTRUE(fit$converged))
  expect_logbin_var_matches_fisher(fit, dat$X, link = "log", j = 3L)
})

test_that("log-binomial variance Fisher reuse respects fixed parameters", {
  dat <- make_logbin_var_data(n = 200L, p = 5L, seed = 5703L)
  fixed_idx <- 2L
  fixed_value <- -0.08

  fit <- EDI:::fast_log_binomial_regression_with_var_cpp(
    dat$X,
    dat$y,
    j = 4L,
    fixed_idx = fixed_idx,
    fixed_values = fixed_value
  )

  expect_true(isTRUE(fit$converged))
  expect_equal(fit$b[fixed_idx], fixed_value, tolerance = 0)
  expect_logbin_var_matches_fisher(
    fit,
    dat$X,
    link = "log",
    j = 4L,
    fixed_idx = fixed_idx
  )
})

test_that("identity-binomial variance path uses the same Fisher reuse helper", {
  dat <- make_identity_var_data()

  fit <- EDI:::fast_identity_binomial_regression_with_var_cpp(dat$X, dat$y, j = 2L)

  expect_true(isTRUE(fit$converged))
  expect_logbin_var_matches_fisher(fit, dat$X, link = "identity", j = 2L)
})
