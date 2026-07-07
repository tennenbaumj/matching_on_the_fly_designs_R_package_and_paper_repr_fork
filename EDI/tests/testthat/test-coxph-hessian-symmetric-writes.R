cox_loglik_reference <- function(X, y, dead, beta) {
  eta <- drop(X %*% beta)
  event_times <- sort(unique(y[dead > 0.5]))
  sum(vapply(event_times, function(tk) {
    events <- dead > 0.5 & y == tk
    risk <- y >= tk
    sum(eta[events]) - sum(events) * log(sum(exp(eta[risk])))
  }, numeric(1L)))
}

stratified_cox_loglik_reference <- function(X, y, dead, strata, beta) {
  sum(vapply(split(seq_along(y), strata), function(idx) {
    cox_loglik_reference(X[idx, , drop = FALSE], y[idx], dead[idx], beta)
  }, numeric(1L)))
}

test_that("coxph Hessian symmetric writes match independent numerical derivatives", {
  skip_if_not_installed("numDeriv")
  set.seed(5301)
  n <- 80L
  p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(0.35, -0.25, 0.15)
  y <- rexp(n, rate = 0.12 * exp(drop(X %*% beta)))
  dead <- as.numeric(seq_len(n) %% 5L != 0L)
  params <- c(0.2, -0.1, 0.05)

  hessian <- EDI:::get_coxph_hessian_cpp(X, y, dead, params)
  hessian_ref <- numDeriv::hessian(
    function(par) cox_loglik_reference(X, y, dead, par),
    params
  )

  expect_equal(unname(hessian), unname(hessian_ref), tolerance = 1e-5)
  expect_equal(unname(hessian), unname(t(hessian)), tolerance = 0)
})

test_that("stratified coxph Hessian symmetric writes match independent numerical derivatives", {
  skip_if_not_installed("numDeriv")
  set.seed(5302)
  n <- 96L
  p <- 2L
  strata <- rep(seq_len(4L), each = n / 4L)
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(0.3, -0.2)
  y <- rexp(n, rate = rep(c(0.08, 0.11, 0.14, 0.17), each = n / 4L) * exp(drop(X %*% beta)))
  dead <- as.numeric(seq_len(n) %% 4L != 0L)
  params <- c(0.15, -0.05)

  hessian <- EDI:::get_stratified_coxph_hessian_cpp(X, y, dead, as.integer(strata), params)
  hessian_ref <- numDeriv::hessian(
    function(par) stratified_cox_loglik_reference(X, y, dead, strata, par),
    params
  )

  expect_equal(unname(hessian), unname(hessian_ref), tolerance = 1e-5)
  expect_equal(unname(hessian), unname(t(hessian)), tolerance = 0)
})

test_that("coxph Hessian symmetric writes are repeatable", {
  set.seed(5303)
  n <- 70L
  p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rexp(n, rate = 0.1)
  dead <- as.numeric(seq_len(n) %% 6L != 0L)
  params <- c(0.2, -0.1, 0.05)

  hessian1 <- EDI:::get_coxph_hessian_cpp(X, y, dead, params)
  hessian2 <- EDI:::get_coxph_hessian_cpp(X, y, dead, params)

  expect_equal(hessian1, hessian2, tolerance = 0)
})
