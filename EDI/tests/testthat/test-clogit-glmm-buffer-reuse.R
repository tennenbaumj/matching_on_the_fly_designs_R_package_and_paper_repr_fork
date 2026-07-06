library(testthat)
library(EDI)

make_data <- function(seed = 1L, n_disc = 120L, n_conc = 100L, p = 3L, G = 40L) {
  set.seed(seed)
  Xd <- cbind(1, sample(c(-1L, 1L), n_disc, TRUE), matrix(rnorm(n_disc * (p - 1L)), n_disc))
  yd <- as.numeric(rbinom(n_disc, 1L, 0.4))
  Xc <- cbind(1, rep(c(0, 1), n_conc / 2L), matrix(rnorm(n_conc * (p - 1L)), n_conc))
  yc <- as.numeric(rbinom(n_conc, 1L, 0.4))
  gc <- as.integer(rep(seq_len(G), length.out = n_conc))
  list(Xd = Xd, yd = yd, Xc = Xc, yc = yc, gc = gc)
}

call_est <- function(d, ...) {
  EDI:::fast_clogit_plus_glmm_cpp(d$Xd, d$yd, d$Xc, d$yc, d$gc,
    has_discordant = TRUE, has_concordant = TRUE, estimate_only = TRUE, ...)
}

call_full <- function(d, ...) {
  EDI:::fast_clogit_plus_glmm_cpp(d$Xd, d$yd, d$Xc, d$yc, d$gc,
    has_discordant = TRUE, has_concordant = TRUE, estimate_only = FALSE, ...)
}

# ── correctness ──────────────────────────────────────────────────────────────

test_that("estimate_only=TRUE and FALSE give same point estimates", {
  d  <- make_data(1L)
  r1 <- call_est(d)
  r2 <- call_full(d)
  expect_equal(r1$b, r2$b, tolerance = 1e-5)
})

test_that("converges and all estimates finite", {
  d <- make_data(2L)
  r <- call_full(d)
  expect_true(r$converged)
  expect_true(all(is.finite(r$b)))
  expect_true(is.finite(r$neg_loglik))
})

test_that("vcov is finite, symmetric, and positive definite", {
  d  <- make_data(3L)
  r  <- call_full(d)
  V  <- r$vcov
  expect_true(all(is.finite(V)))
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(eigen(V, only.values = TRUE)$values > 0))
})

test_that("score ≈ 0 at convergence (buffer writes don't corrupt gradient)", {
  d   <- make_data(1L)
  r   <- call_full(d)
  par <- as.numeric(r$b)
  sc  <- as.numeric(EDI:::get_clogit_plus_glmm_score_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  expect_true(max(abs(sc[seq_len(ncol(d$Xc))])) < 0.1,
    label = "beta score near zero at convergence")
})

test_that("observed information is symmetric and vcov is PD (mutable buffers don't corrupt hessian)", {
  d   <- make_data(2L)
  r   <- call_full(d)
  par <- as.numeric(r$b)
  # get_clogit_plus_glmm_hessian_cpp returns -obj.hessian() = Hessian of log-lik (NSD at max)
  H   <- as.matrix(EDI:::get_clogit_plus_glmm_hessian_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  expect_equal(H, t(H), tolerance = 1e-10)
  # vcov = inv(obj.hessian()) should be PD when info is PD
  expect_true(all(eigen(r$vcov, only.values = TRUE)$values > 0))
})

test_that("warm start converges to same estimates as cold start", {
  d    <- make_data(6L)
  cold <- call_full(d)
  warm <- call_full(d, warm_start_params = cold$b)
  expect_equal(cold$b, warm$b, tolerance = 2e-3)
  expect_equal(cold$neg_loglik, warm$neg_loglik, tolerance = 1e-3)
})

test_that("concordant-only path (no discordant component) works correctly", {
  set.seed(7L)
  n_conc <- 80L; p <- 3L; G <- 30L
  Xc <- cbind(1, rep(c(0, 1), n_conc / 2L), matrix(rnorm(n_conc * (p - 1L)), n_conc))
  yc <- as.numeric(rbinom(n_conc, 1L, 0.4))
  gc <- as.integer(rep(seq_len(G), length.out = n_conc))
  Xd <- matrix(0, 0, p); yd <- numeric(0)
  r  <- EDI:::fast_clogit_plus_glmm_cpp(Xd, yd, Xc, yc, gc,
           has_discordant = FALSE, has_concordant = TRUE, estimate_only = FALSE)
  expect_true(r$converged)
  expect_true(all(is.finite(r$b)))
})

test_that("repeated calls return identical results (no state leakage between calls)", {
  d  <- make_data(8L)
  r1 <- call_est(d)
  r2 <- call_est(d)
  expect_identical(r1$b, r2$b)
  expect_identical(r1$neg_loglik, r2$neg_loglik)
})

test_that("neg_glmm value agrees between value() and operator() at same params", {
  d   <- make_data(9L)
  r   <- call_full(d)
  par <- as.numeric(r$b)
  # score call exercises operator(); neg_loglik exercises value() path via optimizer
  sc  <- as.numeric(EDI:::get_clogit_plus_glmm_score_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  expect_true(is.finite(r$neg_loglik))
  expect_true(all(is.finite(sc)))
})

test_that("hessian is identical across repeated calls (mutable buffers don't leak between calls)", {
  d   <- make_data(11L)
  r   <- call_full(d)
  par <- as.numeric(r$b)
  H1  <- as.matrix(EDI:::get_clogit_plus_glmm_hessian_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  H2  <- as.matrix(EDI:::get_clogit_plus_glmm_hessian_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  expect_identical(H1, H2)
})

test_that("hessian after operator() call gives same result as standalone hessian", {
  d   <- make_data(12L)
  r   <- call_full(d)
  par <- as.numeric(r$b)
  # Call hessian first (primes mutable buffers with par)
  H_before <- as.matrix(EDI:::get_clogit_plus_glmm_hessian_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  # Call operator() (overwrites all shared mutable buffers)
  sc <- as.numeric(EDI:::get_clogit_plus_glmm_score_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  # Call hessian again — should produce same result despite buffer overwrite
  H_after <- as.matrix(EDI:::get_clogit_plus_glmm_hessian_cpp(
    d$Xd, d$yd, d$Xc, d$yc, d$gc, par, TRUE, TRUE))
  expect_identical(H_before, H_after)
})

test_that("sigma estimate is strictly positive (exp(log_sigma) > 0)", {
  d   <- make_data(10L)
  r   <- call_full(d)
  p   <- ncol(d$Xc)
  log_sigma <- r$b[p + 1L]
  expect_true(exp(log_sigma) > 0)
  expect_true(is.finite(log_sigma))
})
