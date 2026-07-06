library(testthat)
library(EDI)

make_hnb_data <- function(seed = 1L, n = 300L, p = 3L,
                           theta_true = 10.0, mu_mean = 200) {
  set.seed(seed)
  beta <- c(log(mu_mean), 0.5, -0.3)[seq_len(p)]
  X    <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  mu   <- exp(X %*% beta)
  y    <- rnbinom(n, size = theta_true, mu = mu)
  y[y == 0L] <- 1L
  list(X = X, y = as.integer(y))
}

# ── correctness ───────────────────────────────────────────────────────────────

test_that("fast-path data: >99% obs have p0 < 1e-7 (fast path active)", {
  d    <- make_hnb_data(seed = 1L, theta_true = 10.0, mu_mean = 200)
  r    <- (10.0 / (10.0 + exp(d$X %*% c(log(200), 0.5, -0.3))))^10
  expect_true(mean(r < 1e-7) > 0.99,
    label = "at least 99% of obs should hit the fast path")
})

test_that("converges and all estimates finite (fast-path data)", {
  d <- make_hnb_data(1L)
  r <- EDI:::fast_hurdle_negbin_cpp(d$X, d$y, d$X,
         smart_cold_start = TRUE, estimate_only = TRUE)
  expect_true(r$converged)
  expect_true(all(is.finite(r$b)))
  expect_true(is.finite(r$neg_ll))
})

test_that("estimate_only=TRUE and FALSE give same point estimates", {
  d  <- make_hnb_data(2L)
  r1 <- EDI:::fast_hurdle_negbin_cpp(d$X, d$y, d$X,
          smart_cold_start = TRUE, estimate_only = TRUE)
  r2 <- EDI:::fast_hurdle_negbin_with_var_cpp(d$X, d$y, d$X, j = 2L,
          smart_cold_start = TRUE)
  expect_equal(r1$b, r2$b, tolerance = 1e-5)
})

test_that("observed information is finite, symmetric, and PD", {
  d  <- make_hnb_data(3L)
  r  <- EDI:::fast_hurdle_negbin_with_var_cpp(d$X, d$y, d$X, j = 2L,
          smart_cold_start = TRUE)
  I  <- r$observed_information
  expect_true(is.matrix(I))
  expect_true(all(is.finite(I)))
  expect_equal(I, t(I), tolerance = 1e-12)
  expect_true(all(eigen(I, only.values = TRUE)$values > 0))
})

test_that("score near zero at convergence (log1p fast path doesn't corrupt gradient)", {
  d   <- make_hnb_data(4L)
  r   <- EDI:::fast_hurdle_negbin_cpp(d$X, d$y, d$X,
           smart_cold_start = TRUE, estimate_only = TRUE)
  par <- c(r$b, log(r$theta_hat))
  sc  <- as.numeric(EDI:::get_hurdle_negbin_count_score_cpp(d$X, d$y, par))
  expect_true(max(abs(sc[seq_len(ncol(d$X))])) < 0.1,
    label = "beta score near zero at convergence")
})

test_that("repeated calls return identical results (no state mutation)", {
  d  <- make_hnb_data(5L)
  r1 <- EDI:::fast_hurdle_negbin_cpp(d$X, d$y, d$X,
          smart_cold_start = TRUE, estimate_only = TRUE)
  r2 <- EDI:::fast_hurdle_negbin_cpp(d$X, d$y, d$X,
          smart_cold_start = TRUE, estimate_only = TRUE)
  expect_identical(r1$b, r2$b)
  expect_identical(r1$neg_ll, r2$neg_ll)
})

test_that("mixed fast/slow path data (theta=1, small mu): also converges", {
  # theta=1, mu=5: p0=(1/6)^1 ≈ 0.17 — all obs use std::log branch
  set.seed(6L)
  n <- 200L; p <- 3L
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y <- pmax(1L, rnbinom(n, size = 1, mu = 5))
  r <- EDI:::fast_hurdle_negbin_cpp(X, y, X,
         smart_cold_start = TRUE, estimate_only = TRUE)
  expect_true(r$converged)
  expect_true(all(is.finite(r$b)))
})

test_that("warm start matches cold start on fast-path data", {
  d    <- make_hnb_data(7L)
  cold <- EDI:::fast_hurdle_negbin_cpp(d$X, d$y, d$X,
            smart_cold_start = TRUE, estimate_only = TRUE)
  warm_par <- c(cold$b, log(cold$theta_hat))
  warm <- EDI:::fast_hurdle_negbin_cpp(d$X, d$y, d$X,
            warm_start_params = warm_par,
            smart_cold_start = FALSE, estimate_only = TRUE)
  expect_equal(cold$b, warm$b, tolerance = 1e-4)
  expect_equal(cold$neg_ll, warm$neg_ll, tolerance = 1e-3)
})
