library(testthat)
library(EDI)

test_that("fast_ordinal_glmm_cpp estimate_only matches estimate from full run", {
  set.seed(42)
  n_grp <- 20L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- rep(c(0L, 0L, 1L, 1L), n_grp)
  x1 <- rnorm(n)
  re <- rnorm(n_grp, sd = 0.4)[group_id]
  eta <- -0.3 + 0.5 * w + 0.2 * x1 + re
  p1 <- plogis(1 - eta); p2 <- plogis(-1 - eta)
  prob <- cbind(p1 - p2, p1 - (p1 - p2), 1 - p1)
  prob <- pmax(prob, 1e-6); prob <- prob / rowSums(prob)
  y <- as.integer(apply(prob, 1L, function(p) sample(1:3, 1L, prob = p)))
  X <- cbind(w = w, x1 = x1)

  fit_est <- fast_ordinal_glmm_cpp(X, y, group_id, K = 3L, j_T = 0L, estimate_only = TRUE)
  fit_var <- fast_ordinal_glmm_cpp(X, y, group_id, K = 3L, j_T = 0L, estimate_only = FALSE)

  # Coefficients from estimate_only must equal those from the full run
  expect_equal(as.numeric(fit_est$b), as.numeric(fit_var$b), tolerance = 1e-8)
  expect_equal(fit_est$log_sigma, fit_var$log_sigma, tolerance = 1e-8)
})

test_that("fast_ordinal_glmm_cpp output is finite and converged", {
  set.seed(7)
  n_grp <- 30L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- rep(c(0L, 1L), n)
  x1 <- rnorm(n)
  re <- rnorm(n_grp, sd = 0.5)[group_id]
  eta <- 0.4 * w + 0.3 * x1 + re
  p1 <- plogis(1 - eta); p2 <- plogis(-1 - eta)
  prob <- cbind(p1-p2, p1-(p1-p2), 1-p1)
  prob <- pmax(prob, 1e-6); prob <- prob / rowSums(prob)
  y <- as.integer(apply(prob, 1L, function(p) sample(1:3, 1L, prob = p)))
  X <- cbind(w = w, x1 = x1)

  fit <- fast_ordinal_glmm_cpp(X, y, group_id, K = 3L, j_T = 0L, estimate_only = FALSE)
  expect_true(all(is.finite(as.numeric(fit$b))))
  expect_true(is.finite(fit$log_sigma))
  expect_true(is.finite(fit$ssq_b_T))
  expect_true(fit$converged)
})

test_that("fast_ordinal_glmm_cpp treatment effect direction is correct", {
  set.seed(99)
  n_grp <- 40L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- as.integer(rep(c(0L, 0L, 1L, 1L), n_grp))
  re <- rnorm(n_grp, sd = 0.3)[group_id]
  # Strong positive treatment effect: higher w -> higher ordinal category
  eta <- 1.5 * w + re
  p1 <- plogis(2 - eta); p2 <- plogis(-2 - eta)
  prob <- cbind(p1-p2, p1-(p1-p2), 1-p1)
  prob <- pmax(prob, 1e-6); prob <- prob / rowSums(prob)
  y <- as.integer(apply(prob, 1L, function(p) sample(1:3, 1L, prob = p)))
  X <- matrix(w, ncol = 1L, dimnames = list(NULL, "w"))

  fit <- fast_ordinal_glmm_cpp(X, y, group_id, K = 3L, j_T = 0L, estimate_only = TRUE)
  # With a strong positive treatment effect the treatment coefficient should be positive
  expect_gt(as.numeric(fit$b)[1L], 0)
})

test_that("fast_ordinal_glmm_cpp matches ordinal::clmm on moderate dataset", {
  skip_if_not_installed("ordinal")
  set.seed(55)
  n_grp <- 25L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- as.integer(rep(c(0L, 0L, 1L, 1L), n_grp))
  re <- rnorm(n_grp, sd = 0.4)[group_id]
  eta <- 0.6 * w + re
  p1 <- plogis(1 - eta); p2 <- plogis(-1 - eta)
  prob <- cbind(p1-p2, p1-(p1-p2), 1-p1)
  prob <- pmax(prob, 1e-6); prob <- prob / rowSums(prob)
  y <- as.integer(apply(prob, 1L, function(p) sample(1:3, 1L, prob = p)))
  X <- matrix(w, ncol = 1L, dimnames = list(NULL, "w"))

  fit_cpp <- fast_ordinal_glmm_cpp(X, y, group_id, K = 3L, j_T = 0L, estimate_only = FALSE)

  df <- data.frame(y = ordered(y), w = w, grp = factor(group_id))
  fit_r <- ordinal::clmm(y ~ w + (1 | grp), data = df)

  beta_cpp <- as.numeric(fit_cpp$b)[1L]  # treatment coefficient
  beta_r   <- unname(coef(fit_r)["w"])
  expect_equal(beta_cpp, beta_r, tolerance = 0.1,
    label = "ordinal GLMM treatment beta vs ordinal::clmm")
})
