library(testthat)
library(EDI)

test_that("logistic GLMM with log1pexp_array_fast gives finite converged estimates", {
  set.seed(77)
  n_grp <- 30L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- as.integer(rep(c(0L, 0L, 1L, 1L), n_grp))
  x1 <- rnorm(n)
  re <- rnorm(n_grp, sd = 0.5)[group_id]
  y  <- as.numeric(rbinom(n, 1L, plogis(-0.5 + 0.8 * w + 0.3 * x1 + re)))
  X  <- cbind(1, w = w, x1 = x1)

  fit <- EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 1L, estimate_only = FALSE)
  expect_true(fit$converged)
  expect_true(all(is.finite(fit$b)))
  expect_true(is.finite(fit$log_sigma))
  expect_true(is.finite(fit$ssq_b_T))
})

test_that("logistic GLMM treatment coefficient direction is correct", {
  set.seed(88)
  n_grp <- 40L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- as.integer(rep(c(0L, 0L, 1L, 1L), n_grp))
  re <- rnorm(n_grp, sd = 0.3)[group_id]
  y  <- as.numeric(rbinom(n, 1L, plogis(2.0 * w + re)))
  X  <- cbind(1, w = w)

  fit <- EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 1L, estimate_only = TRUE)
  expect_gt(as.numeric(fit$b)[2L], 0)
})

test_that("logistic GLMM estimate_only and full run agree on coefficients", {
  set.seed(55)
  n_grp <- 25L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- as.integer(rep(c(0L, 0L, 1L, 1L), n_grp))
  x1 <- rnorm(n)
  re <- rnorm(n_grp, sd = 0.4)[group_id]
  y  <- as.numeric(rbinom(n, 1L, plogis(-0.3 + 0.6 * w + 0.2 * x1 + re)))
  X  <- cbind(1, w = w, x1 = x1)

  fit_est <- EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 1L, estimate_only = TRUE)
  fit_var <- EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 1L, estimate_only = FALSE)
  expect_equal(as.numeric(fit_est$b), as.numeric(fit_var$b), tolerance = 1e-8)
  expect_equal(fit_est$log_sigma, fit_var$log_sigma, tolerance = 1e-8)
})

test_that("logistic GLMM matches lme4::glmer on moderate dataset", {
  skip_if_not_installed("lme4")
  set.seed(33)
  n_grp <- 30L; grp_sz <- 4L; n <- n_grp * grp_sz
  group_id <- as.integer(rep(seq_len(n_grp), each = grp_sz))
  w  <- as.integer(rep(c(0L, 0L, 1L, 1L), n_grp))
  x1 <- rnorm(n)
  re <- rnorm(n_grp, sd = 0.4)[group_id]
  y  <- as.numeric(rbinom(n, 1L, plogis(-0.3 + 0.7 * w + 0.2 * x1 + re)))
  X  <- cbind(1, w = w, x1 = x1)

  fit_cpp <- EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 1L, estimate_only = FALSE)

  df <- data.frame(y = y, w = w, x1 = x1, grp = factor(group_id))
  fit_r <- lme4::glmer(y ~ w + x1 + (1 | grp), data = df, family = binomial)

  beta_cpp_w <- as.numeric(fit_cpp$b)[2L]
  beta_r_w   <- unname(lme4::fixef(fit_r)["w"])
  expect_equal(beta_cpp_w, beta_r_w, tolerance = 0.1,
    label = "logistic GLMM treatment beta vs lme4::glmer")
})
