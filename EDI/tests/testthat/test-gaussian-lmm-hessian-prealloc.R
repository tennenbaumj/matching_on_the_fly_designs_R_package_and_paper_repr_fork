library(testthat)
library(EDI)

make_lmm_data <- function(seed = 1L, G = 60L) {
  set.seed(seed)
  n        <- G * 2L  # matched pairs
  group_id <- as.integer(rep(seq_len(G), each = 2L))
  X        <- cbind(1, rbinom(n, 1, 0.5), rnorm(n))
  b_true   <- c(0.8, 0.4, -0.2)
  u_true   <- rnorm(G, 0, 0.35)
  y        <- as.numeric(X %*% b_true + u_true[group_id] + rnorm(n, 0, 0.7))
  list(X = X, y = y, group_id = group_id)
}

test_that("gaussian_lmm: estimate_only=TRUE and FALSE give same point estimates", {
  d    <- make_lmm_data()
  res1 <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = TRUE)
  res2 <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)
  expect_equal(res1$b[1:3], res2$b[1:3], tolerance = 1e-6,
    label = "beta estimates match across estimate_only flag")
  expect_equal(res1$neg_loglik, res2$neg_loglik, tolerance = 1e-6)
})

test_that("gaussian_lmm: converges and returns finite estimates", {
  d   <- make_lmm_data()
  res <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)
  expect_true(res$converged)
  expect_true(all(is.finite(res$b)))
  expect_true(is.finite(res$neg_loglik))
  expect_true(all(is.finite(res$vcov)))
})

test_that("gaussian_lmm: sigma_e > 0, sigma_b > 0", {
  d   <- make_lmm_data()
  res <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)
  p   <- ncol(d$X)
  expect_true(exp(res$b[p + 1]) > 0)   # sigma_e
  expect_true(exp(res$b[p + 2]) > 0)   # sigma_b
})

test_that("gaussian_lmm: beta estimates close to lme4::lmer", {
  skip_if_not_installed("lme4")
  d   <- make_lmm_data(seed = 7L)
  res <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)

  df <- data.frame(y = d$y, x1 = d$X[,2], x2 = d$X[,3], g = factor(d$group_id))
  suppressMessages(fit4 <- lme4::lmer(y ~ x1 + x2 + (1|g), data = df, REML = FALSE))
  beta4 <- as.numeric(lme4::fixef(fit4))
  expect_equal(unname(res$b[1:3]), unname(beta4), tolerance = 0.05,
    label = "EDI beta vs lme4 fixef (Gaussian LMM)")
})

test_that("gaussian_lmm: vcov is symmetric and positive definite", {
  d    <- make_lmm_data()
  res  <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)
  V    <- res$vcov
  k    <- ncol(d$X) + 2L  # p + log_se + log_sb
  expect_equal(dim(V), c(k, k))
  expect_equal(V, t(V), tolerance = 1e-12, label = "vcov symmetric")
  expect_true(all(eigen(V, only.values = TRUE)$values > 0), label = "vcov positive definite")
})

test_that("gaussian_lmm: hessian symmetric (preallocated SM2 buffers don't corrupt)", {
  d   <- make_lmm_data()
  res <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)
  par <- as.numeric(res$b)
  H   <- as.matrix(EDI:::get_gaussian_lmm_fisher_cpp(d$X, d$y, d$group_id, par))
  expect_equal(H, t(H), tolerance = 1e-10, label = "Hessian symmetric")
  expect_true(all(eigen(H, only.values = TRUE)$values > 0),
    label = "Hessian positive definite at MLE")
})

test_that("gaussian_lmm: score ≈ 0 at convergence (SM2 buffers don't corrupt gradient)", {
  d   <- make_lmm_data()
  res <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)
  par <- as.numeric(res$b)
  sc  <- as.numeric(EDI:::get_gaussian_lmm_score_cpp(d$X, d$y, d$group_id, par))
  expect_true(max(abs(sc[1:ncol(d$X)])) < 0.01,
    label = "beta score near zero at convergence")
})

test_that("gaussian_lmm: warm start matches cold start", {
  d    <- make_lmm_data(seed = 3L)
  cold <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id, estimate_only = FALSE)
  warm <- EDI:::fast_gaussian_lmm_cpp(d$X, d$y, d$group_id,
           warm_start_params = cold$b, estimate_only = FALSE)
  expect_equal(cold$b, warm$b, tolerance = 1e-4)
  expect_equal(cold$neg_loglik, warm$neg_loglik, tolerance = 1e-4)
})

test_that("gaussian_lmm: mixed group sizes (pairs + singletons)", {
  set.seed(42)
  G_pairs <- 40L; G_sing <- 20L
  group_id <- c(rep(seq_len(G_pairs), each = 2L),
                as.integer(G_pairs + seq_len(G_sing)))
  n   <- length(group_id)
  X   <- cbind(1, rnorm(n))
  y   <- as.numeric(X %*% c(0.5, 0.3) + rnorm(n, 0, 0.6))
  res <- EDI:::fast_gaussian_lmm_cpp(X, y, group_id, estimate_only = FALSE)
  expect_true(res$converged)
  expect_true(all(is.finite(res$vcov)))
})
