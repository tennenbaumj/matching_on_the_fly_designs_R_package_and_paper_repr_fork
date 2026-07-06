library(testthat)
library(EDI)

# Shared simulation fixture
make_poisson_glmm_data <- function(seed = 1L, G = 30L, nj = 4L) {
  set.seed(seed)
  n        <- G * nj
  group_id <- as.integer(rep(seq_len(G), each = nj))
  X        <- cbind(1, rnorm(n), rnorm(n))
  b_true   <- c(0.3, 0.4, -0.2)
  u_true   <- rnorm(G, 0, 0.3)
  eta      <- X %*% b_true + u_true[group_id]
  y        <- as.numeric(rpois(n, exp(eta)))
  list(X = X, y = y, group_id = group_id, n = n, G = G)
}

test_that("poisson_glmm_cpp: estimate_only=TRUE and FALSE give same point estimates", {
  d    <- make_poisson_glmm_data()
  res1 <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = TRUE)
  res2 <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)

  expect_equal(res1$b,         res2$b,         tolerance = 1e-6)
  expect_equal(res1$log_sigma, res2$log_sigma,  tolerance = 1e-6)
  expect_equal(res1$neg_loglik, res2$neg_loglik, tolerance = 1e-6)
})

test_that("poisson_glmm_cpp: converges and returns finite estimates", {
  d   <- make_poisson_glmm_data()
  res <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)

  expect_true(res$converged)
  expect_true(all(is.finite(res$b)))
  expect_true(is.finite(res$log_sigma))
  expect_true(is.finite(res$neg_loglik))
  expect_true(all(is.finite(res$vcov)))
})

test_that("poisson_glmm_cpp: sigma > 0 (log_sigma is finite and sensible)", {
  d   <- make_poisson_glmm_data()
  res <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)
  expect_true(exp(res$log_sigma) > 0)
  expect_true(exp(res$log_sigma) < 10)  # reasonable range
})

test_that("poisson_glmm_cpp: beta estimates close to lme4::glmer", {
  skip_if_not_installed("lme4")
  d    <- make_poisson_glmm_data(seed = 7L)
  res  <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)

  df   <- data.frame(y = d$y, x1 = d$X[, 2], x2 = d$X[, 3], g = factor(d$group_id))
  suppressMessages(
    fit4 <- lme4::glmer(y ~ x1 + x2 + (1 | g), data = df, family = poisson)
  )
  beta4 <- as.numeric(lme4::fixef(fit4))

  expect_equal(res$b, beta4, tolerance = 0.05,
    label = "EDI beta vs lme4 fixef (Poisson GLMM)")
})

test_that("poisson_glmm_cpp: variance-covariance matrix is symmetric and positive-definite", {
  d    <- make_poisson_glmm_data()
  res  <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)
  V    <- res$vcov
  p1   <- ncol(d$X) + 1L  # p + 1 (includes log_sigma)

  expect_equal(dim(V), c(p1, p1))
  expect_equal(V, t(V), tolerance = 1e-12, label = "vcov is symmetric")
  evals <- eigen(V, only.values = TRUE)$values
  expect_true(all(evals > 0), label = "vcov is positive definite")
})

test_that("poisson_glmm_cpp: score ≈ 0 at convergence (preallocated buffers don't corrupt gradient)", {
  d    <- make_poisson_glmm_data()
  res  <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)
  par  <- as.numeric(c(res$b, res$log_sigma))
  score <- as.numeric(EDI:::get_poisson_glmm_score_cpp(d$X, d$y, d$group_id, par))
  # Score should be near zero at convergence (ignoring soft barrier penalty on log_sigma)
  expect_true(max(abs(score[seq_len(ncol(d$X))])) < 0.05,
    label = "beta score components near zero at convergence")
})

test_that("poisson_glmm_cpp: hessian symmetric at convergence", {
  d    <- make_poisson_glmm_data()
  res  <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)
  par  <- as.numeric(c(res$b, res$log_sigma))
  H    <- as.matrix(EDI:::get_poisson_glmm_hessian_cpp(d$X, d$y, d$group_id, par))
  expect_equal(H, t(H), tolerance = 1e-10, label = "Hessian is symmetric")
})

test_that("poisson_glmm_cpp: warm-start matches cold-start estimates", {
  d     <- make_poisson_glmm_data(seed = 3L)
  cold  <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L, estimate_only = FALSE)
  warm  <- EDI:::fast_poisson_glmm_cpp(d$X, d$y, d$group_id, 1L,
             warm_start_params = c(cold$b, cold$log_sigma), estimate_only = FALSE)

  expect_equal(cold$b,         warm$b,         tolerance = 1e-4)
  expect_equal(cold$log_sigma, warm$log_sigma,  tolerance = 1e-4)
  expect_equal(cold$neg_loglik, warm$neg_loglik, tolerance = 1e-4)
})
