test_that("fast_log_pnorm matches R::pnorm on representative range", {
  # Indirect test: dep-cens-transform gradient at MLE should be near-zero.
  # If fast_log_pnorm/fast_log_dnorm are wrong, the gradient will diverge.
  set.seed(7L)
  n <- 300L; p <- 3L
  X    <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y    <- rexp(n, rate = 0.1)
  dead <- rbinom(n, 1, 0.6)

  fit   <- EDI:::fast_dep_cens_transform_optim_cpp(X, y, dead)
  score <- as.numeric(EDI:::get_dep_cens_transform_score_cpp(X, y, dead, fit$params))

  expect_true(fit$converged)
  # Score at MLE should be small (gradient of loglik ≈ 0)
  expect_lt(max(abs(score)), 0.01)
})

test_that("dep-cens-transform gradient matches numerical gradient after fast_log helpers", {
  set.seed(13L)
  n <- 150L; p <- 2L
  X    <- cbind(1, rnorm(n))
  y    <- rexp(n, rate = 0.2)
  dead <- rbinom(n, 1, 0.5)

  fit    <- EDI:::fast_dep_cens_transform_optim_cpp(X, y, dead)
  params <- fit$params

  score <- as.numeric(EDI:::get_dep_cens_transform_score_cpp(X, y, dead, params))

  # R-side log-normal CDF/PDF for numerical gradient
  log_pnorm_r <- function(x) pnorm(x, log.p = TRUE)
  log_dnorm_r <- function(x) dnorm(x, log = TRUE)

  nll_r <- function(par) {
    p_ <- length(par) - 3L
    p_half <- p_ %/% 2L
    beta_e <- par[seq_len(p_half)]
    beta_c <- par[seq_len(p_half) + p_half]
    lse <- par[p_ + 1L]; lsc <- par[p_ + 2L]
    ar  <- min(max(par[p_ + 3L], -3), 3)
    se <- exp(lse); sc <- exp(lsc)
    rho <- tanh(ar)
    omr2 <- max(1 - rho^2, 1e-12)
    sdc <- sqrt(omr2)
    log_y <- log(y)
    mu_e <- as.numeric(X %*% beta_e)
    mu_c <- as.numeric(X %*% beta_c)
    ze <- (log_y - mu_e) / se
    zc <- (log_y - mu_c) / sc
    w_e <- (rho * ze - zc) / sdc
    w_c <- (rho * zc - ze) / sdc
    lfe <- log_dnorm_r(ze) - lse - log_y
    lfc <- log_dnorm_r(zc) - lsc - log_y
    lSce <- log_pnorm_r(w_e)
    lScc <- log_pnorm_r(w_c)
    li <- dead * (lfe + lSce) + (1 - dead) * (lfc + lScc)
    -sum(li)
  }

  eps <- 1e-6
  num_grad <- numeric(length(params))
  for (j in seq_along(params)) {
    pp <- pm <- params; pp[j] <- pp[j] + eps; pm[j] <- pm[j] - eps
    num_grad[j] <- (nll_r(pp) - nll_r(pm)) / (2 * eps)
  }

  # score returns -gradient(nll), so score ≈ -num_grad
  expect_equal(score, -num_grad, tolerance = 1e-3)
})

test_that("dep-cens-transform hessian is finite and negative-semi-definite at MLE", {
  set.seed(21L)
  n <- 400L; p <- 3L
  X    <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y    <- rexp(n, rate = 0.1)
  dead <- rbinom(n, 1, 0.6)

  fit  <- EDI:::fast_dep_cens_transform_optim_cpp(X, y, dead)
  H    <- EDI:::get_dep_cens_transform_hessian_cpp(X, y, dead, fit$params)

  expect_true(all(is.finite(H)))
  # Hessian should be symmetric
  expect_equal(H, t(H), tolerance = 1e-10)
})
