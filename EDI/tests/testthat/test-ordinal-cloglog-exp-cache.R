test_that("fast_cloglog_link_eval_cpp fused CDF/PDF/PDF' match closed-form cloglog", {
  # F(z) = 1 - exp(-exp(z)), f(z) = exp(z - exp(z)), f'(z) = f*(1 - exp(z))
  z_vals <- c(-10, -3, -1, 0, 0.5, 1, 2, 3, 5, 6)
  res <- EDI:::fast_cloglog_link_eval_cpp(z_vals)

  F_ref  <- ifelse(z_vals > 5, 1, ifelse(z_vals < -37, 0, 1 - exp(-exp(z_vals))))
  f_ref  <- ifelse(z_vals > 5 | z_vals < -37, 0, exp(z_vals - exp(z_vals)))
  fp_ref <- ifelse(z_vals > 5 | z_vals < -37, 0, f_ref * (1 - exp(z_vals)))

  expect_equal(res$cdf,            F_ref,  tolerance = 1e-12)
  expect_equal(res$pdf,            f_ref,  tolerance = 1e-12)
  expect_equal(res$pdf_derivative, fp_ref, tolerance = 1e-12)
})

test_that("cloglog gradient matches numerical gradient after exp-caching refactor", {
  set.seed(7L)
  n <- 200L; p <- 3L; K <- 4L
  X  <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y  <- sample(seq_len(K), n, replace = TRUE)
  fit    <- EDI:::fast_ordinal_cloglog_regression_cpp(X, y)
  params <- fit$params

  score <- as.numeric(EDI:::get_ordinal_cloglog_regression_score_cpp(X, y, params))

  # R-side NLL for numerical gradient
  cloglog_cdf <- function(z) ifelse(z > 5, 1, ifelse(z < -37, 0, 1 - exp(-exp(z))))
  nll_r <- function(par) {
    n_alpha <- K - 1L
    alpha <- par[seq_len(n_alpha)]
    beta  <- par[seq_len(p) + n_alpha]
    eta   <- as.numeric(X %*% beta)
    lev   <- sort(unique(y))
    total <- 0
    for (i in seq_len(n)) {
      yi   <- which(lev == y[i])
      pu   <- if (yi == K) 1 else cloglog_cdf(alpha[yi] + eta[i])
      pl   <- if (yi == 1) 0 else cloglog_cdf(alpha[yi - 1L] + eta[i])
      prob <- max(1e-12, pu - pl)
      total <- total - log(prob)
    }
    total
  }

  eps      <- 1e-6
  num_grad <- numeric(length(params))
  for (j in seq_along(params)) {
    pp <- pm <- params; pp[j] <- pp[j] + eps; pm[j] <- pm[j] - eps
    num_grad[j] <- (nll_r(pp) - nll_r(pm)) / (2 * eps)
  }

  # score_cpp returns -gradient(nll), so score ≈ -num_grad
  expect_equal(score, -num_grad, tolerance = 1e-4)
})

test_that("cloglog with-var fit produces finite params and SE after exp-caching refactor", {
  set.seed(55L)
  n <- 300L; p <- 3L; K <- 5L
  X   <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y   <- sample(seq_len(K), n, replace = TRUE)
  fit <- EDI:::fast_ordinal_cloglog_regression_with_var_cpp(X, y)

  expect_true(fit$converged)
  expect_true(all(is.finite(fit$params)))
  expect_true(all(is.finite(fit$se)) && all(fit$se > 0))
})
