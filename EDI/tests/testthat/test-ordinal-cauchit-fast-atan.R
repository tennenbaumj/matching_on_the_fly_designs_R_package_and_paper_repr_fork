test_that("fast cauchit atan meets a near-double-precision error guard", {
  tan_pi_over_8 <- sqrt(2) - 1
  tan_3pi_over_8 <- sqrt(2) + 1
  boundaries <- c(
    -Inf, Inf, -100, 100, -1, 1, -0, 0,
    outer(
      c(-tan_3pi_over_8, -tan_pi_over_8, tan_pi_over_8, tan_3pi_over_8),
      c(1 - 4 * .Machine$double.eps, 1, 1 + 4 * .Machine$double.eps)
    )
  )
  dense <- seq(-100, 100, length.out = 400001L)
  tails <- 10^seq(-300, 300, length.out = 20000L)
  x <- c(boundaries, dense, -tails, tails)

  actual <- EDI:::fast_atan_cauchit_cpp(x)
  reference <- atan(x)
  absolute_error <- abs(actual - reference)

  expect_lte(max(absolute_error, na.rm = TRUE), 5e-16)
  expect_lte(max(abs((0.5 + actual / pi) - pcauchy(x)), na.rm = TRUE), 5e-16)
  expect_identical(actual[c(1L, 2L)], reference[c(1L, 2L)])
  expect_identical(1 / actual[[7L]], -Inf)

  special <- EDI:::fast_atan_cauchit_cpp(c(NA_real_, NaN))
  expect_true(is.na(special[[1L]]))
  expect_true(is.nan(special[[2L]]))
})

test_that("fast cauchit CDF preserves exact-likelihood derivatives", {
  skip_if_not_installed("numDeriv")
  set.seed(8801)
  X <- matrix(rnorm(160), ncol = 2L)
  y <- rep(1:4, length.out = nrow(X))
  params <- c(-1.4, -0.1, 1.2, 0.35, -0.2)

  exact_nll <- function(par) {
    alpha <- par[1:3]
    eta <- drop(X %*% par[4:5])
    upper <- ifelse(y == 4L, 1, pcauchy(alpha[pmin(y, 3L)] - eta))
    lower <- ifelse(y == 1L, 0, pcauchy(alpha[pmax(y - 1L, 1L)] - eta))
    -sum(log(pmax(upper - lower, 1e-12)))
  }

  score <- EDI:::get_ordinal_cauchit_regression_score_cpp(X, y, params)
  hessian <- EDI:::get_ordinal_cauchit_regression_hessian_cpp(X, y, params)

  expect_equal(
    as.numeric(score),
    as.numeric(-numDeriv::grad(exact_nll, params)),
    tolerance = 1e-7
  )
  expect_equal(
    unname(hessian),
    unname(-numDeriv::hessian(exact_nll, params)),
    tolerance = 1e-5
  )
})
