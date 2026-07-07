cloglog_reference <- function(z) {
  F <- numeric(length(z))
  f <- numeric(length(z))
  fp <- numeric(length(z))

  high <- z > 5
  low <- z < -37
  mid <- !(high | low)

  F[high] <- 1
  F[low] <- 0
  f[high | low] <- 0
  fp[high | low] <- 0

  ez <- exp(z[mid])
  exp_neg_ez <- exp(-ez)
  F[mid] <- 1 - exp_neg_ez
  f[mid] <- exp(z[mid] - ez)
  fp[mid] <- f[mid] * (1 - ez)

  list(cdf = F, pdf = f, pdf_derivative = fp)
}

test_that("cached cloglog endpoint evaluation matches independent formulas", {
  x <- c(
    -Inf, Inf, -100, 100,
    -37 + c(-1e-10, 0, 1e-10),
    5 + c(-1e-10, 0, 1e-10),
    seq(-37, 5, length.out = 50001L)
  )

  actual <- EDI:::fast_cloglog_link_eval_cpp(x)
  reference <- cloglog_reference(x)

  expect_lte(max(abs(actual$cdf - reference$cdf)), 3e-16)
  expect_lte(max(abs(actual$pdf - reference$pdf)), 1e-15)
  expect_lte(max(abs(actual$pdf_derivative - reference$pdf_derivative)), 1e-15)
  expect_identical(actual$cdf[x == Inf], 1)
  expect_identical(actual$pdf[x == Inf], 0)
  expect_identical(actual$cdf[x == -Inf], 0)
  expect_identical(actual$pdf[x == -Inf], 0)
})

test_that("cached cloglog endpoint evaluation preserves score and Hessian", {
  skip_if_not_installed("numDeriv")
  set.seed(8901)
  n <- 160L
  X <- matrix(rnorm(n * 2L, sd = 0.4), ncol = 2L)
  y <- rep(1:4, length.out = n)
  params <- c(-1.5, -0.25, 1.15, 0.2, -0.15)

  exact_nll <- function(par) {
    alpha <- par[1:3]
    eta <- drop(X %*% par[4:5])
    upper <- numeric(n)
    lower <- numeric(n)

    upper[y == 4L] <- 1
    has_upper <- y < 4L
    upper[has_upper] <- cloglog_reference(alpha[y[has_upper]] + eta[has_upper])$cdf

    lower[y == 1L] <- 0
    has_lower <- y > 1L
    lower[has_lower] <- cloglog_reference(alpha[y[has_lower] - 1L] + eta[has_lower])$cdf

    -sum(log(pmax(upper - lower, 1e-12)))
  }

  score <- EDI:::get_ordinal_cloglog_regression_score_cpp(X, y, params)
  hessian <- EDI:::get_ordinal_cloglog_regression_hessian_cpp(X, y, params)

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
