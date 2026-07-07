robust_ssq_reference <- function(X, y, fit, j, method = "MM", c = 1.345,
                                 fixed_idx = integer()) {
  beta <- as.numeric(fit$coefficients)
  scale <- fit$scale
  n <- nrow(X)
  p <- ncol(X)
  free_idx <- setdiff(seq_len(p), fixed_idx)
  free_j <- match(j, free_idx)
  if (is.na(free_j)) {
    return(NA_real_)
  }

  r <- as.numeric(y - X %*% beta)
  u <- r / scale
  if (method == "M") {
    abs_u <- abs(u)
    psi <- ifelse(abs_u <= c, r, c * scale * ifelse(u > 0, 1, -1))
    sum_psi_prime <- sum(abs_u <= c)
  } else {
    c_b <- 4.685
    abs_u <- abs(u)
    u_scaled_sq <- (u / c_b)^2
    tmp <- 1 - u_scaled_sq
    psi <- ifelse(abs_u <= c_b, r * tmp^2, 0)
    sum_psi_prime <- sum(ifelse(abs_u <= c_b, tmp * (1 - 5 * u_scaled_sq), 0))
  }

  m <- sum_psi_prime / n
  factor <- (n / (n - length(free_idx))) * sum(psi^2) / (n * m^2)
  XtX <- crossprod(X[, free_idx, drop = FALSE])
  factor * solve(XtX)[free_j, free_j]
}

test_that("robust regression cached XtX gives reference ssq_b_j", {
  set.seed(5401)
  n <- 240L
  x1 <- rnorm(n)
  x2 <- x1 * 0.4 + rnorm(n)
  x3 <- rnorm(n)
  x4 <- x2 * -0.3 + rnorm(n)
  X <- cbind(1, x3, x1, x4, x2)
  beta <- c(0.5, -0.2, 0.35, 0.1, -0.25)
  y <- as.numeric(X %*% beta + rt(n, df = 4))
  y[seq(7L, n, by = 31L)] <- y[seq(7L, n, by = 31L)] + 8

  fit_mm <- EDI:::fast_robust_regression_cpp(X, y, method = "MM", j = 3L)
  expect_true(isTRUE(fit_mm$converged))
  expect_equal(
    fit_mm$ssq_b_j,
    robust_ssq_reference(X, y, fit_mm, j = 3L, method = "MM"),
    tolerance = 1e-10
  )

  fit_m <- EDI:::fast_robust_regression_cpp(X, y, method = "M", j = 5L, c = 1.5)
  expect_true(isTRUE(fit_m$converged))
  expect_equal(
    fit_m$ssq_b_j,
    robust_ssq_reference(X, y, fit_m, j = 5L, method = "M", c = 1.5),
    tolerance = 1e-10
  )
})

test_that("robust regression cached XtX respects fixed parameters", {
  set.seed(5402)
  n <- 180L
  X <- cbind(1, matrix(rnorm(n * 3L), n, 3L))
  beta <- c(0.25, 0.4, -0.2, 0.15)
  y <- as.numeric(X %*% beta + rnorm(n))
  fixed_idx <- 2L
  fixed_values <- 0.4

  fit <- EDI:::fast_robust_regression_cpp(
    X,
    y,
    method = "MM",
    j = 4L,
    fixed_idx = fixed_idx,
    fixed_values = fixed_values
  )

  expect_equal(fit$coefficients[fixed_idx], fixed_values, tolerance = 0)
  expect_equal(
    fit$ssq_b_j,
    robust_ssq_reference(X, y, fit, j = 4L, method = "MM", fixed_idx = fixed_idx),
    tolerance = 1e-10
  )
})
