reference_objective_vals <- function(X, indicTs, objective, inv_cov_X = NULL) {
  sum_all <- colSums(X)
  sd_all <- apply(X, 2L, stats::sd)
  apply(indicTs, 1L, function(w) {
    treated <- w == 1L
    nT <- sum(treated)
    nC <- length(w) - nT
    sum_T <- colSums(X[treated, , drop = FALSE])
    diff <- sum_T / nT - (sum_all - sum_T) / nC
    if (objective == "abs_sum_diff") {
      sum(abs(diff / sd_all))
    } else {
      as.numeric(crossprod(diff, inv_cov_X %*% diff))
    }
  })
}

balanced_indicators <- function(r, n, nT) {
  t(replicate(r, as.integer(seq_len(n) %in% sample.int(n, nT))))
}

test_that("compute_objective_vals_cpp matches independent R objective values", {
  set.seed(10701)
  n <- 30L
  p <- 4L
  r <- 400L
  X <- matrix(rnorm(n * p), n, p)
  indicTs <- balanced_indicators(r, n, nT = 11L)

  expect_equal(
    EDI:::compute_objective_vals_cpp(X, indicTs, "abs_sum_diff"),
    reference_objective_vals(X, indicTs, "abs_sum_diff"),
    tolerance = 1e-12
  )

  inv_cov_X <- solve(stats::cov(X) + diag(0.1, p))
  expect_equal(
    EDI:::compute_objective_vals_cpp(X, indicTs, "mahal_dist", inv_cov_X),
    reference_objective_vals(X, indicTs, "mahal_dist", inv_cov_X),
    tolerance = 1e-12
  )
})

test_that("compute_objective_vals_cpp validates rerandomization objective inputs", {
  X <- matrix(rnorm(20L), 10L, 2L)
  indicTs <- balanced_indicators(5L, 10L, nT = 4L)

  expect_error(
    EDI:::compute_objective_vals_cpp(X, indicTs[, -1L], "abs_sum_diff"),
    "ncol\\(indicTs\\) == nrow\\(X\\)"
  )
  expect_error(
    EDI:::compute_objective_vals_cpp(X, indicTs, "mahal_dist"),
    "inv_cov_X required"
  )
  expect_error(
    EDI:::compute_objective_vals_cpp(X, indicTs, "unknown"),
    "objective must be"
  )
})
