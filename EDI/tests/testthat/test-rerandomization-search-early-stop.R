rerandomization_scores_reference <- function(X_raw, w_mat, objective) {
  n <- nrow(X_raw)
  X <- sweep(X_raw, 2L, colMeans(X_raw), check.margin = FALSE)

  if (objective == "mahal_dist") {
    cov_mat <- crossprod(X) / max(1L, n - 1L)
    L <- t(chol(cov_mat))
    M <- forwardsolve(L, t(X)) / n
    score <- function(d) sum(d * d)
  } else {
    inv_sd <- vapply(
      seq_len(ncol(X)),
      function(j) {
        var_j <- sum(X[, j] * X[, j]) / max(1L, n - 1L)
        if (var_j < 1e-24) 1 else 1 / sqrt(var_j)
      },
      numeric(1L)
    )
    M <- t(sweep(X, 2L, inv_sd, `*`)) / n
    score <- function(d) sum(abs(d))
  }

  apply(w_mat, 2L, function(w) score(drop(M %*% (2 * as.numeric(w) - 1))))
}

expect_valid_rerandomization_search <- function(X, w_mat, r, cutoff, objective) {
  expect_true(is.matrix(w_mat))
  expect_equal(nrow(w_mat), nrow(X))
  expect_equal(ncol(w_mat), r)
  expect_true(all(w_mat %in% c(0L, 1L)))
  expect_true(all(colSums(w_mat) == nrow(X) / 2L))
  expect_lte(max(rerandomization_scores_reference(X, w_mat, objective)), cutoff + 1e-12)
}

test_that("rerandomization_search_cpp returns valid accepted allocations with early stopping", {
  set.seed(10801)
  n <- 40L
  p <- 3L
  r <- 20L
  max_draws <- 5000L
  X <- matrix(rnorm(n * p), n, p)

  cutoff_abs <- 0.8
  abs_w <- EDI:::rerandomization_search_cpp(
    X,
    r,
    "abs_sum_diff",
    cutoff_abs,
    max_draws
  )
  expect_valid_rerandomization_search(X, abs_w, r, cutoff_abs, "abs_sum_diff")

  cutoff_mahal <- 0.5
  mahal_w <- EDI:::rerandomization_search_cpp(
    X,
    r,
    "mahal_dist",
    cutoff_mahal,
    max_draws
  )
  expect_valid_rerandomization_search(X, mahal_w, r, cutoff_mahal, "mahal_dist")
})
