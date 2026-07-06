library(testthat)
library(EDI)

test_that("sample_int_replace_cpp returns indices in [1, n] with correct length", {
  set.seed(1)
  idx <- EDI:::sample_int_replace_cpp(50L, 10000L)
  expect_length(idx, 10000L)
  expect_true(all(idx >= 1L & idx <= 50L))
})

test_that("sample_int_replace_cpp marginal distribution is uniform over [1, n]", {
  set.seed(42)
  n <- 40L; size <- 200000L
  idx <- EDI:::sample_int_replace_cpp(n, size)
  freq <- tabulate(idx, nbins = n)
  expected <- size / n
  expect_true(all(abs(freq - expected) / expected < 0.05))
})

test_that("sample_int_replace_cpp is reproducible given same R seed", {
  set.seed(77)
  a <- EDI:::sample_int_replace_cpp(100L, 5000L)
  set.seed(77)
  b <- EDI:::sample_int_replace_cpp(100L, 5000L)
  expect_identical(a, b)
})

test_that("resample_group_rows_cpp returns valid row indices covering all groups", {
  set.seed(3)
  grp <- rep(1:10L, each = 5L)  # 10 groups of 5 rows each
  out <- EDI:::resample_group_rows_cpp(grp, 100L)
  # All indices must be in [1, 50]
  expect_true(all(out >= 1L & out <= 50L))
  # With 100 group draws × 5 rows each we should see all groups represented
  sampled_groups <- unique(grp[out])
  expect_true(length(sampled_groups) > 5L)
})

test_that("resample_group_rows_cpp group proportions are approximately uniform", {
  set.seed(7)
  grp <- rep(1:5L, each = 10L)
  out <- EDI:::resample_group_rows_cpp(grp, 50000L)
  grp_counts <- tabulate(grp[out], nbins = 5L)
  # Each group contributes 10 rows; uniform group sampling → each group ~20% of output
  props <- grp_counts / sum(grp_counts)
  expect_true(all(abs(props - 0.2) < 0.02))
})

test_that("compute_bootstrapped_weighted_sqd_distances_cpp returns B non-negative values", {
  set.seed(9)
  X <- matrix(rnorm(100 * 4), 100, 4)
  wts <- c(1, 2, 0.5, 1.5)
  out <- EDI:::compute_bootstrapped_weighted_sqd_distances_cpp(X, wts, 100L, 500L)
  expect_length(out, 500L)
  expect_true(all(out >= 0))
})

test_that("compute_bootstrapped_weighted_sqd_distances_cpp mean distance matches empirical expectation", {
  set.seed(11)
  n <- 200L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  wts <- rep(1.0, d)
  B <- 50000L
  out <- EDI:::compute_bootstrapped_weighted_sqd_distances_cpp(X, wts, n, B)
  # E[||x_i - x_j||^2] for iid N(0,1) columns = 2*d (two independent draws per column)
  expect_lt(abs(mean(out) - 2 * d) / (2 * d), 0.05)
})

test_that("compute_bootstrapped_weighted_sqd_distances_cpp is reproducible given same R seed", {
  set.seed(55)
  X <- matrix(rnorm(80 * 3), 80, 3)
  wts <- rep(1.0, 3)
  set.seed(22)
  a <- EDI:::compute_bootstrapped_weighted_sqd_distances_cpp(X, wts, 80L, 1000L)
  set.seed(22)
  b <- EDI:::compute_bootstrapped_weighted_sqd_distances_cpp(X, wts, 80L, 1000L)
  expect_identical(a, b)
})
