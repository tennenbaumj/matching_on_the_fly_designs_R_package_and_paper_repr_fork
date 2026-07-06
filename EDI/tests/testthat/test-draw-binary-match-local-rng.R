library(testthat)
library(EDI)

make_pairs <- function(n) matrix(c(seq(1L, n, 2L), seq(2L, n, 2L)), ncol = 2L)

test_that("draw_binary_match_assignments_cpp output is binary and correct dimensions", {
  set.seed(1)
  n <- 100L; r <- 50L
  pairs <- make_pairs(n)
  w <- EDI:::draw_binary_match_assignments_cpp(pairs, n, r, 1L)
  expect_equal(dim(w), c(n, r))
  expect_true(all(w %in% c(0, 1)))
})

test_that("draw_binary_match_assignments_cpp: each matched pair has exactly one treatment per sim", {
  set.seed(2)
  n <- 80L; r <- 200L
  pairs <- make_pairs(n)
  w <- EDI:::draw_binary_match_assignments_cpp(pairs, n, r, 2L)
  num_pairs <- nrow(pairs)
  for (b in seq_len(r)) {
    pair_sums <- w[pairs[, 1], b] + w[pairs[, 2], b]
    expect_true(all(pair_sums == 1L),
      label = sprintf("sim %d: every pair sums to 1", b))
  }
})

test_that("draw_binary_match_assignments_cpp: Bernoulli marginals ≈ 0.5", {
  set.seed(3)
  n <- 100L; r <- 5000L
  pairs <- make_pairs(n)
  w <- EDI:::draw_binary_match_assignments_cpp(pairs, n, r, 2L)
  subject_means <- rowMeans(w)
  expect_true(all(abs(subject_means - 0.5) < 0.03),
    label = "all subject marginals within 0.03 of 0.5")
})

test_that("draw_binary_match_assignments_cpp: each subject assigned in every sim", {
  set.seed(4)
  n <- 60L; r <- 100L
  pairs <- make_pairs(n)
  w <- EDI:::draw_binary_match_assignments_cpp(pairs, n, r, 2L)
  # w contains 0/1 for all subjects; no NAs, no zeros from unassigned
  expect_false(anyNA(w))
  # Each subject should not be all-zero across sims (would be nearly impossible)
  expect_true(all(rowSums(w) > 0))
})

test_that("draw_binary_match_assignments_cpp is reproducible given same R seed", {
  set.seed(5)
  n <- 80L; r <- 100L
  pairs <- make_pairs(n)
  set.seed(42)
  w1 <- EDI:::draw_binary_match_assignments_cpp(pairs, n, r, 2L)
  set.seed(42)
  w2 <- EDI:::draw_binary_match_assignments_cpp(pairs, n, r, 2L)
  expect_identical(w1, w2)
})

test_that("draw_binary_match_assignments_cpp: within-pair treatment is symmetric across sims", {
  set.seed(6)
  n <- 100L; r <- 2000L
  pairs <- make_pairs(n)
  w <- EDI:::draw_binary_match_assignments_cpp(pairs, n, r, 2L)
  # For each pair, fraction of sims where first member is treated should be near 0.5
  frac_first <- rowMeans(w[pairs[, 1], , drop = FALSE])
  expect_true(all(abs(frac_first - 0.5) < 0.05),
    label = "within-pair first-member treatment fraction near 0.5")
})
