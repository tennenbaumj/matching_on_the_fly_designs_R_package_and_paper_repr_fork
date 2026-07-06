library(testthat)
library(EDI)

test_that("generate_permutations_bernoulli_cpp marginals match Bernoulli(prob_T)", {
  set.seed(42)
  prob_T <- 0.4
  n <- 100L; nsim <- 5000L
  w <- EDI:::generate_permutations_bernoulli_cpp(n, nsim, prob_T)$w_mat
  # Row-wise means should all be ~prob_T
  row_means <- rowMeans(w)
  expect_lt(max(abs(row_means - prob_T)), 0.05)
  # Column sums: each column is iid Bernoulli(prob_T) for each subject
  grand_mean <- mean(w)
  expect_lt(abs(grand_mean - prob_T), 0.01)
})

test_that("generate_permutations_ibcrd_cpp preserves exact treatment balance", {
  set.seed(7)
  n <- 200L; nsim <- 500L; prob_T <- 0.5
  w <- EDI:::generate_permutations_ibcrd_cpp(n, nsim, prob_T)$w_mat
  col_sums <- colSums(w)
  n_T_expected <- round(n * prob_T)
  expect_true(all(col_sums == n_T_expected))
})

test_that("generate_permutations_efron_cpp: symmetric (prob_T=0.5) gives grand mean 0.5", {
  set.seed(99)
  n <- 80L; nsim <- 4000L
  w <- EDI:::generate_permutations_efron_cpp(n, nsim, 0.5, 2/3)$w_mat
  grand_mean <- mean(w)
  expect_lt(abs(grand_mean - 0.5), 0.02)
  # Efron reduces variance: per-simulation column sums should cluster tighter than Bernoulli
  col_sums_sd <- sd(colSums(w))
  bernoulli_sd <- sqrt(n * 0.5 * 0.5)  # iid Bernoulli SD
  expect_lt(col_sums_sd, bernoulli_sd)
})

test_that("generate_permutations_matching_cpp: each matched pair sums to 1", {
  set.seed(5)
  n_pairs <- 300L
  m_vec <- rep(seq_len(n_pairs), each = 2L)
  nsim <- 200L
  w <- EDI:::generate_permutations_matching_cpp(m_vec, nsim, 0.5)$w_mat
  # For each pair, the two rows should sum to 1 in every simulation
  for (pair in seq_len(n_pairs)) {
    rows <- which(m_vec == pair)
    pair_sums <- colSums(w[rows, , drop = FALSE])
    expect_true(all(pair_sums == 1L))
  }
})

test_that("bootstrap_indices_cpp returns indices in [1, n] with correct dimensions", {
  set.seed(123)
  n <- 150L; B <- 300L
  idx <- EDI:::bootstrap_indices_cpp(n, B)
  expect_equal(dim(idx), c(B, n))
  expect_true(all(idx >= 1L & idx <= n))
  # Each row should have at least some repeated indices (sampling with replacement)
  row_unique_counts <- apply(idx, 1, function(r) length(unique(r)))
  expect_true(all(row_unique_counts < n))
})

test_that("bootstrap_indices_cpp marginal distribution is uniform over [1,n]", {
  set.seed(777)
  n <- 50L; B <- 10000L
  idx <- EDI:::bootstrap_indices_cpp(n, B)
  freq <- tabulate(as.vector(idx), nbins = n)
  expected <- B * n / n  # each value expected B times
  # Chi-sq goodness-of-fit: all frequencies within 15% of expected
  expect_true(all(abs(freq - expected) / expected < 0.15))
})

test_that("bootstrap_m_indices_cpp returns valid pair indices", {
  set.seed(42)
  n_pairs <- 100L
  m_vec <- rep(seq_len(n_pairs), each = 2L)
  i_res <- integer(0)
  B <- 200L
  result <- EDI:::bootstrap_m_indices_cpp(m_vec, i_res, 0L, n_pairs, B)
  # Each pair slot should reference a valid original index (1-indexed into m_vec)
  valid_indices <- seq_len(2L * n_pairs)
  expect_true(all(result[, seq_len(2L * n_pairs)] %in% valid_indices))
})

test_that("generate_permutations_bernoulli_cpp is reproducible given same R seed", {
  set.seed(2025)
  w1 <- EDI:::generate_permutations_bernoulli_cpp(50L, 200L, 0.5)$w_mat
  set.seed(2025)
  w2 <- EDI:::generate_permutations_bernoulli_cpp(50L, 200L, 0.5)$w_mat
  expect_identical(w1, w2)
})
