library(testthat)
library(EDI)

test_that("generate_permutations_spbr_cpp: each subject gets exactly one assignment per simulation", {
  set.seed(1)
  strata_keys <- as.character(rep(paste0("S", 1:5), each = 20L))
  n <- length(strata_keys); nsim <- 100L
  w <- EDI:::generate_permutations_spbr_cpp(strata_keys, 4L, 0.5, nsim)$w_mat
  expect_equal(dim(w), c(n, nsim))
  expect_true(all(w %in% c(0L, 1L)))
})

test_that("generate_permutations_spbr_cpp: within each stratum each block sums to n_T_block", {
  set.seed(2)
  block_size <- 4L; prob_T <- 0.5
  n_T_block  <- round(block_size * prob_T)
  strata_keys <- as.character(rep(paste0("S", 1:3), each = 40L))
  nsim <- 200L
  w <- EDI:::generate_permutations_spbr_cpp(strata_keys, block_size, prob_T, nsim)$w_mat

  for (s in 1:3) {
    rows_s <- which(strata_keys == paste0("S", s))  # 40 subjects
    n_blocks <- length(rows_s) / block_size
    for (b in seq_len(nsim)) {
      col <- w[rows_s, b]
      block_sums <- vapply(seq_len(n_blocks), function(k) sum(col[((k-1)*block_size+1):(k*block_size)]), 0L)
      expect_true(all(block_sums == n_T_block),
        label = sprintf("stratum %d sim %d: block sums should all equal %d", s, b, n_T_block))
    }
  }
})

test_that("generate_permutations_spbr_cpp: grand mean ≈ prob_T across simulations", {
  set.seed(3)
  strata_keys <- as.character(rep(paste0("S", 1:4), each = 100L))
  prob_T <- 0.5
  w <- EDI:::generate_permutations_spbr_cpp(strata_keys, 4L, prob_T, 1000L)$w_mat
  expect_lt(abs(mean(w) - prob_T), 0.01)
})

test_that("generate_permutations_spbr_cpp: works with a single stratum", {
  set.seed(4)
  strata_keys <- rep("only", 20L)
  w <- EDI:::generate_permutations_spbr_cpp(strata_keys, 4L, 0.5, 50L)$w_mat
  # Each column: 5 blocks of 4, each block summing to 2
  for (b in seq_len(50L)) {
    col <- w[, b]
    block_sums <- vapply(1:5, function(k) sum(col[((k-1)*4+1):(k*4)]), 0L)
    expect_true(all(block_sums == 2L))
  }
})

test_that("generate_permutations_spbr_cpp: many strata produce valid output", {
  set.seed(5)
  strata_keys <- as.character(rep(paste0("S", 1:20), each = 10L))
  w <- EDI:::generate_permutations_spbr_cpp(strata_keys, 4L, 0.5, 200L)$w_mat
  expect_equal(dim(w), c(200L, 200L))
  expect_true(all(w %in% c(0L, 1L)))
  expect_lt(abs(mean(w) - 0.5), 0.02)
})

test_that("generate_permutations_spbr_cpp is reproducible given same R seed", {
  set.seed(6)
  strata_keys <- as.character(rep(paste0("S", 1:3), each = 20L))
  set.seed(42)
  w1 <- EDI:::generate_permutations_spbr_cpp(strata_keys, 4L, 0.5, 100L)$w_mat
  set.seed(42)
  w2 <- EDI:::generate_permutations_spbr_cpp(strata_keys, 4L, 0.5, 100L)$w_mat
  expect_identical(w1, w2)
})
