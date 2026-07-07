test_that("bootstrap_m_indices_cpp is reproducible and preserves sampling units", {
  m_vec <- as.integer(c(0, 0, rep(seq_len(4L), each = 2L)))
  i_reservoir <- as.integer(c(1, 2))

  set.seed(8601)
  actual <- EDI:::bootstrap_m_indices_cpp(
    m_vec, i_reservoir, length(i_reservoir), 4L, 5000L
  )
  set.seed(8601)
  repeated <- EDI:::bootstrap_m_indices_cpp(
    m_vec, i_reservoir, length(i_reservoir), 4L, 5000L
  )

  expect_identical(actual, repeated)
  expect_equal(dim(actual), c(5000L, 10L))
  expect_true(all(actual[, 1:2] %in% i_reservoir))

  first_in_pair <- actual[, seq.int(3L, 9L, by = 2L)]
  second_in_pair <- actual[, seq.int(4L, 10L, by = 2L)]
  expect_true(all(second_in_pair == first_in_pair + 1L))
  expect_true(all(first_in_pair %in% c(3L, 5L, 7L, 9L)))

  reservoir_counts <- tabulate(actual[, 1:2], nbins = 2L)
  expect_true(all(abs(reservoir_counts - 5000L) / 5000L < 0.08))

  sampled_pair_ids <- (as.vector(first_in_pair) - 1L) %/% 2L
  pair_counts <- tabulate(sampled_pair_ids, nbins = 4L)
  expect_true(all(abs(pair_counts - 5000L) / 5000L < 0.08))
})

test_that("bootstrap_m_indices_cpp handles reservoir-only and pair-only samples", {
  set.seed(8602)
  reservoir_only <- EDI:::bootstrap_m_indices_cpp(
    integer(0), as.integer(11:13), 3L, 0L, 25L
  )
  expect_equal(dim(reservoir_only), c(25L, 3L))
  expect_true(all(reservoir_only %in% 11:13))

  set.seed(8603)
  pair_only <- EDI:::bootstrap_m_indices_cpp(
    as.integer(rep(1:3, each = 2L)), integer(0), 0L, 3L, 25L
  )
  expect_equal(dim(pair_only), c(25L, 6L))
  expect_true(all(pair_only[, seq.int(2L, 6L, by = 2L)] ==
                    pair_only[, seq.int(1L, 5L, by = 2L)] + 1L))
})
