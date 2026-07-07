library(testthat)
library(EDI)

make_pocock_simon_levels <- function(n, levels_per_factor) {
	offset <- c(0L, cumsum(head(levels_per_factor, -1L)))
	vapply(seq_along(levels_per_factor), function(j) {
		sample.int(levels_per_factor[j], n, replace = TRUE) + offset[j]
	}, integer(n))
}

test_that("Pocock-Simon permutation buffers preserve seeded output", {
	set.seed(81)
	n <- 80L
	levels_per_factor <- c(2L, 3L, 2L)
	level_rows <- make_pocock_simon_levels(n, levels_per_factor)
	weights <- c(1, 0.75, 1.25)

	set.seed(811)
	first <- EDI:::generate_permutations_pocock_simon_cpp(
		level_rows, sum(levels_per_factor), weights, 0.7, 0.5, 120L
	)
	set.seed(811)
	second <- EDI:::generate_permutations_pocock_simon_cpp(
		level_rows, sum(levels_per_factor), weights, 0.7, 0.5, 120L
	)

	expect_identical(first$w_mat, second$w_mat)
	expect_equal(dim(first$w_mat), c(n, 120L))
	expect_true(all(first$w_mat %in% 0:1))
	expect_null(first$m_mat)
	expect_lt(abs(mean(first$w_mat) - 0.5), 0.03)
})

test_that("Pocock-Simon permutation indices are validated once", {
	valid <- matrix(as.integer(c(1, 2, 3, 4)), nrow = 2L)
	expect_error(
		EDI:::generate_permutations_pocock_simon_cpp(valid - 1L, 4L, c(1, 1), 0.7, 0.5, 2L),
		"outside 1..num_levels_total"
	)
	expect_error(
		EDI:::generate_permutations_pocock_simon_cpp(valid, 4L, 1, 0.7, 0.5, 2L),
		"weights length"
	)
})
