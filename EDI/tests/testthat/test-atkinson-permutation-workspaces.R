library(testthat)
library(EDI)

test_that("Atkinson permutation workspace reuse preserves seeded output", {
	set.seed(80)
	n <- 40L
	p <- 4L
	X <- matrix(rnorm(n * p), n, p)

	set.seed(801)
	first <- EDI:::generate_permutations_atkinson_cpp(X, n, p, 0.5, 60L)
	set.seed(801)
	second <- EDI:::generate_permutations_atkinson_cpp(X, n, p, 0.5, 60L)

	expect_identical(first$w_mat, second$w_mat)
	expect_equal(dim(first$w_mat), c(n, 60L))
	expect_true(all(first$w_mat %in% 0:1))
	expect_null(first$m_mat)
})

test_that("Atkinson permutation workspace handles rank deficiency", {
	set.seed(802)
	n <- 35L
	x <- rnorm(n)
	X <- cbind(constant = 1, x = x, duplicate = x, scaled = 2 * x)
	result <- EDI:::generate_permutations_atkinson_cpp(X, n, 4L, 0.5, 80L)$w_mat

	expect_equal(dim(result), c(n, 80L))
	expect_true(all(result %in% 0:1))
})

test_that("Atkinson permutations remain symmetric at probability one half", {
	set.seed(803)
	n <- 40L
	p <- 4L
	X <- matrix(rnorm(n * p), n, p)
	result <- EDI:::generate_permutations_atkinson_cpp(X, n, p, 0.5, 500L)$w_mat

	expect_lt(abs(mean(result) - 0.5), 0.03)
})
