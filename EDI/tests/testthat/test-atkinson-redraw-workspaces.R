library(testthat)
library(EDI)

test_that("Atkinson redraw workspace reuse preserves seeded assignments", {
	set.seed(83)
	n <- 60L
	p <- 4L
	X <- matrix(rnorm(n * p), n, p)

	set.seed(831)
	first <- EDI:::atkinson_redraw_batch_cpp(X, n, p, 0.5)
	set.seed(831)
	second <- EDI:::atkinson_redraw_batch_cpp(X, n, p, 0.5)

	expect_identical(first, second)
	expect_length(first, n)
	expect_true(all(first %in% 0:1))
})

test_that("Atkinson redraw workspace handles rank-deficient covariates", {
	set.seed(832)
	n <- 50L
	x <- rnorm(n)
	X <- cbind(constant = 1, x = x, duplicate = x, scaled = 2 * x)
	result <- EDI:::atkinson_redraw_batch_cpp(X, n, 4L, 0.5)

	expect_length(result, n)
	expect_true(all(is.finite(result)))
	expect_true(all(result %in% 0:1))
})

test_that("Atkinson redraw remains symmetric at probability one half", {
	set.seed(833)
	n <- 40L
	p <- 4L
	X <- matrix(rnorm(n * p), n, p)
	draws <- replicate(300L, EDI:::atkinson_redraw_batch_cpp(X, n, p, 0.5))

	expect_lt(abs(mean(draws) - 0.5), 0.03)
})
