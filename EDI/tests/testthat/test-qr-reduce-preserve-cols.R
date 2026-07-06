library(testthat)
library(EDI)

test_that("qr_reduce_preserve_cols_cpp preserves required columns when possible", {
	X = cbind(
		1,
		c(0, 1, 0, 1, 0, 1),
		c(10, 11, 10, 11, 10, 11),
		c(1, 1, 1, 1, 1, 1),
		c(0, 1, 0, 1, 0, 1)
	)

	reduced = EDI:::qr_reduce_preserve_cols_cpp(X, c(1L, 2L))

	expect_true(all(c(1L, 2L) %in% reduced$keep))
	expect_equal(qr(reduced$X_reduced)$rank, ncol(reduced$X_reduced))
	expect_equal(qr(X)$rank, qr(reduced$X_reduced)$rank)
})

test_that("qr_reduce_preserve_cols_cpp drops treatment when it is linearly dependent", {
	X = cbind(
		1,
		rep(1, 6),
		c(1, 2, 3, 4, 5, 6)
	)

	reduced = EDI:::qr_reduce_preserve_cols_cpp(X, c(1L, 2L))

	expect_true(1L %in% reduced$keep)
	expect_false(2L %in% reduced$keep)
	expect_equal(qr(reduced$X_reduced)$rank, ncol(reduced$X_reduced))
})
