# Extracted from test-qr-reduce-preserve-cols.R:10

# test -------------------------------------------------------------------------
X = cbind(
		1,
		c(0, 1, 0, 1, 0, 1),
		c(10, 11, 10, 11, 10, 11),
		c(1, 1, 1, 1, 1, 1),
		c(0, 1, 0, 1, 0, 1)
	)
reduced = qr_reduce_preserve_cols_cpp(X, c(1L, 2L))
