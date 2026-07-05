# Extracted from test-qr-reduce-preserve-cols.R:24

# test -------------------------------------------------------------------------
X = cbind(
		1,
		rep(1, 6),
		c(1, 2, 3, 4, 5, 6)
	)
reduced = qr_reduce_preserve_cols_cpp(X, c(1L, 2L))
