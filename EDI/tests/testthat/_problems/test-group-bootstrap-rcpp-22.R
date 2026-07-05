# Extracted from test-group-bootstrap-rcpp.R:22

# test -------------------------------------------------------------------------
expect_error(
		resample_group_rows_cpp(c(1L, 3L, 3L), 2L),
		"consecutive positive integers"
	)
expect_error(
		resample_group_rows_cpp(c(1L, 0L, 1L), 2L),
		"positive integers"
	)
