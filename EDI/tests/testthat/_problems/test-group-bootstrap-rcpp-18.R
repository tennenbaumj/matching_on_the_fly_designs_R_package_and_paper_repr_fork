# Extracted from test-group-bootstrap-rcpp.R:18

# test -------------------------------------------------------------------------
expect_error(
		resample_group_rows_cpp(c(1L, 3L, 3L), 2L),
		"consecutive positive integers"
	)
