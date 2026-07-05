# Extracted from test-group-bootstrap-rcpp.R:5

# test -------------------------------------------------------------------------
group_id = c(1L, 1L, 2L, 3L, 3L, 3L)
set.seed(20260510)
actual = resample_group_rows_cpp(group_id, 3L)
