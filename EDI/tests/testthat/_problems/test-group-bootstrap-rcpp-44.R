# Extracted from test-group-bootstrap-rcpp.R:44

# test -------------------------------------------------------------------------
des = DesignFixedBlocking$new(
		n = 6,
		response_type = "continuous",
		strata_cols = "x1",
		verbose = FALSE
	)
des$add_all_subjects_to_experiment(
		data.frame(x1 = factor(c("a", "a", "b", "b", "c", "c")))
	)
des$overwrite_all_subject_assignments(c(1, -1, 1, -1, 1, -1))
des$add_all_subject_responses(c(1, 2, 3, 4, 5, 6))
set.seed(17)
actual = des$.__enclos_env__$private$draw_bootstrap_indices("by_blocks")$i_b
strata_keys = des$.__enclos_env__$private$get_strata_keys()
group_id = match(strata_keys, unique(strata_keys))
set.seed(17)
expected = resample_group_rows_cpp(as.integer(group_id), length(unique(group_id)))
