library(testthat)
library(EDI)

test_that("resample_group_rows_cpp matches explicit group expansion", {
	group_id = c(1L, 1L, 2L, 3L, 3L, 3L)

	set.seed(20260510)
	actual = EDI:::resample_group_rows_cpp(group_id, 3L)

	set.seed(20260510)
	sampled_groups = EDI:::sample_int_replace_cpp(3L, 3L)
	expected = unlist(lapply(sampled_groups, function(g) which(group_id == g)), use.names = FALSE)

	expect_identical(actual, as.integer(expected))
})

test_that("resample_group_rows_cpp validates group ids", {
	expect_error(
		EDI:::resample_group_rows_cpp(c(1L, 3L, 3L), 2L),
		"consecutive positive integers"
	)
	expect_error(
		EDI:::resample_group_rows_cpp(c(1L, 0L, 1L), 2L),
		"positive integers"
	)
})

test_that("DesignFixedBlocking block bootstrap uses group resampling helper semantics", {
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
	expected = EDI:::resample_group_rows_cpp(as.integer(group_id), length(unique(group_id)))

	expect_identical(actual, expected)
})
