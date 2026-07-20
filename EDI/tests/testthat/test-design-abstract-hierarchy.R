test_that("Design hierarchy supports both fixed and sequential designs", {
	seq_des = DesignSeqOneByOneBernoulli$new(n = 4, response_type = "continuous", verbose = FALSE)
	fixed_des = DesignFixedBernoulli$new(n = 4, response_type = "continuous", verbose = FALSE)

	expect_true(is(seq_des, "Design"))
	expect_true(is(fixed_des, "Design"))
	expect_true(is(seq_des, "DesignBlocking"))
	expect_true(is(fixed_des, "DesignBlocking"))
	expect_true(is(seq_des, "DesignMatching"))
	expect_true(is(fixed_des, "DesignMatching"))
	expect_false(is(seq_des, "DesignFixed"))
	expect_true(is(fixed_des, "DesignFixed"))
	expect_false(seq_des$is_blocking_design())
	expect_false(fixed_des$is_blocking_design())
	expect_false(seq_des$is_matching_design())
	expect_false(fixed_des$is_matching_design())
})

test_that("plain DesignFixed supports analysis but not redraw-based resampling", {
	des = EDI:::DesignFixed$new(n = 4, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = 1:4))
	des$overwrite_all_subject_assignments(c(-1, 1, -1, 1))
	des$add_all_subject_responses(c(1, 3, 2, 4))

	expect_false(des$supports_resampling())
	expect_error(des$assign_w_to_all_subjects(), "draw_ws_raw must be implemented")

	inf = InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	expect_equal(inf$compute_estimate(), 2)
	expect_length(inf$compute_asymp_confidence_interval(), 2)
	expect_true(is.finite(inf$compute_asymp_two_sided_pval()))
	expect_error(
		inf$compute_bootstrap_two_sided_pval(B = 11),
		"Bootstrap inference is not available for plain DesignFixed objects"
	)
	expect_error(
		inf$compute_rand_two_sided_pval(r = 11),
		"Randomization inference is not available for plain DesignFixed objects"
	)

	# Test matched pair ID setter
	des$add_all_subject_matched_pair_ids(c(1, 1, 2, 2))
	expect_equal(des$get_t(), 4)
	expect_true(des$is_blocking_design())
	expect_true(des$is_complete_blocking_design())
})

test_that("DesignFixed batch ingest validates input shape and type", {
	des = EDI:::DesignFixed$new(n = 4, response_type = "continuous", verbose = FALSE)
	expect_error(des$add_all_subjects_to_experiment(matrix(1:4, ncol = 1)), "data.frame")
	expect_error(
		des$add_all_subjects_to_experiment(data.frame(x1 = 1:3)),
		"exactly 4 rows"
	)
})
