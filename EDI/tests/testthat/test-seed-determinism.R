library(testthat)
library(EDI)

### Helpers ###################################################################

make_fixed_design = function(seed = NULL, n = 10) {
	des = DesignFixedBernoulli$new(n = n, response_type = "continuous", seed = seed)
	set.seed(1)
	des$add_all_subjects_to_experiment(data.frame(x = rnorm(n)))
	des
}

make_completed_fixed_design = function(seed = NULL, n = 10) {
	des = make_fixed_design(seed = seed, n = n)
	des$assign_w_to_all_subjects()
	set.seed(2)
	des$add_all_subject_responses(rnorm(n))
	des
}

make_sequential_design = function(seed = NULL, n = 8) {
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", seed = seed)
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(data.frame(x = i))
	des
}

make_completed_sequential_design = function(seed = NULL, n = 8) {
	des = make_sequential_design(seed = seed, n = n)
	set.seed(4)
	des$add_all_subject_responses(rnorm(n))
	des
}

### Design: fixed ############################################################

test_that("DesignFixed seed: same seed produces same draw_ws_according_to_design", {
	des = make_fixed_design(seed = 42)
	w1 = des$draw_ws_according_to_design(r = 10)
	w2 = des$draw_ws_according_to_design(r = 10)
	expect_identical(w1, w2)
})

test_that("DesignFixed seed: different seeds produce different draws", {
	des1 = make_fixed_design(seed = 42)
	des2 = make_fixed_design(seed = 99)
	w1 = des1$draw_ws_according_to_design(r = 20)
	w2 = des2$draw_ws_according_to_design(r = 20)
	expect_false(identical(w1, w2))
})

test_that("DesignFixed seed: NULL seed is non-deterministic across calls", {
	des = make_fixed_design(seed = NULL)
	w1 = des$draw_ws_according_to_design(r = 20)
	w2 = des$draw_ws_according_to_design(r = 20)
	# With no seed, successive calls should (almost certainly) differ
	expect_false(identical(w1, w2))
})

test_that("DesignFixed seed: duplicate() clears seed from clone", {
	des = make_completed_fixed_design(seed = 42)
	clone = des$duplicate()
	expect_null(clone$.__enclos_env__$private$seed)
	# Original retains seed
	expect_equal(des$.__enclos_env__$private$seed, 42L)
})

test_that("DesignFixed seed: two objects with same seed give identical draws", {
	des1 = make_fixed_design(seed = 7)
	des2 = make_fixed_design(seed = 7)
	expect_identical(
		des1$draw_ws_according_to_design(r = 15),
		des2$draw_ws_according_to_design(r = 15)
	)
})

### Design: sequential #######################################################

test_that("DesignSeqOneByOne seed: same seed produces same assignment sequence", {
	des1 = make_sequential_design(seed = 11)
	des2 = make_sequential_design(seed = 11)
	expect_identical(des1$get_w(), des2$get_w())
})

test_that("DesignSeqOneByOne seed: different seeds produce different assignments", {
	# Use n=30 so assignment vectors are long enough that collision is negligible
	des1 = make_sequential_design(seed = 11, n = 30)
	des2 = make_sequential_design(seed = 22, n = 30)
	expect_false(identical(des1$get_w(), des2$get_w()))
})

test_that("DesignSeqOneByOne seed: duplicate() clears seed from clone", {
	des = make_completed_sequential_design(seed = 55)
	clone = des$duplicate()
	expect_null(clone$.__enclos_env__$private$seed)
	expect_equal(des$.__enclos_env__$private$seed, 55L)
})

### Inference: serial ########################################################

test_that("Inference seed: same seed gives same rand p-value (serial)", {
	des = make_completed_fixed_design(seed = NULL)
	inf1 = InferenceAllSimpleMeanDiff$new(des); inf1$set_seed(42)
	inf2 = InferenceAllSimpleMeanDiff$new(des); inf2$set_seed(42)
	p1 = inf1$compute_rand_two_sided_pval(r = 49, show_progress = FALSE)
	p2 = inf2$compute_rand_two_sided_pval(r = 49, show_progress = FALSE)
	expect_identical(p1, p2)
})

test_that("Inference seed: different seeds give different bootstrap distributions (serial)", {
	# Use bootstrap distribution (continuous) instead of rand p-value (discrete count)
	# to avoid coincidental equality; exact numerical match across seeds is essentially impossible
	des = make_completed_fixed_design(seed = NULL, n = 20)
	inf1 = InferenceAllSimpleMeanDiff$new(des); inf1$set_seed(100)
	inf2 = InferenceAllSimpleMeanDiff$new(des); inf2$set_seed(200)
	d1 = inf1$approximate_bootstrap_distribution_beta_hat_T(B = 99, show_progress = FALSE)
	d2 = inf2$approximate_bootstrap_distribution_beta_hat_T(B = 99, show_progress = FALSE)
	expect_false(identical(d1, d2))
})

test_that("Inference seed: same seed gives same bootstrap CI (serial)", {
	des = make_completed_fixed_design(seed = NULL, n = 20)
	inf1 = InferenceAllSimpleMeanDiff$new(des); inf1$set_seed(13)
	inf2 = InferenceAllSimpleMeanDiff$new(des); inf2$set_seed(13)
	ci1 = inf1$compute_bootstrap_confidence_interval(B = 49, show_progress = FALSE)
	ci2 = inf2$compute_bootstrap_confidence_interval(B = 49, show_progress = FALSE)
	expect_identical(ci1, ci2)
})

### Inference: multi-core (fork cluster, Unix only) ##########################

test_that("Inference seed: same seed + same num_cores gives same rand p-value (fork cluster)", {
	skip_on_os("windows")
	skip_if(
		isTRUE(EDI:::edi_env$mirai_has_been_used),
		"fork clusters cannot be started safely after mirai has been used in this R session"
	)
	des = make_completed_fixed_design(seed = NULL, n = 20)
	set_num_cores(2)
	on.exit(unset_num_cores(), add = TRUE)
	inf1 = InferenceAllSimpleMeanDiff$new(des); inf1$set_seed(42)
	inf2 = InferenceAllSimpleMeanDiff$new(des); inf2$set_seed(42)
	p1 = inf1$compute_rand_two_sided_pval(r = 99, show_progress = FALSE)
	p2 = inf2$compute_rand_two_sided_pval(r = 99, show_progress = FALSE)
	expect_identical(p1, p2)
})

### SimulationFramework: seed ################################################

test_that("SimulationFramework seed: same seed gives same estimates (serial)", {
	run_sim = function() {
		f = tempfile(fileext = ".csv")
		sim = SimulationFramework$new(
			response_type = "continuous",
			design_classes_and_params = list(DesignFixedBernoulli),
			inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
				inference_types_and_params = list(asymp_pval = list()),
				n = 10L, Nrep_W = 1L, Nrep_Y_w = 3L, seed = 321,
				num_cores = 1L,
				results_filename = f, verbose = FALSE,
				continue_from_last_result_row = FALSE
		)
		sim$run()
		SimulationFrameworkReport$new(sim)$get_results()
	}
	res1 = run_sim()
	res2 = run_sim()
	expect_equal(res1$estimate, res2$estimate)
})

test_that("SimulationFramework seed: different seeds give different estimates", {
	run_sim = function(seed) {
		f = tempfile(fileext = ".csv")
		sim = SimulationFramework$new(
			response_type = "continuous",
			design_classes_and_params = list(DesignFixedBernoulli),
			inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
				inference_types_and_params = list(asymp_pval = list()),
				n = 10L, Nrep_W = 1L, Nrep_Y_w = 5L, seed = seed,
				num_cores = 1L,
				results_filename = f, verbose = FALSE,
				continue_from_last_result_row = FALSE
		)
		sim$run()
		SimulationFrameworkReport$new(sim)$get_results()
	}
	res1 = run_sim(100)
	res2 = run_sim(200)
	expect_false(identical(res1$estimate, res2$estimate))
})
