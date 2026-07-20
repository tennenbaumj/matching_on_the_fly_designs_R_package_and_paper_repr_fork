library(testthat)
library(EDI)

test_that("SimulationFramework restores parallelism settings", {
	# Setup initial state
	set_num_cores(2L, force_mirai = isTRUE(EDI:::edi_env$mirai_has_been_used))
	initial_cores = EDI:::get_num_cores()
	initial_threads = getOption(".edi_last_set_threads")
	
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		n = 10,
		Nrep_W = 2, Nrep_Y_w = 1L,
		num_cores = 1, # Run simulation serially
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	
	expect_equal(EDI:::get_num_cores(), initial_cores)
	expect_equal(getOption(".edi_last_set_threads"), initial_threads)
	
	# Cleanup
	set_num_cores(1L)
})

test_that("SimulationFramework restores num_cores_override", {
	ns = asNamespace("EDI")
	initial_override = ns$edi_env$num_cores_override
	initial_mirai_cores = ns$edi_env$global_mirai_num_cores
		
	# Set a manual override
	assign("num_cores_override", 5L, envir = ns$edi_env)
	if (isTRUE(ns$edi_env$mirai_has_been_used)) {
		assign("global_mirai_num_cores", 5L, envir = ns$edi_env)
	}
	on.exit({
		assign("num_cores_override", initial_override, envir = ns$edi_env)
		assign("global_mirai_num_cores", initial_mirai_cores, envir = ns$edi_env)
	})
	
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		n = 10,
		Nrep_W = 1, Nrep_Y_w = 1L,
		num_cores = 1,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	
	expect_equal(ns$edi_env$num_cores_override, 5L)
})

test_that("mirai use blocks later fork clusters in the same R session", {
	skip_if_not_installed("mirai")
	ns = asNamespace("EDI")
	old_mirai_used = ns$edi_env$mirai_has_been_used
	on.exit({
		unset_num_cores()
		assign("mirai_has_been_used", old_mirai_used, envir = ns$edi_env)
	}, add = TRUE)

	set_num_cores(2L, force_mirai = TRUE)
	expect_true(isTRUE(ns$edi_env$mirai_has_been_used))

	unset_num_cores()
	expect_null(ns$edi_env$global_mirai_num_cores)
	expect_true(isTRUE(ns$edi_env$mirai_has_been_used))
	expect_error(
		set_num_cores(2L),
		"Cannot switch from mirai-backed parallelism to fork-based parallelism"
	)
})
