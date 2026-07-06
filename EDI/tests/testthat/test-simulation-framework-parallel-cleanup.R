library(testthat)
library(EDI)

test_that("SimulationFramework restores parallelism settings", {
	# Setup initial state
	set_num_cores(2L)
	initial_cores = EDI:::get_num_cores()
	initial_threads = getOption(".edi_last_set_threads")
	
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		n = 10,
		Nrep = 2,
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
	
	# Set a manual override
	assign("num_cores_override", 5L, envir = ns$edi_env)
	on.exit(assign("num_cores_override", initial_override, envir = ns$edi_env))
	
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		n = 10,
		Nrep = 1,
		num_cores = 1,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	
	expect_equal(ns$edi_env$num_cores_override, 5L)
})
