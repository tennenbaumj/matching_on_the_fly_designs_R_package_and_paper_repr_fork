library(testthat)
library(EDI)

test_that("SimulationFramework handles multiple cells and summarizes correctly", {
	# Test grid of parameters
	results_file <- tempfile(fileext = ".csv")
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		n = c(10, 20),
		p = c(1, 2),
		betaT = c(0, 0.5),
		Nrep_W = 2, Nrep_Y_w = 1L,
		inference_types_and_params = list(asymp_pval = list(delta = 0)),
		results_filename = results_file,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	sm <- SimulationFrameworkReport$new(sim)$summarize()
	
	# Grid: 2(n) * 2(p) * 2(betaT) = 8 cells.
	expect_equal(nrow(sm), 8)
	expect_true(all(c("n", "p", "betaT", "MSE", "power") %in% names(sm)))
	
	expect_true(all(is.finite(sm$MSE)))
})

test_that("SimulationFramework continue logic works for both csv and csv.bz2", {
	for (ext in c(".csv", ".csv.bz2")) {
		results_file <- tempfile(fileext = ext)
		
		# Specify ONE inference type to get exactly 1 row per rep
		inf_types <- list(asymp_pval = list())

		# First run: 1 rep
		sim1 <- SimulationFramework$new(
			response_type = "continuous",
			design_classes_and_params = list(DesignFixedBernoulli),
			inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
			inference_types_and_params = inf_types,
			n = 10L, Nrep_W = 1L, Nrep_Y_w = 1L,
			results_filename = results_file,
			verbose = FALSE,
			continue_from_last_result_row = FALSE
		)
		sim1$run()
		expect_equal(nrow(SimulationFrameworkReport$new(sim1)$get_results()), 1)
		
		# Second run: continue, total Nrep_W = 3, Nrep_Y_w = 1L
		sim2 <- SimulationFramework$new(
			response_type = "continuous",
			design_classes_and_params = list(DesignFixedBernoulli),
			inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
			inference_types_and_params = inf_types,
			n = 10L, Nrep_W = 3L, Nrep_Y_w = 1L,
			results_filename = results_file,
			verbose = FALSE,
			continue_from_last_result_row = TRUE
		)
		sim2$run()
		res <- SimulationFrameworkReport$new(sim2)$get_results()
		expect_equal(nrow(res), 3)
		expect_equal(sort(res$rep), 1:3)
	}
})

test_that("SimulationFramework handles seed for reproducibility", {
	results_file1 <- tempfile(fileext = ".csv")
	sim1 <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = 10L, Nrep_W = 2L, Nrep_Y_w = 1L,
		seed = 12345,
		results_filename = results_file1,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	sim1$run()
	res1 <- SimulationFrameworkReport$new(sim1)$get_results()
	
	results_file2 <- tempfile(fileext = ".csv")
	sim2 <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = 10L, Nrep_W = 2L, Nrep_Y_w = 1L,
		seed = 12345,
		results_filename = results_file2,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	sim2$run()
	res2 <- SimulationFrameworkReport$new(sim2)$get_results()
	
	expect_equal(res1$estimate, res2$estimate)
})

test_that("SimulationFramework validates response_type", {
	expect_error(
		SimulationFramework$new(response_type = "invalid"),
		"must be one of"
	)
})

test_that("SimulationFramework handles factor covariate in design initialization", {
	set.seed(55)
	n <- 20
	X <- data.frame(
		x1 = rnorm(n),
		cat1 = factor(sample(letters[1:3], n, replace = TRUE))
	)
	X_mat = as.matrix(model.matrix(~ x1 + cat1, data = X)[, -1])
	
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = n,
		p = ncol(X_mat),
		X_mat = X_mat,
		cov_draw_method = NULL,
		Nrep_W = 1, Nrep_Y_w = 1L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	expect_equal(nrow(SimulationFrameworkReport$new(sim)$get_results()), 1)
})

test_that("SimulationFramework summarize handles all-NA results", {
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		Nrep_W = 1, Nrep_Y_w = 1L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	priv <- sim$.__enclos_env__$private
	# Initialize grid and combos so summarize() doesn't return early
	sim$run() 
	
	# Inject a result with NA estimate
	priv$raw_results <- list(data.table::data.table(
		response_type = "continuous",
		cond_exp_func_model = "linear",
		n = 100L, p = 5L, betaT = 1, rep = 1L,
		design = "DesignFixedBernoulli",
		inference = "InferenceAllSimpleMeanDiff",
		inference_type = "asymp_pval",
		estimate = NA_real_,
		pval = NA_real_,
		true_estimand = 1.0
	))
	priv$results_idx = 1L
	
	sm <- SimulationFrameworkReport$new(sim)$summarize()
	expect_true(is.na(sm$MSE[1]))
	expect_equal(sm$n_est[1], 0L)
})

test_that("SimulationFramework respects clamping parameters", {
	# Force an extreme signal that would normally result in 0/1 probabilities
	sim <- SimulationFramework$new(
		response_type = "incidence",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceIncidRiskDiff), # USE RD to ensure te is true_mean_diff_ate
		n = 10,
		betaT = 100, # huge effect
		incidence_clamp = 0.01,
		Nrep_W = 1, Nrep_Y_w = 1L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	res <- SimulationFrameworkReport$new(sim)$get_results()
	
	# true_estimand should be mean(clamped_p_t - clamped_p_c) <= 0.99 - 0.01 = 0.98
	expect_true(all(res$true_estimand <= 0.99))
})
