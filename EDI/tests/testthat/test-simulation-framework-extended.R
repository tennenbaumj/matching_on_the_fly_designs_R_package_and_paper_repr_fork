library(testthat)
library(EDI)

test_that("SimulationFramework supports parallel execution", {
	results_file <- tempfile(fileext = ".csv")
	set.seed(123)
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = 20L,
		Nrep = 4L,
		num_cores = 2L,
		results_filename = results_file,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	results <- SimulationFrameworkReport$new(sim)$get_results()
	expect_equal(nrow(results), 4) 
	expect_true(all(results$rep %in% 1:4))
})

test_that("SimulationFramework supports mirai-backed replication parallelism", {
	skip_if_not_installed("mirai")

	on.exit(unset_num_cores(), add = TRUE)

	set_num_cores(2L, force_mirai = TRUE)

	results_file <- tempfile(fileext = ".csv")
	set.seed(124)
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = 20L,
		Nrep = 4L,
		num_cores = 2L,
		results_filename = results_file,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)

	sim$run()
	results <- SimulationFrameworkReport$new(sim)$get_results()
	expect_equal(nrow(results), 4L)
	expect_true(all(results$rep %in% 1:4))
})

test_that("SimulationFramework supports sequential KK designs and specific inference", {
	set.seed(456)
	sim <- SimulationFramework$new(
		response_type = "count",
		design_classes_and_params = list(
			DesignSeqOneByOneKK14 = list(lambda = 0.5, t_0_pct = 0.2)
		),
		inference_classes_and_params = list(
			InferenceAllKKMeanDiffIVWC
		),
		inference_types_and_params = list(asymp_pval = list()),
		n = 10L,
		Nrep = 2L,
		betaT = 0.5,
		results_filename = tempfile(fileext = ".csv"),
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)

	sim$run()
	results <- SimulationFrameworkReport$new(sim)$get_results()
	expect_true(any(grepl("KK14", results$design)))
	expect_true(any(grepl("IVWC", results$inference)))
	expect_true(nrow(results) > 0L)
})

test_that("SimulationFramework supports randomization and bootstrap inference", {
	set.seed(789)
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(
			boot_ci = list(B = 20),
			rand_pval = list(r = 20)
		),
		n = 10L,
		Nrep = 2L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	results <- SimulationFrameworkReport$new(sim)$get_results()
	
	expect_true("boot_ci" %in% results$inference_type)
	expect_true("rand_pval" %in% results$inference_type)
	
	boot_res <- results[results$inference_type == "boot_ci", ]
	expect_true(all(is.finite(boot_res$ci_lo)))
	expect_true(all(is.finite(boot_res$ci_hi)))
})

test_that("SimulationFramework supports custom X_mat", {
	set.seed(101)
	n <- 15L
	p <- 2L
	X_custom <- matrix(rnorm(n * p), n, p)
	colnames(X_custom) <- c("feat1", "feat2")
	
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = n,
		p = p,
		X_mat = X_custom,
		cov_draw_method = NULL,
		Nrep = 2L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	results <- SimulationFrameworkReport$new(sim)$get_results()
	expect_equal(nrow(results), 2)
})

test_that("SimulationFramework can keep and retrieve intermediate data", {
	set.seed(202)
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = 10L,
		Nrep = 3L,
		keep_all_intermediate_data = TRUE,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	inter <- sim$get_all_intermediate_data()
	
	expect_length(inter, 3)
	expect_true("designs" %in% names(inter[[1]]))
	expect_true("inferences" %in% names(inter[[1]]))
	expect_true(is(inter[[1]]$designs[[1]], "Design"))
	expect_true(is(inter[[1]]$inferences[[1]][[1]], "Inference"))
	
	sim$clear_all_intermediate_data_and_gc()
	expect_null(sim$get_all_intermediate_data())
})

test_that("SimulationFramework supports custom replication and response hooks", {
	custom_data <- function(state, rep) {
		expect_equal(state$n, 6L)
		expect_true(is.null(rep) || identical(rep, 1L))
		list(
			X = data.frame(x1 = seq_len(state$n)),
			y_linear_model = seq(0, 1, length.out = state$n)
		)
	}

	custom_apply <- function(y_linear_model, w, state) {
		list(
			y = as.numeric(y_linear_model) + ifelse(as.integer(w) == 1L, 5, 0),
			dead = rep(1L, length(w))
		)
	}

	make_te <- function(beta_T) {
		function(y_linear_model, state) {
			999
		}
	}

	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list()),
		n = 6L,
		p = 1L,
		Nrep = 1L,
		betaT = 2,
		custom_replication_data_generator = custom_data,
		custom_apply_treatment_and_noise = custom_apply,
		make_estimand_fn = make_te,
		results_filename = tempfile(fileext = ".csv"),
		continue_from_last_result_row = FALSE,
		verbose = FALSE
	)

	sim$run()
	res <- SimulationFrameworkReport$new(sim)$get_results()

	expect_equal(nrow(res), 1L)
	expect_equal(res$true_estimand[[1L]], 999)
	expect_true(is.finite(res$estimate[[1L]]))
})

test_that("SimulationFramework handles Friedman nonlinear model", {
	set.seed(303)
	sim <- SimulationFramework$new(
		response_type = "continuous",
		cond_exp_func_model = "nonlinear",
		inference_types_and_params = list(asymp_pval = list()),
		p = 6L,
		n = 20L,
		Nrep = 1L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE,
		stop_on_error = FALSE
	)
	
	sim$run()
	results <- SimulationFrameworkReport$new(sim)$get_results()
	expect_equal(results$cond_exp_func_model[1], "nonlinear")
})

test_that("SimulationFramework survival responses respect censoring", {
	set.seed(404)
	# High censoring
	sim_high <- SimulationFramework$new(
		response_type = "survival",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceSurvivalCoxPHRegr),
		inference_types_and_params = list(asymp_pval = list()),
		prob_censoring = 0.99, # Force heavy censoring
		n = 50L,
		Nrep = 1L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	sim_high$run()
	
	# Low censoring
	sim_low <- SimulationFramework$new(
		response_type = "survival",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceSurvivalCoxPHRegr),
		inference_types_and_params = list(asymp_pval = list()),
		prob_censoring = 0.0,
		n = 50L,
		Nrep = 1L,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	sim_low$run()
	
	# Deep check with intermediate data
	sim_high_data <- SimulationFramework$new(
		response_type = "survival",
		design_classes_and_params = list(DesignFixedBernoulli),
		inference_classes_and_params = list(InferenceSurvivalCoxPHRegr),
		inference_types_and_params = list(asymp_pval = list()),
		prob_censoring = 1.0, 
		n = 20L,
		Nrep = 1L,
		keep_all_intermediate_data = TRUE,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	sim_high_data$run()
	inter <- sim_high_data$get_all_intermediate_data()
	des <- inter[[1]]$designs[[1]]
	
	# With prob_censoring=1.0, we expect some censored observations (dead=0)
	expect_true(any(des$.__enclos_env__$private$dead == 0))
})
