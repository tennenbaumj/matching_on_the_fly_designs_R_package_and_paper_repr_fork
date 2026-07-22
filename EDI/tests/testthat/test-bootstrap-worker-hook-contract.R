make_bootstrap_worker_contract_inference <- function(){
	set.seed(20260722)
	n = 48L
	X = data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	des$add_all_subject_responses(stats::rpois(n, exp(0.2 + 0.5 * w + 0.15 * X$x1)))
	InferenceCountPoisson$new(des, verbose = FALSE)
}

test_that("reusable bootstrap worker hooks create, load, and estimate a draw", {
	inf = make_bootstrap_worker_contract_inference()
	priv = inf$.__enclos_env__$private

	expect_true(priv$use_reusable_bootstrap_worker())
	worker_state = priv$create_reusable_bootstrap_worker()
	expect_true(is.list(worker_state))
	expect_false(is.null(worker_state$worker))

	draw = priv$bootstrap_sample_indices(priv$n)
	expect_invisible(priv$load_bootstrap_draw_into_worker(worker_state, draw))
	theta = priv$estimate_bootstrap_worker(worker_state)
	expect_type(theta, "double")
	expect_length(theta, 1L)
	expect_true(is.finite(theta))
})

test_that("generic reusable worker executor serves bootstrap and jackknife wrappers", {
	inf = make_bootstrap_worker_contract_inference()
	priv = inf$.__enclos_env__$private
	draws = replicate(3L, priv$bootstrap_sample_indices(priv$n), simplify = FALSE)

	executor_values = priv$compute_reusable_bootstrap_worker_distribution(
		draws = draws,
		actual_cores = 1L
	)
	bootstrap_values = priv$compute_bootstrap_distribution_with_reused_workers(
		boot_draws = draws,
		actual_cores = 1L
	)
	jackknife_values = priv$compute_jackknife_distribution_with_reused_workers(
		deletion_draws = draws,
		actual_cores = 1L
	)

	expect_type(executor_values, "double")
	expect_length(executor_values, length(draws))
	expect_equal(bootstrap_values, executor_values)
	expect_equal(jackknife_values, executor_values)
	expect_identical(
		priv$compute_reusable_bootstrap_worker_distribution(draws = list(), actual_cores = 1L),
		numeric(0)
	)

	FailingWorkerPoisson = R6::R6Class(
		"FailingWorkerPoisson",
		inherit = InferenceCountPoisson,
		private = list(
			load_bootstrap_sample_into_worker = function(worker_state, indices) invisible(NULL),
			compute_bootstrap_worker_estimate = function(worker_state) stop("intentional worker failure")
		)
	)
	failing_inf = FailingWorkerPoisson$new(inf$get_design_object(), verbose = FALSE)
	failing_priv = failing_inf$.__enclos_env__$private
	expect_equal(
		failing_priv$compute_reusable_bootstrap_worker_distribution(draws = draws, actual_cores = 1L),
		rep(NA_real_, length(draws))
	)
})

test_that("generic reusable worker executor supports non-default resampling contracts", {
	ContractOperationPoisson = R6::R6Class(
		"ContractOperationPoisson",
		inherit = InferenceCountPoisson,
		private = list(
			load_rand_bootstrap_draw_into_worker = function(worker_state, draw, multiplier){
				worker_state$worker$.__enclos_env__$private$cached_values$contract_test_value =
					as.numeric(draw$value) * multiplier
				invisible(worker_state)
			},
			compute_bootstrap_worker_estimate = function(worker_state){
				worker_state$worker$.__enclos_env__$private$cached_values$contract_test_value
			}
		)
	)

	base_inf = make_bootstrap_worker_contract_inference()
	inf = ContractOperationPoisson$new(base_inf$get_design_object(), verbose = FALSE)
	priv = inf$.__enclos_env__$private
	values = priv$compute_reusable_bootstrap_worker_distribution(
		draws = list(list(value = 2), list(value = 4)),
		actual_cores = 1L,
		operation = "rand_bootstrap",
		loader_args = list(multiplier = 10)
	)

	expect_equal(values, c(20, 40))
})

test_that("reusable bootstrap worker contract rejects invalid worker state and estimate output", {
	MissingWorkerStatePoisson = R6::R6Class(
		"MissingWorkerStatePoisson",
		inherit = InferenceCountPoisson,
		private = list(create_bootstrap_worker_state = function() list())
	)
	MissingEstimatePoisson = R6::R6Class(
		"MissingEstimatePoisson",
		inherit = InferenceCountPoisson,
		private = list(compute_bootstrap_worker_estimate = function(worker_state) numeric(0))
	)

	base_inf = make_bootstrap_worker_contract_inference()
	missing_state_inf = MissingWorkerStatePoisson$new(base_inf$get_design_object(), verbose = FALSE)
	missing_state_priv = missing_state_inf$.__enclos_env__$private
	expect_error(
		missing_state_priv$create_reusable_bootstrap_worker(),
		"non-NULL lists containing `worker`"
	)

	missing_estimate_inf = MissingEstimatePoisson$new(base_inf$get_design_object(), verbose = FALSE)
	missing_estimate_priv = missing_estimate_inf$.__enclos_env__$private
	worker_state = missing_estimate_priv$create_reusable_bootstrap_worker()
	expect_error(
		missing_estimate_priv$estimate_bootstrap_worker(worker_state),
		"one numeric treatment estimate"
	)
})
