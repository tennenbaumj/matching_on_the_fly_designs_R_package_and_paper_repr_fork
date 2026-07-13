.libPaths(c("Rlib", .libPaths()))
library(testthat)
library(EDI)

set.seed(1914)
n_brt = 40
X_brt = data.frame(x1 = rnorm(n_brt), x2 = runif(n_brt))

build_brt_design = function(design_gen, effect, response_type = "continuous"){
	des = design_gen()
	if (inherits(des, "DesignSeqOneByOne")){
		for (t in 1 : n_brt){
			w_t = des$add_one_subject_to_experiment_and_assign(X_brt[t, , drop = FALSE])
			y_t = 2 * X_brt$x1[t] + effect * (w_t == 1) + rnorm(1)
			des$add_one_subject_response(t, y_t)
		}
	} else {
		des$add_all_subjects_to_experiment(X_brt)
		des$assign_w_to_all_subjects()
		w = des$get_w()
		for (t in 1 : n_brt){
			y_t = 2 * X_brt$x1[t] + effect * (w[t] == 1) + rnorm(1)
			des$add_one_subject_response(t, y_t)
		}
	}
	des
}

new_brt_inference = function(des, seed = 42L){
	inf = InferenceAllSimpleMeanDiff$new(des)
	inf$.__enclos_env__$private$seed = seed
	inf$num_cores = 1L
	inf
}

test_that("all concrete inference classes inherit the BRT layer in the right order", {
	des = build_brt_design(function() DesignFixedBernoulli$new(response_type = "continuous", n = n_brt), 0)
	inf = new_brt_inference(des)
	expect_true(is(inf, "InferenceNonParamBootstrap"))
	expect_true(is(inf, "InferenceRandBootstrap"))
	expect_true(is(inf, "InferenceRandBootstrapCI"))
	expect_true(is(inf, "InferenceBayesianBootstrap"))
	expect_true(is.function(inf$compute_rand_bootstrap_two_sided_pval))
	expect_true(is.function(inf$compute_rand_bootstrap_confidence_interval))
	expect_true(is.function(inf$approximate_rand_bootstrap_distribution_beta_hat_T))
})

test_that("BRT null distribution has B finite draws and is centered near zero under the null", {
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(response_type = "continuous", n = n_brt), 0)
	inf = new_brt_inference(des)
	t0s = inf$approximate_rand_bootstrap_distribution_beta_hat_T(B = 61, show_progress = FALSE)
	expect_length(t0s, 61)
	expect_true(all(is.finite(t0s)))
	expect_lt(abs(mean(t0s)), 1)
	# cached on repeat call
	t0s_again = inf$approximate_rand_bootstrap_distribution_beta_hat_T(B = 61, show_progress = FALSE)
	expect_identical(t0s, t0s_again)
})

test_that("BRT p-value is valid and behaves correctly under null and alternative", {
	for (design_gen in list(
		function() DesignSeqOneByOneKK14$new(response_type = "continuous", n = n_brt),
		function() DesignFixedBernoulli$new(response_type = "continuous", n = n_brt)
	)) {
		des_null = build_brt_design(design_gen, 0)
		p_null = new_brt_inference(des_null)$compute_rand_bootstrap_two_sided_pval(B = 61, show_progress = FALSE)
		expect_true(is.finite(p_null))
		expect_gte(p_null, 0); expect_lte(p_null, 1)
		expect_gt(p_null, 0.05)

		des_alt = build_brt_design(design_gen, 3)
		p_alt = new_brt_inference(des_alt)$compute_rand_bootstrap_two_sided_pval(B = 61, show_progress = FALSE)
		expect_lt(p_alt, 0.05)
	}
})

test_that("BRT p-value is deterministic given the inference seed", {
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(response_type = "continuous", n = n_brt), 1)
	p1 = new_brt_inference(des, seed = 7L)$compute_rand_bootstrap_two_sided_pval(B = 41, show_progress = FALSE)
	p2 = new_brt_inference(des, seed = 7L)$compute_rand_bootstrap_two_sided_pval(B = 41, show_progress = FALSE)
	expect_identical(p1, p2)
})

test_that("BRT p-value at delta equal to the true effect is not rejected", {
	des = build_brt_design(function() DesignFixedBernoulli$new(response_type = "continuous", n = n_brt), 2)
	inf = new_brt_inference(des)
	p_at_truth = inf$compute_rand_bootstrap_two_sided_pval(B = 61, delta = 2, show_progress = FALSE)
	expect_true(is.finite(p_at_truth))
	expect_gt(p_at_truth, 0.05)
})

test_that("BRT debug mode returns per-iteration diagnostics", {
	des = build_brt_design(function() DesignFixedBernoulli$new(response_type = "continuous", n = n_brt), 0)
	inf = new_brt_inference(des)
	dbg = inf$approximate_rand_bootstrap_distribution_beta_hat_T(B = 11, debug = TRUE, show_progress = FALSE)
	expect_named(dbg, c("values", "errors", "warnings", "num_errors", "num_warnings",
		"prop_iterations_with_errors", "prop_iterations_with_warnings", "prop_illegal_values"),
		ignore.order = TRUE)
	expect_length(dbg$values, 11)
	expect_equal(dbg$prop_illegal_values, 0)
})

test_that("BRT confidence interval covers the truth and inverts the test correctly", {
	des_alt = build_brt_design(function() DesignSeqOneByOneKK14$new(response_type = "continuous", n = n_brt), 3)
	inf_alt = new_brt_inference(des_alt)
	# B must satisfy 2/B < alpha/2 or the p-value floor prevents bracketing the bounds
	ci_alt = suppressMessages(inf_alt$compute_rand_bootstrap_confidence_interval(B = 101, pval_epsilon = 0.02, show_progress = FALSE))
	expect_length(ci_alt, 2)
	expect_true(all(is.finite(ci_alt)))
	expect_lt(ci_alt[1], ci_alt[2])
	expect_lt(ci_alt[1], 3); expect_gt(ci_alt[2], 3)  # covers truth
	expect_gt(ci_alt[1], 0)                           # excludes 0 (matches p_alt < 0.05)
	expect_named(ci_alt, c("2.5%", "97.5%"))

	des_null = build_brt_design(function() DesignFixedBernoulli$new(response_type = "continuous", n = n_brt), 0)
	inf_null = new_brt_inference(des_null)
	ci_null = suppressMessages(inf_null$compute_rand_bootstrap_confidence_interval(B = 101, pval_epsilon = 0.02, show_progress = FALSE))
	expect_true(all(is.finite(ci_null)))
	expect_lt(ci_null[1], 0); expect_gt(ci_null[2], 0)  # covers 0 under the null
})

test_that("BRT reusable-worker fast path matches the standard path exactly", {
	set.seed(310)
	des = DesignSeqOneByOneKK14$new(response_type = "count", n = n_brt)
	for (t in 1 : n_brt){
		w_t = des$add_one_subject_to_experiment_and_assign(X_brt[t, , drop = FALSE])
		lam = exp(0.3 + 0.5 * X_brt$x1[t] + 0.5 * (w_t == 1))
		des$add_one_subject_response(t, rpois(1, lam))
	}
	inf_fast = InferenceCountPoisson$new(des)
	inf_fast$.__enclos_env__$private$seed = 99L
	inf_fast$num_cores = 1L
	inf_slow = InferenceCountPoisson$new(des)
	inf_slow$.__enclos_env__$private$seed = 99L
	inf_slow$num_cores = 1L
	inf_slow$.__enclos_env__$private$reusable_bootstrap_worker_enabled = FALSE
	skip_if_not(isTRUE(inf_fast$.__enclos_env__$private$use_reusable_bootstrap_worker()))

	p_fast = inf_fast$compute_rand_bootstrap_two_sided_pval(B = 41, show_progress = FALSE)
	p_slow = inf_slow$compute_rand_bootstrap_two_sided_pval(B = 41, show_progress = FALSE)
	expect_true(is.finite(p_fast))
	expect_equal(p_fast, p_slow, tolerance = 1e-10)
})
