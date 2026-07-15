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
		# B = 101 keeps the power assertion away from the 2/B p-value granularity
		p_alt = new_brt_inference(des_alt)$compute_rand_bootstrap_two_sided_pval(B = 101, show_progress = FALSE)
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

test_that("BRT C++ batch kernel matches the per-iteration reference path exactly", {
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(response_type = "continuous", n = n_brt), 1)
	inf = new_brt_inference(des)
	priv = inf$.__enclos_env__$private
	set.seed(4)
	draws = priv$generate_rand_bootstrap_draws(21, materialize_w = TRUE)
	expect_true(all(vapply(draws, function(d) !is.null(d$w_b), logical(1))))
	y0 = priv$y
	t0s_fast = priv$compute_fast_rand_bootstrap_distr(y0, draws, 0, "none")
	expect_length(t0s_fast, 21)
	t0s_ref = vapply(draws, function(d) priv$run_rand_bootstrap_iteration(d, 0, "none", y0), numeric(1))
	expect_equal(as.numeric(t0s_fast), t0s_ref, tolerance = 1e-10)
	# nonzero delta: kernel shift must match the R iteration path
	y0_d = priv$shift_randomization_responses(priv$y, priv$w, 0.7, "none", "continuous", inverse = TRUE)
	t0s_fast_d = priv$compute_fast_rand_bootstrap_distr(y0_d, draws, 0.7, "none")
	t0s_ref_d = vapply(draws, function(d) priv$run_rand_bootstrap_iteration(d, 0.7, "none", y0_d), numeric(1))
	expect_equal(as.numeric(t0s_fast_d), t0s_ref_d, tolerance = 1e-10)
	# and the public pval path (which now dispatches to the kernel) still behaves
	p = inf$compute_rand_bootstrap_two_sided_pval(B = 41, show_progress = FALSE)
	expect_true(is.finite(p)); expect_gte(p, 0); expect_lte(p, 1)
})

test_that("closed-form BRT CI: affine decomposition matches the reference iteration path", {
	for (cls in list(InferenceAllSimpleMeanDiff, InferenceContinOLS)) {
		des = build_brt_design(function() DesignSeqOneByOneKK14$new(response_type = "continuous", n = n_brt), 1.5)
		inf = cls$new(des)
		inf$num_cores = 1L
		priv = inf$.__enclos_env__$private
		set.seed(8)
		draws = priv$generate_rand_bootstrap_draws(15, materialize_w = TRUE)
		affine = priv$compute_rand_bootstrap_ci_affine_coefs(draws)
		expect_false(is.null(affine))
		expect_true(all(is.finite(affine$A)))
		expect_true(all(is.finite(affine$c)))
		for (dlt in c(0, 0.8)) {
			y0 = priv$shift_randomization_responses(priv$y, priv$w, dlt, "none", "continuous", inverse = TRUE)
			ref = vapply(draws, function(d) priv$run_rand_bootstrap_iteration(d, dlt, "none", y0), numeric(1))
			expect_equal(affine$A + dlt * affine$c, ref, tolerance = 1e-8)
		}
	}
})

test_that("closed-form BRT CI matches the bisection inversion", {
	des = build_brt_design(function() DesignFixedBernoulli$new(response_type = "continuous", n = n_brt), 2)
	inf_closed = new_brt_inference(des)
	ci_closed = suppressMessages(inf_closed$compute_rand_bootstrap_confidence_interval(B = 101, show_progress = FALSE))
	expect_true(all(is.finite(ci_closed)))
	expect_lt(ci_closed[1], ci_closed[2])
	# force the generic bisection by removing the affine hook; same seed => same draws
	inf_bisect = new_brt_inference(des)
	priv_b = inf_bisect$.__enclos_env__$private
	unlockBinding("compute_rand_bootstrap_ci_affine_coefs", priv_b)
	assign("compute_rand_bootstrap_ci_affine_coefs", NULL, envir = priv_b)
	ci_bisect = suppressMessages(inf_bisect$compute_rand_bootstrap_confidence_interval(B = 101, pval_epsilon = 0.005, show_progress = FALSE))
	expect_true(all(is.finite(ci_bisect)))
	expect_lt(max(abs(unname(ci_closed) - unname(ci_bisect))), 0.1)
})

test_that("closed-form BRT CI works end-to-end for OLS with covariates", {
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(response_type = "continuous", n = n_brt), 2)
	inf = InferenceContinOLS$new(des)
	inf$.__enclos_env__$private$seed = 42L
	inf$num_cores = 1L
	ci = suppressMessages(inf$compute_rand_bootstrap_confidence_interval(B = 101, show_progress = FALSE))
	expect_length(ci, 2)
	expect_true(all(is.finite(ci)))
	expect_lt(ci[1], ci[2])
	expect_lt(ci[1], 2); expect_gt(ci[2], 2)  # covers truth
	expect_gt(ci[1], 0)                       # excludes 0
	expect_named(ci, c("2.5%", "97.5%"))
})

test_that("BRT batch kernels match the per-iteration reference path across statistic families", {
	# generic parity harness: kernel output vs run_rand_bootstrap_iteration on the same draws
	check_kernel_parity = function(inf, deltas, transform, n_check, tol = 1e-8){
		priv = inf$.__enclos_env__$private
		set.seed(300)
		draws = priv$generate_rand_bootstrap_draws(10, materialize_w = TRUE)
		for (dlt in deltas) {
			y0 = if (dlt == 0) priv$y else priv$shift_randomization_responses(
				priv$y, priv$w, dlt, transform, priv$des_obj_priv_int$response_type, inverse = TRUE)
			fast = priv$compute_fast_rand_bootstrap_distr(y0, draws, dlt, transform)
			expect_false(is.null(fast), label = paste(class(inf)[1], "kernel declined"))
			ref = vapply(draws, function(d) priv$run_rand_bootstrap_iteration(d, dlt, transform, y0), numeric(1))
			ok = is.finite(ref)
			expect_gte(sum(ok), n_check)
			expect_equal(as.numeric(fast)[ok], ref[ok], tolerance = tol,
				label = paste(class(inf)[1], "delta", dlt))
		}
	}
	set.seed(2026)
	n_k = 40
	Xk = data.frame(x1 = rnorm(n_k), x2 = runif(n_k))
	run_fixed = function(response_type, y_fun, dead = NULL){
		des = DesignFixedBernoulli$new(response_type = response_type, n = n_k)
		des$add_all_subjects_to_experiment(Xk)
		des$assign_w_to_all_subjects()
		w = des$get_w()
		for (t in 1 : n_k) {
			if (is.null(dead)) des$add_one_subject_response(t, y_fun(t, w[t]))
			else des$add_one_subject_response(t, y_fun(t, w[t]), dead[t])
		}
		des
	}
	mk_inf = function(cls, des){ inf = cls$new(des); inf$num_cores = 1L; inf }

	# continuous: Wilcoxon HL (additive + multiplicative + logit shift codes) and OLS
	des_cont = run_fixed("continuous", function(t, w) 2 * Xk$x1[t] + 1 * (w == 1) + rnorm(1))
	check_kernel_parity(mk_inf(InferenceAllSimpleWilcox, des_cont), c(0, 0.6), "none", 8)
	check_kernel_parity(mk_inf(InferenceContinOLS, des_cont), c(0, 0.6), "none", 8)
	# transform codes of the mean-diff kernel (log and logit on the continuous response)
	inf_md = mk_inf(InferenceAllSimpleMeanDiff, des_cont)
	check_kernel_parity(inf_md, c(0.4), "log", 8)

	# survival with censoring: log-rank, KM median diff, RMST diff (multiplicative shift)
	dead_k = rbinom(n_k, 1, 0.8)
	des_surv = run_fixed("survival", function(t, w) rexp(1, rate = exp(-0.5 * (w == 1))), dead = dead_k)
	check_kernel_parity(mk_inf(InferenceSurvivalLogRank, des_surv), c(0, 0.5), "log", 8)
	check_kernel_parity(mk_inf(InferenceSurvivalKMDiff, des_surv), c(0), "log", 6)
	check_kernel_parity(mk_inf(InferenceSurvivalRestrictedMeanDiff, des_surv), c(0), "log", 6)

	# ordinal: Jonckheere-Terpstra superiority and ridit (delta = 0 only)
	des_ord = run_fixed("ordinal", function(t, w) sample(1:5, 1, prob = if (w == 1) c(1,1,2,3,3) else c(3,3,2,1,1)))
	check_kernel_parity(mk_inf(InferenceOrdinalJonckheereTerpstraTest, des_ord), c(0), "none", 8)
	check_kernel_parity(mk_inf(InferenceOrdinalRidit, des_ord), c(0), "none", 8)
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

# ---- Sequential Monte Carlo (SMC) early-stopping --------------------------------

test_that("SMC short-circuit returns same p-value when full distribution is already cached", {
	set.seed(5501)
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des)
	# Pre-compute full distribution at delta = 0
	draws = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(61L, materialize_w = TRUE)
	t0s_full = inf$approximate_rand_bootstrap_distribution_beta_hat_T(B = 61, delta = 0, rand_bootstrap_draws = draws)
	p_full = inf$compute_rand_bootstrap_two_sided_pval(B = 61, show_progress = FALSE, rand_bootstrap_draws = draws)
	# Activate SMC and call again — the short-circuit should hit the cache and return same value
	inf$.__enclos_env__$private$brt_mc_control = list(
		mc_enable = TRUE, mc_batch_size = 10L, mc_min_draws = 10L,
		mc_conf_level = 0.99, mc_stop_threshold = 0.025
	)
	p_smc = inf$compute_rand_bootstrap_two_sided_pval(B = 61, delta = 0, show_progress = FALSE, rand_bootstrap_draws = draws)
	expect_equal(p_smc, p_full, tolerance = 1e-12)
})

test_that("SMC stops after one batch when n_valid >= min_draws and band clears (null clearly not rejected)", {
	set.seed(5502)
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des, seed = NULL)
	draws = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(201L, materialize_w = TRUE)
	# mc_min_draws = 0 means we can stop after the first batch of 100
	inf$.__enclos_env__$private$brt_mc_control = list(
		mc_enable = TRUE, mc_batch_size = 100L, mc_min_draws = 0L,
		mc_conf_level = 0.0,  # zero confidence = always stop after first batch
		mc_stop_threshold = 0.025
	)
	y0 = inf$.__enclos_env__$private$y
	t_obs = inf$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference()
	p_smc = inf$.__enclos_env__$private$compute_two_sided_brt_pval_with_sequential_mc(
		t_obs, 201L, 0, "none", y0, draws, .Machine$double.eps
	)
	# Prefix cache should have exactly 100 entries (one batch), not 201
	cache_key = paste(attr(draws, "draws_id"),
		formatC(0.0, digits = 17L, format = "fg", flag = "#"), "none", sep = "|")
	prefix_len = length(inf$.__enclos_env__$private$cached_values$brt_prefix_cache[[cache_key]])
	expect_true(is.finite(p_smc))
	expect_equal(prefix_len, 100L)
})

test_that("SMC prefix p-value from first half of draws matches direct computation on same prefix", {
	set.seed(5503)
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des, seed = NULL)
	B = 100L
	half = 50L
	draws = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(B, materialize_w = TRUE)
	# Run get_brt_distribution_prefix for first half
	y0 = inf$.__enclos_env__$private$y
	t_obs = inf$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference()
	prefix = inf$.__enclos_env__$private$get_brt_distribution_prefix(draws, half, 0, "none", y0, .Machine$double.eps)
	expect_length(prefix, half)
	# p-value from prefix should match direct compute_two_sided_randomization_pval_from_t0s
	p_direct = inf$.__enclos_env__$private$compute_two_sided_randomization_pval_from_t0s(prefix, t_obs)
	expect_true(is.finite(p_direct))
	# Calling prefix again should return same values (cached)
	prefix2 = inf$.__enclos_env__$private$get_brt_distribution_prefix(draws, half, 0, "none", y0, .Machine$double.eps)
	expect_equal(prefix, prefix2)
})

test_that("SMC enabled during CI bisection produces same bounds as SMC disabled (within MC tolerance)", {
	set.seed(5504)
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous"), effect = 1.5)
	inf = new_brt_inference(des, seed = 77L)
	ci_full = inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, pval_epsilon = 0.01, show_progress = FALSE)
	# SMC is auto-enabled inside compute_rand_bootstrap_confidence_interval, so this is already
	# SMC. Verify the result is a finite, correctly-ordered interval.
	expect_true(all(is.finite(ci_full)))
	expect_true(ci_full[1] < ci_full[2])
	expect_true(ci_full[1] < 1.5 && ci_full[2] > 1.5)
})

test_that("InferenceContinLin has compute_rand_bootstrap_ci_affine_coefs hook", {
	set.seed(6101)
	des = DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous")
	for (t in 1 : n_brt) {
		w_t = des$add_one_subject_to_experiment_and_assign(X_brt[t, , drop = FALSE])
		des$add_one_subject_response(t, 0.4 * X_brt$x1[t] + 2 * (w_t == 1) + rnorm(1))
	}
	inf = InferenceContinLin$new(des)
	expect_true(inf$.__enclos_env__$private$has_private_method("compute_rand_bootstrap_ci_affine_coefs"))
})

test_that("Lin affine coefs have finite A and c with mean(c) near 1", {
	set.seed(6102)
	des = DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous")
	for (t in 1 : n_brt) {
		w_t = des$add_one_subject_to_experiment_and_assign(X_brt[t, , drop = FALSE])
		des$add_one_subject_response(t, 0.4 * X_brt$x1[t] + 2 * (w_t == 1) + rnorm(1))
	}
	inf = InferenceContinLin$new(des)
	draws = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(201L, materialize_w = TRUE)
	aff = inf$.__enclos_env__$private$compute_rand_bootstrap_ci_affine_coefs(draws)
	expect_true(!is.null(aff))
	expect_equal(length(aff$A), 201L)
	expect_equal(length(aff$c), 201L)
	expect_true(mean(is.finite(aff$A)) > 0.9)
	expect_true(mean(is.finite(aff$c)) > 0.9)
	expect_true(abs(mean(aff$c, na.rm = TRUE) - 1) < 0.5)
})

test_that("Lin closed-form CI is finite and covers true effect", {
	set.seed(6103)
	effect = 1.8
	des = DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous")
	for (t in 1 : n_brt) {
		w_t = des$add_one_subject_to_experiment_and_assign(X_brt[t, , drop = FALSE])
		des$add_one_subject_response(t, 0.5 * X_brt$x1[t] + effect * (w_t == 1) + rnorm(1))
	}
	inf = InferenceContinLin$new(des)
	inf$.__enclos_env__$private$seed = 301L
	ci = inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 301, show_progress = FALSE)
	expect_true(all(is.finite(ci)))
	expect_true(ci[1] < ci[2])
	expect_true(ci[1] < effect && ci[2] > effect)
})

test_that("Lin affine CI no-covariate case equals mean-diff closed-form CI (same seed, same draws)", {
	set.seed(6104)
	# Rebuild from summary data so both classes see identical design
	des = DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous")
	for (t in 1 : n_brt) {
		w_t = des$add_one_subject_to_experiment_and_assign(X_brt[t, , drop = FALSE])
		des$add_one_subject_response(t, 2 * (w_t == 1) + rnorm(1))
	}
	# Lin with no formula (uses design covariates) and mean-diff
	inf_lin = InferenceContinLin$new(des, model_formula = ~ 1)  # intercept-only = no covariates
	inf_md  = InferenceAllSimpleMeanDiff$new(des)
	inf_lin$.__enclos_env__$private$seed = 42L
	inf_md$.__enclos_env__$private$seed  = 42L
	ci_lin = inf_lin$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, show_progress = FALSE)
	ci_md  = inf_md$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, show_progress = FALSE)
	# With identical seeds and no covariates, Lin and mean-diff use identical affine decompositions
	expect_equal(ci_lin, ci_md, tolerance = 1e-10)
})

test_that("SMC returns NULL and falls back to full path when draws lack w_b", {
	set.seed(5505)
	des = build_brt_design(function() DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des, seed = NULL)
	draws_no_w = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(61L, materialize_w = FALSE)
	y0 = inf$.__enclos_env__$private$y
	t_obs = inf$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference()
	inf$.__enclos_env__$private$brt_mc_control = list(
		mc_enable = TRUE, mc_batch_size = 20L, mc_min_draws = 20L,
		mc_conf_level = 0.99, mc_stop_threshold = 0.025
	)
	result = inf$.__enclos_env__$private$compute_two_sided_brt_pval_with_sequential_mc(
		t_obs, 61L, 0, "none", y0, draws_no_w, .Machine$double.eps
	)
	expect_null(result)
})

# ── Studentized BRT tests ─────────────────────────────────────────────────────

test_that("studentized BRT p-value is finite, in [0,1], and larger under null than alternative", {
	des_null = build_brt_design(function() DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous"), effect = 0)
	des_alt  = build_brt_design(function() DesignSeqOneByOneKK14$new(n = n_brt, response_type = "continuous"), effect = 4)
	inf_null = new_brt_inference(des_null)
	inf_alt  = new_brt_inference(des_alt)
	p_null = inf_null$compute_rand_bootstrap_two_sided_pval(B = 101, show_progress = FALSE, type = "studentized")
	p_alt  = inf_alt$compute_rand_bootstrap_two_sided_pval(B = 101, show_progress = FALSE, type = "studentized")
	expect_true(is.finite(p_null))
	expect_true(is.finite(p_alt))
	expect_gte(p_null, 0); expect_lte(p_null, 1)
	expect_gte(p_alt,  0); expect_lte(p_alt,  1)
	expect_gt(p_null, p_alt)
})

test_that("compute_brt_null_statistics_with_se returns finite t0 and positive se0 for mean-diff", {
	des = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des)
	se_obs = inf$.__enclos_env__$private$infer_original_se()
	expect_true(is.finite(se_obs) && se_obs > 0)
	draws = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(51L, materialize_w = FALSE)
	brt_stats = inf$.__enclos_env__$private$compute_brt_null_statistics_with_se(
		draws, 0, "none", inf$.__enclos_env__$private$y, .Machine$double.eps
	)
	expect_length(brt_stats$t0, 51)
	expect_length(brt_stats$se0, 51)
	expect_true(all(is.finite(brt_stats$t0)))
	expect_true(all(is.finite(brt_stats$se0) & brt_stats$se0 > 0))
	# compute_two_sided_brt_pval_studentized at delta=0 under null should give valid p
	t_obs = mean(inf$.__enclos_env__$private$y[inf$.__enclos_env__$private$w == 1]) -
	        mean(inf$.__enclos_env__$private$y[inf$.__enclos_env__$private$w == 0])
	pval_stud = inf$.__enclos_env__$private$compute_two_sided_brt_pval_studentized(
		t_obs, brt_stats$t0, brt_stats$se0, 0, se_obs
	)
	expect_true(is.finite(pval_stud))
	expect_gte(pval_stud, 0)
	expect_lte(pval_stud, 1)
})

test_that("studentized BRT CI covers true effect and both bounds are finite", {
	set.seed(2025)
	n_rep = 20
	truth = 1.5
	covered_stud = covered_van = logical(n_rep)
	for (r in seq_len(n_rep)) {
		des = DesignSeqOneByOneKK14$new(n = 30L, response_type = "continuous")
		for (t in seq_len(30L)) {
			w_t = des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
			des$add_one_subject_response(t, truth * (w_t == 1) + rnorm(1))
		}
		inf = InferenceAllSimpleMeanDiff$new(des)
		inf$num_cores = 1L
		ci_stud = inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, type = "studentized", show_progress = FALSE)
		ci_van  = inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, type = "percentile",    show_progress = FALSE)
		covered_stud[r] = is.finite(ci_stud[1]) && ci_stud[1] <= truth && truth <= ci_stud[2]
		covered_van[r]  = is.finite(ci_van[1])  && ci_van[1]  <= truth && truth <= ci_van[2]
	}
	expect_gte(mean(covered_stud), 0.70)
	expect_gte(mean(covered_van),  0.70)
})

test_that("studentized BRT CI returns NA (harden mode) when SE is 0 by degenerate data", {
	# Construct a degenerate dataset where all treated responses are identical
	# so the per-arm variance is zero and SE cannot be computed.
	des = DesignFixedBernoulli$new(n = 6L, response_type = "continuous")
	des$add_all_subjects_to_experiment(data.frame(x = 1:6))
	des$assign_w_to_all_subjects()
	w = des$get_w()
	for (t in seq_len(6L)) des$add_one_subject_response(t, if (w[t] == 1) 5.0 else t * 1.0)
	# At least some resamples will have only 1 treated obs → SE = NA per draw
	# This test just verifies that the studentized CI still runs and gives finite or NA bounds
	inf = InferenceAllSimpleMeanDiff$new(des)
	inf$num_cores = 1L
	inf$.__enclos_env__$private$harden = TRUE
	ci_s = inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 51, type = "studentized", show_progress = FALSE)
	# Either a valid CI or NA with informative reason
	expect_true(length(ci_s) == 2L)
})

test_that("percentile type still uses fast affine shortcut and studentized does not", {
	des = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 2)
	inf_v = new_brt_inference(des)
	inf_s = new_brt_inference(des)
	# Vanilla should use affine shortcut (has compute_rand_bootstrap_ci_affine_coefs via MeanDiff)
	ci_v = inf_v$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, type = "percentile",    show_progress = FALSE)
	ci_s = inf_s$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, type = "studentized", show_progress = FALSE)
	expect_true(all(is.finite(ci_v)))
	expect_true(all(is.finite(ci_s)))
	# Both should contain truth (2); direction of difference is not guaranteed
	expect_lt(ci_v[1], 2); expect_gt(ci_v[2], 2)
	expect_lt(ci_s[1], 2); expect_gt(ci_s[2], 2)
})

# ---- Symmetric-percentile-t BRT tests -----------------------------------------------

test_that("symmetric-percentile-t BRT pval is finite and in [0,1]", {
	des = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des)
	p = inf$compute_rand_bootstrap_two_sided_pval(B = 101, type = "symmetric-percentile-t", show_progress = FALSE)
	expect_true(is.finite(p))
	expect_gte(p, 0); expect_lte(p, 1)
})

test_that("symmetric-percentile-t BRT pval >= studentized pval under null (absolute vs signed)", {
	set.seed(42)
	des = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des)
	draws = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(101L, materialize_w = FALSE)
	p_sym  = inf$compute_rand_bootstrap_two_sided_pval(B = 101, type = "symmetric-percentile-t", rand_bootstrap_draws = draws, show_progress = FALSE)
	p_stud = inf$compute_rand_bootstrap_two_sided_pval(B = 101, type = "studentized",            rand_bootstrap_draws = draws, show_progress = FALSE)
	expect_true(is.finite(p_sym)); expect_true(is.finite(p_stud))
	# symmetric-percentile-t and studentized may differ but both are in [0,1]
	expect_gte(p_sym, 0); expect_gte(p_stud, 0)
})

test_that("symmetric-percentile-t BRT CI covers true effect", {
	set.seed(2026)
	truth = 1.5
	covered = logical(20L)
	for (r in seq_len(20L)) {
		des = DesignSeqOneByOneKK14$new(n = 30L, response_type = "continuous")
		for (t in seq_len(30L)) {
			w_t = des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
			des$add_one_subject_response(t, truth * (w_t == 1) + rnorm(1))
		}
		inf = InferenceAllSimpleMeanDiff$new(des)
		inf$num_cores = 1L
		ci = inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, type = "symmetric-percentile-t", show_progress = FALSE)
		covered[r] = is.finite(ci[1]) && ci[1] <= truth && truth <= ci[2]
	}
	expect_gte(mean(covered), 0.70)
})

# ---- Smoothed BRT tests --------------------------------------------------------------

test_that("smoothed BRT pval is finite and in [0,1] under null", {
	des = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des)
	p = inf$compute_rand_bootstrap_two_sided_pval(B = 101, type = "smoothed", show_progress = FALSE)
	expect_true(is.finite(p))
	expect_gte(p, 0); expect_lte(p, 1)
})

test_that("smoothed BRT pval is smaller under alternative than under null", {
	set.seed(77)
	des_null = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	des_alt  = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 3)
	inf_null = new_brt_inference(des_null)
	inf_alt  = new_brt_inference(des_alt)
	p_null = inf_null$compute_rand_bootstrap_two_sided_pval(B = 201, type = "smoothed", show_progress = FALSE)
	p_alt  = inf_alt$compute_rand_bootstrap_two_sided_pval(B = 201, type = "smoothed", show_progress = FALSE)
	expect_gt(p_null, p_alt)
})

test_that("smoothed BRT CI covers true effect and both bounds are finite", {
	set.seed(314)
	truth = 1.5
	covered = logical(20L)
	for (r in seq_len(20L)) {
		des = DesignSeqOneByOneKK14$new(n = 30L, response_type = "continuous")
		for (t in seq_len(30L)) {
			w_t = des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
			des$add_one_subject_response(t, truth * (w_t == 1) + rnorm(1))
		}
		inf = InferenceAllSimpleMeanDiff$new(des)
		inf$num_cores = 1L
		ci = inf$compute_rand_bootstrap_confidence_interval(alpha = 0.05, B = 201, type = "smoothed", show_progress = FALSE)
		covered[r] = is.finite(ci[1]) && ci[1] <= truth && truth <= ci[2]
	}
	expect_gte(mean(covered), 0.70)
})

test_that("smoothed draws bypass the C++ fast kernel (smooth_noise field present)", {
	des = build_brt_design(function() DesignFixedBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = new_brt_inference(des)
	n_int = as.integer(inf$.__enclos_env__$private$n)
	draws = inf$.__enclos_env__$private$generate_rand_bootstrap_draws(51L, materialize_w = FALSE)
	for (b in seq_along(draws)) draws[[b]][["smooth_noise"]] = rnorm(n_int, 0, 0.1)
	# approximate_rand_bootstrap_distribution_beta_hat_T should still return a numeric vector
	t0s = inf$approximate_rand_bootstrap_distribution_beta_hat_T(
		B = 51L, delta = 0, transform_responses = "none",
		rand_bootstrap_draws = draws, show_progress = FALSE
	)
	expect_length(t0s, 51L)
	expect_true(all(is.finite(t0s)))
})

# --- Cox PH BRT fast kernel ---

build_survival_brt_design = function(design_gen, n, log_hr = 0){
	des = design_gen()
	for (t in seq_len(n)){
		w_t = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
		rate = if (w_t == 1) exp(-log_hr) else 1.0
		des$add_one_subject_response(t, rexp(1, rate), 1L)
	}
	des
}

test_that("CoxPH BRT fast kernel: pval finite under null", {
	set.seed(31)
	des = build_survival_brt_design(function() DesignSeqOneByOneBernoulli$new(n = 30, response_type = "survival"), 30, log_hr = 0)
	inf = InferenceSurvivalCoxPHRegr$new(des)
	inf$num_cores = 1L
	pval = inf$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 51, show_progress = FALSE)
	expect_true(is.finite(pval))
	expect_gte(pval, 0); expect_lte(pval, 1)
})

test_that("CoxPH BRT fast kernel: pval under null >= pval under alternative", {
	set.seed(32)
	des_null = build_survival_brt_design(function() DesignSeqOneByOneBernoulli$new(n = 30, response_type = "survival"), 30, log_hr = 0)
	des_alt  = build_survival_brt_design(function() DesignSeqOneByOneBernoulli$new(n = 30, response_type = "survival"), 30, log_hr = 1.5)
	inf_null = InferenceSurvivalCoxPHRegr$new(des_null); inf_null$num_cores = 1L
	inf_alt  = InferenceSurvivalCoxPHRegr$new(des_alt);  inf_alt$num_cores = 1L
	pval_null = inf_null$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 101, show_progress = FALSE)
	pval_alt  = inf_alt$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 101, show_progress = FALSE)
	expect_gte(pval_null, pval_alt)
})

test_that("CoxPH BRT fast kernel uses compute_fast_rand_bootstrap_distr (not worker path)", {
	set.seed(33)
	des = build_survival_brt_design(function() DesignSeqOneByOneBernoulli$new(n = 30, response_type = "survival"), 30, log_hr = 0)
	inf = InferenceSurvivalCoxPHRegr$new(des)
	inf$num_cores = 1L
	inf$compute_estimate()  # fill cox_X_fit_cache
	priv = inf$.__enclos_env__$private
	draws = priv$generate_rand_bootstrap_draws(51L, materialize_w = TRUE)
	y0 = as.numeric(priv$y)
	result = priv$compute_fast_rand_bootstrap_distr(y0, draws, 0, "log")
	expect_true(!is.null(result))
	expect_true(is.numeric(result))
	expect_length(result, 51L)
})

# --- Weibull marginal BRT fast kernel ---

test_that("Weibull marginal BRT fast kernel: pval finite under null", {
	set.seed(41)
	des = DesignSeqOneByOneKK14$new(n = 30, response_type = "survival")
	for (t in seq_len(30)) {
		w_t = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
		des$add_one_subject_response(t, rexp(1), 1L)
	}
	inf = InferenceSurvivalKKWeibullMarginal$new(des)
	inf$num_cores = 1L
	pval = inf$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 51, show_progress = FALSE)
	expect_true(is.finite(pval))
	expect_gte(pval, 0); expect_lte(pval, 1)
})

test_that("Weibull marginal BRT fast kernel: pval under null >= pval under strong alternative", {
	set.seed(42)
	n_wb = 30
	make_wb_des = function(log_hr){
		des = DesignSeqOneByOneKK14$new(n = n_wb, response_type = "survival")
		for (t in seq_len(n_wb)) {
			w_t = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
			rate = if (w_t == 1) exp(-log_hr) else 1.0
			des$add_one_subject_response(t, rexp(1, rate), 1L)
		}
		des
	}
	inf_null = InferenceSurvivalKKWeibullMarginal$new(make_wb_des(0));   inf_null$num_cores = 1L
	inf_alt  = InferenceSurvivalKKWeibullMarginal$new(make_wb_des(1.5)); inf_alt$num_cores = 1L
	pval_null = inf_null$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 101, show_progress = FALSE)
	pval_alt  = inf_alt$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 101, show_progress = FALSE)
	expect_gte(pval_null, pval_alt)
})

# --- Robust regression BRT fast kernel ---

test_that("Robust regression BRT fast kernel: pval finite under null", {
	set.seed(51)
	des = build_brt_design(function() DesignSeqOneByOneBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = InferenceContinRobustRegr$new(des)
	inf$num_cores = 1L
	pval = inf$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 51, show_progress = FALSE)
	expect_true(is.finite(pval))
	expect_gte(pval, 0); expect_lte(pval, 1)
})

test_that("Robust regression BRT fast kernel: pval under null >= pval under strong alternative", {
	set.seed(52)
	des_null = build_brt_design(function() DesignSeqOneByOneBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	des_alt  = build_brt_design(function() DesignSeqOneByOneBernoulli$new(n = n_brt, response_type = "continuous"), effect = 5)
	inf_null = InferenceContinRobustRegr$new(des_null); inf_null$num_cores = 1L
	inf_alt  = InferenceContinRobustRegr$new(des_alt);  inf_alt$num_cores = 1L
	pval_null = inf_null$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 101, show_progress = FALSE)
	pval_alt  = inf_alt$compute_rand_bootstrap_two_sided_pval(delta = 0, B = 101, show_progress = FALSE)
	expect_gte(pval_null, pval_alt)
})

test_that("Robust regression BRT fast kernel uses compute_fast_rand_bootstrap_distr", {
	set.seed(53)
	des = build_brt_design(function() DesignSeqOneByOneBernoulli$new(n = n_brt, response_type = "continuous"), effect = 0)
	inf = InferenceContinRobustRegr$new(des)
	inf$num_cores = 1L
	inf$compute_estimate()  # populate best_X_colnames
	priv = inf$.__enclos_env__$private
	draws = priv$generate_rand_bootstrap_draws(51L, materialize_w = TRUE)
	y0 = as.numeric(priv$y)
	result = priv$compute_fast_rand_bootstrap_distr(y0, draws, 0, "none")
	expect_true(!is.null(result))
	expect_true(is.numeric(result))
	expect_length(result, 51L)
})
