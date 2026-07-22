library(testthat)
library(EDI)

make_gcomp_cache_fixed_design = function(){
	set.seed(20260722)
	n = 96L
	X = data.frame(
		x1 = rep(seq(-1.5, 1.5, length.out = 12L), length.out = n),
		x2 = rep(c(-1, 0, 1), length.out = n)
	)
	des = DesignFixediBCRD$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$overwrite_all_subject_assignments(rep(c(1, -1), length.out = n))
	w = des$get_w()
	p = stats::plogis(-0.35 + 0.45 * w + 0.3 * X$x1 - 0.2 * X$x2)
	des$add_all_subject_responses(stats::rbinom(n, 1, p))
	des
}

make_gcomp_cache_kk_design = function(){
	set.seed(20260723)
	n = 96L
	X = data.frame(
		x1 = stats::rnorm(n),
		x2 = rep(c(-1, 0, 1), length.out = n)
	)
	des = DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	w = des$get_w()
	p = stats::plogis(-0.4 + 0.55 * w + 0.35 * X$x1 - 0.2 * X$x2)
	des$add_all_subject_responses(stats::rbinom(n, 1, p))
	des
}

expect_estimate_only_then_inference_ready = function(inf){
	est = as.numeric(inf$compute_estimate(estimate_only = TRUE))[1L]
	expect_true(is.finite(est))
	priv = inf$.__enclos_env__$private
	expect_true(EDI:::gcomp_standardized_effect_cache_is_ready(priv$cached_values, estimate_only = TRUE))
	expect_false(EDI:::gcomp_standardized_effect_cache_is_ready(priv$cached_values, estimate_only = FALSE))

	se = inf$get_standard_error()
	expect_true(is.finite(se))
	expect_gt(se, 0)
	expect_true(EDI:::gcomp_standardized_effect_cache_is_ready(priv$cached_values, estimate_only = FALSE))

	ci = inf$compute_asymp_confidence_interval(alpha = 0.05)
	wald_ci = inf$compute_wald_confidence_interval(alpha = 0.05)
	expect_true(all(is.finite(ci)))
	expect_equal(ci, wald_ci)
}

test_that("incidence g-computation full inference remains available after estimate-only cache", {
	des = make_gcomp_cache_fixed_design()
	expect_estimate_only_then_inference_ready(InferenceIncidGCompRiskDiff$new(des, verbose = FALSE))
	expect_estimate_only_then_inference_ready(InferenceIncidGCompRiskRatio$new(des, verbose = FALSE))
})

test_that("KK incidence g-computation full inference remains available after estimate-only cache", {
	des = make_gcomp_cache_kk_design()
	expect_estimate_only_then_inference_ready(InferenceIncidKKGCompRiskDiff$new(des, verbose = FALSE))
	expect_estimate_only_then_inference_ready(InferenceIncidKKGCompRiskRatio$new(des, verbose = FALSE))
})
