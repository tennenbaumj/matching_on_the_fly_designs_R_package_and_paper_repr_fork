library(testthat)
library(EDI)

make_logit_inference <- function(seed = 1, n = 120){
	set.seed(seed)
	X <- data.frame(x1 = rnorm(n))
	y <- rbinom(n, 1, plogis(0.4 + 0.6 * X$x1))
	des <- DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	InferenceIncidLogRegr$new(des, verbose = FALSE)
}

test_that("InferenceIncidLogRegr auto-inherits approx Bartlett support (via InferenceParamBootstrap) and advertises the testing type", {
	inf <- make_logit_inference(seed = 201)
	expect_true("lik_ratio_bartlett_approx" %in% inf$get_supported_testing_types())
})

test_that("B defaults to 99 on the public wrapper methods", {
	inf <- make_logit_inference(seed = 208)
	expect_equal(formals(inf$compute_lik_ratio_bartlett_approx_two_sided_pval)$B, 99)
	expect_equal(formals(inf$compute_lik_ratio_bartlett_approx_confidence_interval)$B, 99)
})

test_that("Bartlett factor: is finite, positive, and reproducible given a fixed seed and B", {
	inf <- make_logit_inference(seed = 202)
	priv <- inf$.__enclos_env__$private
	spec <- priv$get_likelihood_test_spec()

	inf$set_seed(777)
	factor1 <- priv$get_bartlett_factor_approx(spec = spec, delta = 0, full_fit = spec$full_fit, null_fit = spec$fit_null(0), B = 99)
	factor2 <- priv$get_bartlett_factor_approx(spec = spec, delta = 0, full_fit = spec$full_fit, null_fit = spec$fit_null(0), B = 99)

	expect_true(is.finite(factor1))
	expect_gt(factor1, 0)
	expect_equal(factor1, factor2, tolerance = 0)

	# Sanity: with a genuine null-plausible delta and enough replicates, the MC
	# estimate of E[LR] shouldn't be wildly far from its chi-sq(1) target of 1.
	expect_true(factor1 > 0.1 && factor1 < 10)
})

test_that("Bartlett factor changes with a different seed (not a hard-coded constant)", {
	inf <- make_logit_inference(seed = 203)
	priv <- inf$.__enclos_env__$private
	spec <- priv$get_likelihood_test_spec()
	null_fit <- spec$fit_null(0)

	inf$set_seed(111)
	factor_a <- priv$get_bartlett_factor_approx(spec = spec, delta = 0, full_fit = spec$full_fit, null_fit = null_fit, B = 99)
	inf$set_seed(222)
	factor_b <- priv$get_bartlett_factor_approx(spec = spec, delta = 0, full_fit = spec$full_fit, null_fit = null_fit, B = 99)

	expect_true(is.finite(factor_a) && is.finite(factor_b))
	expect_false(isTRUE(all.equal(factor_a, factor_b)))
})

test_that("Bartlett factor changes with a different B (same seed), since more/fewer replicates are averaged", {
	inf <- make_logit_inference(seed = 209)
	priv <- inf$.__enclos_env__$private
	spec <- priv$get_likelihood_test_spec()
	null_fit <- spec$fit_null(0)

	inf$set_seed(321)
	factor_small_B <- priv$get_bartlett_factor_approx(spec = spec, delta = 0, full_fit = spec$full_fit, null_fit = null_fit, B = 25)
	inf$set_seed(321)
	factor_large_B <- priv$get_bartlett_factor_approx(spec = spec, delta = 0, full_fit = spec$full_fit, null_fit = null_fit, B = 150)

	expect_true(is.finite(factor_small_B) && is.finite(factor_large_B))
	expect_false(isTRUE(all.equal(factor_small_B, factor_large_B)))
})

test_that("Bartlett p-value matches pchisq(raw LR statistic / factor, df=1) exactly, using an explicit B", {
	inf <- make_logit_inference(seed = 204)
	priv <- inf$.__enclos_env__$private
	spec <- priv$get_likelihood_test_spec()

	inf$set_seed(999)
	raw_pval <- inf$compute_lik_ratio_two_sided_pval(delta = 0)
	expect_true(is.finite(raw_pval))
	raw_stat <- qchisq(raw_pval, df = 1, lower.tail = FALSE)

	null_fit_for_factor <- spec$fit_null(0)
	inf$set_seed(999)
	factor <- priv$get_bartlett_factor_approx(spec = spec, delta = 0, full_fit = spec$full_fit, null_fit = null_fit_for_factor, B = 60)
	expected_bartlett_pval <- pchisq(raw_stat / factor, df = 1, lower.tail = FALSE)

	inf$set_seed(999)
	bartlett_pval <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 60)

	expect_true(is.finite(bartlett_pval))
	expect_equal(bartlett_pval, expected_bartlett_pval, tolerance = 1e-4)
})

test_that("Bartlett p-value is reproducible for the same (delta, B, seed) and differs when B changes", {
	inf <- make_logit_inference(seed = 210)
	inf$set_seed(654)
	pval_B60_a <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 60)
	inf$set_seed(654)
	pval_B60_b <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 60)
	inf$set_seed(654)
	pval_B120 <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 120)

	expect_true(is.finite(pval_B60_a) && is.finite(pval_B60_b) && is.finite(pval_B120))
	expect_equal(pval_B60_a, pval_B60_b, tolerance = 0)
	expect_false(isTRUE(all.equal(pval_B60_a, pval_B120)))
})

test_that("Bartlett p-value differs from the raw LR p-value whenever the factor is not 1", {
	inf <- make_logit_inference(seed = 205)
	inf$set_seed(555)
	raw_pval <- inf$compute_lik_ratio_two_sided_pval(delta = 0)
	inf$set_seed(555)
	bartlett_pval <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0)

	expect_true(is.finite(raw_pval))
	expect_true(is.finite(bartlett_pval))
	expect_true(bartlett_pval >= 0 && bartlett_pval <= 1)
})

test_that("Bartlett confidence interval is finite, ordered, and brackets the estimate", {
	inf <- make_logit_inference(seed = 206)
	inf$set_seed(333)
	est <- inf$compute_estimate()
	ci <- inf$compute_lik_ratio_bartlett_approx_confidence_interval(alpha = 0.2)

	expect_length(ci, 2)
	expect_true(all(is.finite(ci)))
	expect_true(ci[[1]] <= est && est <= ci[[2]])
})

test_that("Bartlett confidence interval honors an explicit B and stays finite/ordered", {
	inf <- make_logit_inference(seed = 211)
	inf$set_seed(1357)
	est <- inf$compute_estimate()
	ci <- inf$compute_lik_ratio_bartlett_approx_confidence_interval(alpha = 0.2, B = 40)

	expect_length(ci, 2)
	expect_true(all(is.finite(ci)))
	expect_true(ci[[1]] <= est && est <= ci[[2]])
})

test_that("Bartlett confidence interval via configured testing_type dispatch matches the direct wrapper", {
	inf <- make_logit_inference(seed = 207)
	inf$set_testing_type("lik_ratio_bartlett_approx")
	inf$set_seed(444)
	ci_dispatch <- inf$compute_asymp_confidence_interval(alpha = 0.2)
	inf$set_seed(444)
	ci_direct <- inf$compute_lik_ratio_bartlett_approx_confidence_interval(alpha = 0.2)

	expect_equal(ci_dispatch, ci_direct, tolerance = 1e-6)
})

test_that("The seed is inherited from the inference object's own set_seed(), with no separate seed argument", {
	inf <- make_logit_inference(seed = 212)
	expect_false("seed" %in% names(formals(inf$compute_lik_ratio_bartlett_approx_two_sided_pval)))
	expect_false("seed" %in% names(formals(inf$compute_lik_ratio_bartlett_approx_confidence_interval)))

	inf$set_seed(2468)
	pval_a <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 30)
	inf$set_seed(2468)
	pval_b <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 30)
	expect_equal(pval_a, pval_b, tolerance = 0)
})
