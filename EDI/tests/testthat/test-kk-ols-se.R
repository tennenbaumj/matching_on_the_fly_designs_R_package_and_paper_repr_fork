library(testthat)
library(EDI)

make_kk_ols_binary_match_fixture = function(){
	des = DesignFixedBinaryMatch$new(
		response_type = "continuous",
		n = 6,
		m = c(1L, 1L, 2L, 2L, 3L, 3L),
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(x = c(-2, -1, 0, 1, 2, 3)))
	des$assign_w_to_all_subjects(c(1, -1, -1, 1, 1, -1))
	des$add_all_subject_responses(c(-1, 0, 1, 1, 4, 9))
	des
}

test_that("KK OLS one-lik uses HC2 post-fit standard errors and finite df", {
	des = make_kk_ols_binary_match_fixture()

	inf = InferenceContinKKOLSOneLik$new(des, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	X = rbind(
		cbind(0, 1, c(-2, 0, 2)),
		cbind(1, c(0, 1, 0), c(-1, 1, 3))
	)
	y = c(1, 1, 5, 0, 1, 9)
	j_treat = 2L

	fit = priv$fit_ols(X, y, j_treat = j_treat, estimate_only = FALSE)
	post_fit = EDI:::ols_hc2_post_fit_cpp(X, y, fit$b, j_treat)

	expect_equal(fit$se, as.numeric(post_fit$std_err[j_treat]), tolerance = 1e-10)
	expect_equal(fit$vcov, post_fit$vcov, tolerance = 1e-10)
	expect_equal(fit$df, nrow(X) - ncol(X))

	inf$compute_estimate()
	expect_true(is.finite(priv$cached_values$df))
	expect_gt(priv$cached_values$df, 0)
})

test_that("KK OLS IVWC uses HC2 component standard errors and finite df", {
	des = make_kk_ols_binary_match_fixture()

	inf = InferenceContinKKOLSIVWC$new(des, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	X = cbind(1, c(0, 1, 0, 1, 0, 1), c(-2, -1, 0, 1, 2, 3))
	y = c(-1, 0, 1, 1, 4, 9)
	j_treat = 2L

	fit = priv$fit_ols_with_treatment(X, y, j_treat = j_treat, estimate_only = FALSE)
	post_fit = EDI:::ols_hc2_post_fit_cpp(X, y, as.numeric(stats::coef(stats::lm.fit(X, y))), j_treat)

	expect_equal(sqrt(fit$ssq), as.numeric(post_fit$std_err[j_treat]), tolerance = 1e-10)
	expect_equal(fit$df, nrow(X) - ncol(X))
})

test_that("KK OLS one-lik supports gradient p-values with warm null-fit cache reuse", {
	des = make_kk_ols_binary_match_fixture()

	inf = InferenceContinKKOLSOneLik$new(des, verbose = FALSE)
	expect_true("gradient" %in% inf$get_supported_testing_types())
	inf$set_testing_type("gradient")

	p1 = inf$compute_asymp_two_sided_pval(0.2)
	p2 = inf$compute_asymp_two_sided_pval(0.25)
	cache = inf$.__enclos_env__$private$likelihood_null_warm_cache[["likelihood_test:gradient"]]
	ci = inf$compute_asymp_confidence_interval(alpha = 0.2)
	ci_cache = inf$.__enclos_env__$private$likelihood_null_warm_cache[["gradient_ci"]]

	expect_true(is.finite(p1))
	expect_true(is.finite(p2))
	expect_true(all(is.finite(ci)))
	expect_gte(p1, 0)
	expect_lte(p1, 1)
	expect_gte(p2, 0)
	expect_lte(p2, 1)
	expect_true(!is.null(cache))
	expect_true(!is.null(cache$start))
	expect_true(length(cache$start) > 0L)
	expect_true(!is.null(ci_cache))
	expect_true(!is.null(ci_cache$start))
	expect_true(length(ci_cache$start) > 0L)
})
