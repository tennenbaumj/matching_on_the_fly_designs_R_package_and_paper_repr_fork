library(testthat)
library(EDI)

make_ci_inversion_logit_inference <- function(seed = 1L, n = 80L){
	set.seed(seed)
	x = rnorm(n)
	w = rep(c(1, -1), length.out = n)
	des = EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(rbinom(n, 1L, plogis(-0.2 + 0.5 * ((w + 1) / 2) + 0.3 * x)))
	InferenceIncidLogRegr$new(des, verbose = FALSE)
}

test_that("CI inversion mixin is composed into likelihood inference classes", {
	expect_true(is.list(EDI:::InferenceMixinCIInversion))
	expect_named(
		EDI:::InferenceMixinCIInversion$private,
		c(
			"finalize_inverted_ci",
			"invert_test_pval_confidence_interval",
			"invert_gradient_ci_uniroot",
			"invert_lik_ratio_ci_newton"
		)
	)

	inf = make_ci_inversion_logit_inference()
	priv = inf$.__enclos_env__$private
	expect_true(all(names(EDI:::InferenceMixinCIInversion$private) %in% names(priv)))
})

test_that("CI inversion mixin accepts a valid Wald fallback and records an unusable one", {
	inf = make_ci_inversion_logit_inference(seed = 2L)
	priv = inf$.__enclos_env__$private

	ci = priv$finalize_inverted_ci(
		ci_vals = numeric(0),
		alpha = 0.1,
		est = 0,
		wald_ci = c(-1, 1),
		unavailable_reason = "ci_inversion_test_unavailable"
	)
	expect_equal(unname(ci), c(-1, 1))
	expect_named(ci, c("5%", "95%"))
	expect_false(inf$is_nonestimable("se"))

	ci = priv$finalize_inverted_ci(
		ci_vals = numeric(0),
		alpha = 0.1,
		est = 0,
		wald_ci = c(NA_real_, NA_real_),
		unavailable_reason = "ci_inversion_test_unavailable"
	)
	expect_true(all(is.na(ci)))
	expect_true(inf$is_nonestimable("se"))
	expect_equal(inf$get_nonestimable_reason(), "ci_inversion_test_unavailable")
})
