library(testthat)
library(EDI)

make_memoization_logit_inference <- function(seed = 1L, n = 80L){
	set.seed(seed)
	x = rnorm(n)
	w = rep(c(1, -1), length.out = n)
	des = EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(rbinom(n, 1L, plogis(-0.2 + 0.5 * ((w + 1) / 2) + 0.3 * x)))
	InferenceIncidLogRegr$new(des, verbose = FALSE)
}

test_that("likelihood-test memoization mixin is composed into likelihood inference classes", {
	expect_true(is.list(EDI:::InferenceMixinLikelihoodTestMemoization))
	expect_named(
		EDI:::InferenceMixinLikelihoodTestMemoization$private,
		c(
			"make_warm_fit_null_wrapper",
			"get_memoized_likelihood_test_eval",
			"get_memoized_likelihood_test_pval",
			"compute_likelihood_test_two_sided_pval"
		)
	)

	inf = make_memoization_logit_inference()
	priv = inf$.__enclos_env__$private
	expect_true(all(names(EDI:::InferenceMixinLikelihoodTestMemoization$private) %in% names(priv)))

	pval = inf$compute_score_two_sided_pval(0)
	cache_key = priv$likelihood_test_delta_key("score", 0)
	expect_true(is.finite(pval))
	expect_true(!is.null(priv$cached_values$likelihood_test_eval_cache[[cache_key]]))
})
