library(testthat)
library(EDI)

make_public_contract_logit_design <- function(seed = 1L, n = 100L){
	set.seed(seed)
	x = rnorm(n)
	w = rep(c(1, -1), length.out = n)
	des = EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(rbinom(n, 1L, plogis(-0.2 + 0.5 * ((w + 1) / 2) + 0.3 * x)))
	des
}

make_public_contract_poisson_inference <- function(seed = 1L, n = 100L){
	set.seed(seed)
	x = rnorm(n)
	w = rep(c(1, -1), length.out = n)
	des = EDI:::DesignFixed$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(rpois(n, exp(0.2 + 0.4 * ((w + 1) / 2) + 0.2 * x)))
	InferenceCountPoisson$new(des, verbose = FALSE)
}

make_constant_score_pval_logit_inference <- function(des){
	ext_env = new.env(parent = globalenv())
	ext_env$R6Class = R6::R6Class
	ext_env$InferenceIncidLogRegr = InferenceIncidLogRegr
	evalq({
		ConstantScorePvalLogit = R6Class(
			"ConstantScorePvalLogit",
			inherit = InferenceIncidLogRegr,
			private = list(
				get_likelihood_test_spec = function() list(j = 2L),
				get_memoized_likelihood_test_pval = function(...) 1
			)
		)
	}, envir = ext_env)
	ext_env$ConstantScorePvalLogit$new(des, verbose = FALSE)
}

make_no_likelihood_spec_logit_inference <- function(des){
	ext_env = new.env(parent = globalenv())
	ext_env$R6Class = R6::R6Class
	ext_env$InferenceIncidLogRegr = InferenceIncidLogRegr
	evalq({
		NoLikelihoodSpecLogit = R6Class(
			"NoLikelihoodSpecLogit",
			inherit = InferenceIncidLogRegr,
			private = list(
				get_likelihood_test_spec = function() NULL,
				compute_likelihood_test_two_sided_pval = function(delta, testing_type, bartlett_B = NULL){
					eval(body(EDI:::InferenceMixinLikelihoodTestMemoization$private$compute_likelihood_test_two_sided_pval))
				}
			)
		)
	}, envir = ext_env)
	ext_env$NoLikelihoodSpecLogit$new(des, verbose = FALSE)
}

test_that("public score CI falls back to the Wald interval when inversion cannot bracket", {
	inf = make_constant_score_pval_logit_inference(make_public_contract_logit_design(seed = 11L))

	ci = inf$compute_score_confidence_interval(alpha = 0.1)
	wald_ci = inf$compute_wald_confidence_interval(alpha = 0.1)

	expect_true(all(is.finite(ci)))
	expect_equal(as.numeric(ci), as.numeric(wald_ci), tolerance = 1e-12)
	expect_false(inf$is_nonestimable("se"))
})

test_that("public information preferences select available score-test information", {
	inf = make_public_contract_poisson_inference(seed = 12L)

	inf$set_information_preference("observed")
	p_observed = inf$compute_score_two_sided_pval(delta = 0)
	expect_true(is.finite(p_observed))
	expect_equal(inf$get_information_source_used(), "observed")

	inf$set_information_preference("auto")
	p_auto = inf$compute_score_two_sided_pval(delta = 0.1)
	expect_true(is.finite(p_auto))
	expect_equal(inf$get_information_source_used(), "fisher")
	expect_error(inf$set_information_preference("fisher"), "does not support information_preference")
	expect_error(inf$set_information_preference("invalid"), "information_preference must be one of")
})

test_that("public likelihood tests report an unavailable specification as non-estimable", {
	inf = make_no_likelihood_spec_logit_inference(make_public_contract_logit_design(seed = 13L))

	p_value = inf$compute_score_two_sided_pval(delta = 0)

	expect_true(is.na(p_value))
	expect_true(inf$is_nonestimable("estimate"))
	expect_equal(inf$get_nonestimable_reason(), "likelihood_test_spec_unavailable")
})
