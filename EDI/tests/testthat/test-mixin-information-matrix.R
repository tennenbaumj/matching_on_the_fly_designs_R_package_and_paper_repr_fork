library(testthat)
library(EDI)

make_information_matrix_logit_inference <- function(seed = 1L, n = 80L){
	set.seed(seed)
	x = rnorm(n)
	w = rep(c(1, -1), length.out = n)
	des = EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(rbinom(n, 1L, plogis(-0.2 + 0.5 * ((w + 1) / 2) + 0.3 * x)))
	InferenceIncidLogRegr$new(des, verbose = FALSE)
}

test_that("information-matrix mixin is composed into likelihood inference classes", {
	expect_true(is.list(EDI:::InferenceMixinInformationMatrix))
	expect_named(
		EDI:::InferenceMixinInformationMatrix$private,
		c(
			"get_information_matrix",
			"compute_variance_from_information_matrix",
			"compute_standard_error_from_information_matrix",
			"get_score_test_information_matrix"
		)
	)

	inf = make_information_matrix_logit_inference()
	priv = inf$.__enclos_env__$private
	expect_true(all(names(EDI:::InferenceMixinInformationMatrix$private) %in% names(priv)))
})

test_that("information-matrix mixin selects the requested source and derives standard errors", {
	inf = make_information_matrix_logit_inference(seed = 2L)
	priv = inf$.__enclos_env__$private
	fisher = diag(c(2, 4))
	observed = diag(c(3, 9))
	spec = list(full_fit = list(fisher_information = fisher, observed_information = observed), j = 2L)

	priv$information_preference = "fisher"
	expect_equal(priv$get_information_matrix(spec = spec), fisher)
	expect_equal(inf$get_information_source_used(), "fisher")

	priv$information_preference = "observed"
	expect_equal(priv$get_information_matrix(spec = spec), observed)
	expect_equal(inf$get_information_source_used(), "observed")

	priv$information_preference = "fisher"
	expect_equal(priv$compute_variance_from_information_matrix(fisher, 2L), 0.25)
	expect_equal(priv$compute_standard_error_from_information_matrix(spec = spec), 0.5)
	expect_true(is.na(priv$compute_variance_from_information_matrix(fisher, 3L)))
})
