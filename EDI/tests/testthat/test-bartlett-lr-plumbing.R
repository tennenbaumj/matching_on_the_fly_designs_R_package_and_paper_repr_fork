library(testthat)
library(EDI)

make_logit_inference <- function(seed = 1, n = 100){
	set.seed(seed)
	X <- data.frame(x1 = rnorm(n))
	y <- rbinom(n, 1, plogis(0.5 + X$x1))
	des <- DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	InferenceIncidLogRegr$new(des, verbose = FALSE)
}

# InferenceIncidLogRegr auto-inherits approx Bartlett support from
# InferenceParamBootstrap (supports_bartlett_likelihood_ratio_approx() delegates
# to supports_lik_ratio_param_bootstrap()). To exercise the "no math implemented"
# defaults generically, force that delegation off for a dummy subclass; no
# concrete family implements the exact path yet, so it is always off regardless.
make_no_param_boot_logit_inference <- function(seed = 1, n = 100){
	ext_env <- new.env(parent = globalenv())
	ext_env$R6Class <- R6::R6Class
	ext_env$InferenceIncidLogRegr <- InferenceIncidLogRegr
	evalq({
		NoParamBootLogit <- R6Class(
			"NoParamBootLogit",
			inherit = InferenceIncidLogRegr,
			private = list(
				supports_lik_ratio_param_bootstrap = function() FALSE
			)
		)
	}, envir = ext_env)
	ext_env$NoParamBootLogit$new(make_logit_inference(seed = seed, n = n)$get_design_object(), verbose = FALSE)
}

test_that("neither lik_ratio_bartlett_approx nor lik_ratio_bartlett_exact is advertised without opt-in", {
	inf <- make_no_param_boot_logit_inference(seed = 101)
	supported <- inf$get_supported_testing_types()
	expect_false("lik_ratio_bartlett_approx" %in% supported)
	expect_false("lik_ratio_bartlett_exact" %in% supported)
	expect_error(
		inf$set_testing_type("lik_ratio_bartlett_approx"),
		"does not support testing_type"
	)
	expect_error(
		inf$set_testing_type("lik_ratio_bartlett_exact"),
		"does not support testing_type"
	)
})

test_that("Bartlett wrapper methods exist and return NA when no family math is implemented", {
	inf <- make_no_param_boot_logit_inference(seed = 102)

	pval_approx <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0)
	ci_approx <- inf$compute_lik_ratio_bartlett_approx_confidence_interval(alpha = 0.05)
	pval_exact <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0)
	ci_exact <- inf$compute_lik_ratio_bartlett_exact_confidence_interval(alpha = 0.05)

	expect_true(is.na(pval_approx))
	expect_length(ci_approx, 2)
	expect_true(all(is.na(ci_approx)))

	expect_true(is.na(pval_exact))
	expect_length(ci_exact, 2)
	expect_true(all(is.na(ci_exact)))
})

test_that("normalize_testing_type accepts bartlett_approx and bartlett_exact aliases", {
	inf <- make_logit_inference(seed = 103)
	priv <- inf$.__enclos_env__$private

	expect_equal(priv$normalize_testing_type("lik_ratio_bartlett_approx"), "lik_ratio_bartlett_approx")
	expect_equal(priv$normalize_testing_type("bartlett_approx"), "lik_ratio_bartlett_approx")
	expect_equal(priv$normalize_testing_type("lr_bartlett_approx"), "lik_ratio_bartlett_approx")
	expect_equal(priv$normalize_testing_type("lrb_approx"), "lik_ratio_bartlett_approx")

	expect_equal(priv$normalize_testing_type("lik_ratio_bartlett_exact"), "lik_ratio_bartlett_exact")
	expect_equal(priv$normalize_testing_type("bartlett_exact"), "lik_ratio_bartlett_exact")
	expect_equal(priv$normalize_testing_type("lr_bartlett_exact"), "lik_ratio_bartlett_exact")
	expect_equal(priv$normalize_testing_type("lrb_exact"), "lik_ratio_bartlett_exact")
})

test_that("subclass opting in to exact Bartlett correction with factor = 1 reproduces the raw LR p-value and CI", {
	ext_env <- new.env(parent = globalenv())
	ext_env$R6Class <- R6::R6Class
	ext_env$InferenceIncidLogRegr <- InferenceIncidLogRegr

	evalq({
		BartlettExactUnityLogit <- R6Class(
			"BartlettExactUnityLogit",
			inherit = InferenceIncidLogRegr,
			public = list(
				supports_bartlett_exact = function() private$supports_bartlett_likelihood_ratio_exact(),
				bartlett_factor_exact_public = function(delta) {
					spec <- private$get_likelihood_test_spec()
					private$get_bartlett_factor_exact(spec = spec, delta = delta, full_fit = spec$full_fit, null_fit = NULL)
				}
			),
			private = list(
				supports_bartlett_likelihood_ratio_exact = function() TRUE,
				get_bartlett_factor_exact = function(spec, delta, full_fit, null_fit) 1
			)
		)
	}, envir = ext_env)

	inf <- ext_env$BartlettExactUnityLogit$new(make_logit_inference(seed = 104)$get_design_object(), verbose = FALSE)

	expect_true(inf$supports_bartlett_exact())
	expect_equal(inf$bartlett_factor_exact_public(0), 1)
	expect_true("lik_ratio_bartlett_exact" %in% inf$get_supported_testing_types())

	raw_pval <- inf$compute_lik_ratio_two_sided_pval(delta = 0)
	bartlett_pval <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0)
	expect_true(is.finite(raw_pval))
	expect_equal(bartlett_pval, raw_pval, tolerance = 1e-8)

	raw_ci <- inf$compute_lik_ratio_confidence_interval(alpha = 0.1)
	bartlett_ci <- inf$compute_lik_ratio_bartlett_exact_confidence_interval(alpha = 0.1)
	expect_true(all(is.finite(raw_ci)))
	expect_true(all(is.finite(bartlett_ci)))
	expect_equal(as.numeric(bartlett_ci), as.numeric(raw_ci), tolerance = 1e-3)
})

test_that("Exact Bartlett correction factor > 1 shrinks the LR statistic and inflates the p-value", {
	ext_env <- new.env(parent = globalenv())
	ext_env$R6Class <- R6::R6Class
	ext_env$InferenceIncidLogRegr <- InferenceIncidLogRegr

	evalq({
		BartlettExactQuadrupleLogit <- R6Class(
			"BartlettExactQuadrupleLogit",
			inherit = InferenceIncidLogRegr,
			private = list(
				supports_bartlett_likelihood_ratio_exact = function() TRUE,
				get_bartlett_factor_exact = function(spec, delta, full_fit, null_fit) 4
			)
		)
	}, envir = ext_env)

	inf <- ext_env$BartlettExactQuadrupleLogit$new(make_logit_inference(seed = 105)$get_design_object(), verbose = FALSE)

	raw_pval <- inf$compute_lik_ratio_two_sided_pval(delta = 0)
	bartlett_pval <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0)

	raw_stat <- qchisq(raw_pval, df = 1, lower.tail = FALSE)
	expected_bartlett_pval <- pchisq(raw_stat / 4, df = 1, lower.tail = FALSE)

	expect_equal(bartlett_pval, expected_bartlett_pval, tolerance = 1e-6)
	expect_gt(bartlett_pval, raw_pval)
})

test_that("compute_asymp_two_sided_pval and compute_asymp_confidence_interval dispatch to the exact Bartlett path when selected", {
	ext_env <- new.env(parent = globalenv())
	ext_env$R6Class <- R6::R6Class
	ext_env$InferenceIncidLogRegr <- InferenceIncidLogRegr

	evalq({
		BartlettExactUnityLogit2 <- R6Class(
			"BartlettExactUnityLogit2",
			inherit = InferenceIncidLogRegr,
			private = list(
				supports_bartlett_likelihood_ratio_exact = function() TRUE,
				get_bartlett_factor_exact = function(spec, delta, full_fit, null_fit) 1
			)
		)
	}, envir = ext_env)

	inf <- ext_env$BartlettExactUnityLogit2$new(make_logit_inference(seed = 106)$get_design_object(), verbose = FALSE)
	inf$set_testing_type("lik_ratio_bartlett_exact")

	expect_equal(inf$get_testing_type(), "lik_ratio_bartlett_exact")
	pval <- inf$compute_asymp_two_sided_pval(delta = 0)
	ci <- inf$compute_asymp_confidence_interval(alpha = 0.1)
	expect_true(is.finite(pval))
	expect_true(all(is.finite(ci)))
	expect_equal(pval, inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0), tolerance = 1e-8)
})

test_that("compute_asymp_two_sided_pval and compute_asymp_confidence_interval dispatch to the approx Bartlett path when selected", {
	inf <- make_logit_inference(seed = 107)
	inf$set_testing_type("lik_ratio_bartlett_approx")
	inf$set_seed(4242)

	expect_equal(inf$get_testing_type(), "lik_ratio_bartlett_approx")
	pval <- inf$compute_asymp_two_sided_pval(delta = 0)
	ci <- inf$compute_asymp_confidence_interval(alpha = 0.2)
	expect_true(is.finite(pval))
	expect_true(all(is.finite(ci)))
})

# ── "Best available" smart wrapper (compute_lik_ratio_bartlett_two_sided_pval /
# compute_lik_ratio_bartlett_confidence_interval): exact wins over approx, then
# approx, then errors if neither is supported. ──────────────────────────────

test_that("Smart wrapper prefers exact over approx and matches the exact path exactly, with no warning at default B", {
	ext_env <- new.env(parent = globalenv())
	ext_env$R6Class <- R6::R6Class
	ext_env$InferenceIncidLogRegr <- InferenceIncidLogRegr

	evalq({
		BartlettExactPreferred <- R6Class(
			"BartlettExactPreferred",
			inherit = InferenceIncidLogRegr,
			private = list(
				supports_bartlett_likelihood_ratio_exact = function() TRUE,
				get_bartlett_factor_exact = function(spec, delta, full_fit, null_fit) 2.5
			)
		)
	}, envir = ext_env)

	inf <- ext_env$BartlettExactPreferred$new(make_logit_inference(seed = 301)$get_design_object(), verbose = FALSE)

	# Both exact and approx are structurally available here (InferenceIncidLogRegr
	# auto-inherits approx); the smart wrapper must still prefer exact.
	expect_true(inf$.__enclos_env__$private$supports_bartlett_likelihood_ratio_exact())
	expect_true(inf$.__enclos_env__$private$supports_bartlett_likelihood_ratio_approx())

	exact_pval <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0)
	smart_pval <- expect_no_warning(inf$compute_lik_ratio_bartlett_two_sided_pval(delta = 0))
	expect_equal(smart_pval, exact_pval, tolerance = 1e-8)

	exact_ci <- inf$compute_lik_ratio_bartlett_exact_confidence_interval(alpha = 0.1)
	smart_ci <- expect_no_warning(inf$compute_lik_ratio_bartlett_confidence_interval(alpha = 0.1))
	expect_equal(as.numeric(smart_ci), as.numeric(exact_ci), tolerance = 1e-8)
})

test_that("Smart wrapper warns when B is explicitly supplied but the exact path is used", {
	ext_env <- new.env(parent = globalenv())
	ext_env$R6Class <- R6::R6Class
	ext_env$InferenceIncidLogRegr <- InferenceIncidLogRegr

	evalq({
		BartlettExactPreferred2 <- R6Class(
			"BartlettExactPreferred2",
			inherit = InferenceIncidLogRegr,
			private = list(
				supports_bartlett_likelihood_ratio_exact = function() TRUE,
				get_bartlett_factor_exact = function(spec, delta, full_fit, null_fit) 2.5
			)
		)
	}, envir = ext_env)

	inf <- ext_env$BartlettExactPreferred2$new(make_logit_inference(seed = 302)$get_design_object(), verbose = FALSE)

	expect_warning(
		inf$compute_lik_ratio_bartlett_two_sided_pval(delta = 0, B = 50),
		"B is ignored"
	)
	expect_warning(
		inf$compute_lik_ratio_bartlett_confidence_interval(alpha = 0.1, B = 50),
		"B is ignored"
	)
})

test_that("Smart wrapper falls back to approx when exact is unsupported", {
	inf <- make_logit_inference(seed = 303)
	priv <- inf$.__enclos_env__$private
	expect_false(priv$supports_bartlett_likelihood_ratio_exact())
	expect_true(priv$supports_bartlett_likelihood_ratio_approx())

	inf$set_seed(9999)
	approx_pval <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 99)
	inf$set_seed(9999)
	smart_pval <- expect_no_warning(inf$compute_lik_ratio_bartlett_two_sided_pval(delta = 0, B = 99))
	expect_equal(smart_pval, approx_pval, tolerance = 0)

	inf$set_seed(8888)
	approx_ci <- inf$compute_lik_ratio_bartlett_approx_confidence_interval(alpha = 0.2, B = 40)
	inf$set_seed(8888)
	smart_ci <- expect_no_warning(inf$compute_lik_ratio_bartlett_confidence_interval(alpha = 0.2, B = 40))
	expect_equal(as.numeric(smart_ci), as.numeric(approx_ci), tolerance = 0)
})

test_that("Smart wrapper errors when neither exact nor approx is supported", {
	inf <- make_no_param_boot_logit_inference(seed = 304)
	expect_error(
		inf$compute_lik_ratio_bartlett_two_sided_pval(delta = 0),
		"does not support Bartlett-corrected likelihood-ratio inference"
	)
	expect_error(
		inf$compute_lik_ratio_bartlett_confidence_interval(alpha = 0.05),
		"does not support Bartlett-corrected likelihood-ratio inference"
	)
})
