library(testthat)
library(EDI)

test_that("g-computation design-matrix normalization supplies the treatment schema", {
	name_requests = 0L
	name_supplier = function(){
		name_requests <<- name_requests + 1L
		c("age", "risk")
	}
	X = EDI:::gcomp_normalize_treatment_design_matrix(c(1, 1, 1, 0), name_supplier)
	expect_equal(dim(X), c(2L, 2L))
	expect_identical(colnames(X), c("(Intercept)", "treatment"))
	expect_identical(name_requests, 0L)

	X_cov = EDI:::gcomp_normalize_treatment_design_matrix(
		matrix(seq_len(12L), ncol = 4L),
		covariate_names = name_supplier
	)
	expect_identical(colnames(X_cov), c("(Intercept)", "treatment", "age", "risk"))
	expect_identical(name_requests, 1L)
})

test_that("g-computation design-matrix normalization preserves existing names", {
	X = matrix(seq_len(6L), ncol = 3L)
	colnames(X) = c("intercept", "w", "x")
	expect_identical(
		EDI:::gcomp_normalize_treatment_design_matrix(X, covariate_names = "replacement"),
		X
	)
})

test_that("g-computation covariate dropping preserves intercept and treatment", {
	X = matrix(seq_len(20L), ncol = 5L)
	expect_identical(EDI:::gcomp_select_covariate_to_drop(X, c(0.1, 0.2, 0.3, -1.5, 0.8)), 4L)
	expect_identical(EDI:::gcomp_select_covariate_to_drop(X, c(0.1, 0.2, NA_real_, Inf, NA_real_)), 5L)
	expect_identical(EDI:::gcomp_select_covariate_to_drop(X[, 1:2, drop = FALSE], c(0.1, 0.2)), NA_integer_)
})

test_that("logistic potential outcomes use all-treated and all-control counterfactuals", {
	X = cbind(`(Intercept)` = 1, treatment = c(0, 1, 0), age = c(-1, 0, 2))
	coef_hat = c(-0.3, 0.8, 0.4)
	potential = EDI:::gcomp_logistic_potential_outcomes(X, coef_hat, j_treat = 2L)

	expect_identical(X[, 2L], c(0, 1, 0))
	expect_identical(potential$X1[, 2L], rep(1, 3L))
	expect_identical(potential$X0[, 2L], rep(0, 3L))
	expect_equal(potential$risk1_i, stats::plogis(as.numeric(potential$X1 %*% coef_hat)))
	expect_equal(potential$risk0_i, stats::plogis(as.numeric(potential$X0 %*% coef_hat)))
	expect_equal(potential$risk1, mean(potential$risk1_i))
	expect_equal(potential$risk0, mean(potential$risk0_i))
})

test_that("standardized-effect cache readiness separates point and inference caches", {
	point_only = list(rd = 0.1, se_rd = NA_real_, rr = 1.2, se_log_rr = NA_real_)
	expect_true(EDI:::gcomp_standardized_effect_cache_is_ready(point_only, estimate_only = TRUE))
	expect_false(EDI:::gcomp_standardized_effect_cache_is_ready(point_only, estimate_only = FALSE))

	full_attempt = point_only
	full_attempt$gcomp_standardized_effects_inference_ready = TRUE
	expect_true(EDI:::gcomp_standardized_effect_cache_is_ready(full_attempt, estimate_only = FALSE))
})

test_that("standardized-effect cache mapper updates incidence fields narrowly", {
	effects = list(
		summary_table = matrix(1, nrow = 1L),
		full_coefficients = c(`(Intercept)` = -0.1, treatment = 0.5),
		full_vcov = diag(2),
		risk1 = 0.6,
		risk0 = 0.4,
		rd = 0.2,
		se_rd = 0.05,
		log_rr = log(1.5),
		rr = 1.5,
		se_log_rr = 0.1
	)
	cache = EDI:::gcomp_cache_standardized_effects(
		list(rand_distr_cache = list(existing = TRUE)),
		effects,
		inference_ready = TRUE
	)

	expect_equal(cache$rd, effects$rd)
	expect_equal(cache$rr, effects$rr)
	expect_equal(cache$se_rd, effects$se_rd)
	expect_equal(cache$se_log_rr, effects$se_log_rr)
	expect_true(cache$gcomp_standardized_effects_inference_ready)
	expect_identical(cache$rand_distr_cache, list(existing = TRUE))

	failed = EDI:::gcomp_cache_failed_standardized_effects(cache, inference_ready = FALSE)
	expect_true(is.na(failed$rd))
	expect_true(is.na(failed$rr))
	expect_false(failed$gcomp_standardized_effects_inference_ready)
	expect_identical(failed$rand_distr_cache, list(existing = TRUE))
})
