library(testthat)
library(EDI)

reference_normalize_treatment_design_matrix = function(X, covariate_names = NULL){
	if (is.null(dim(X))) {
		X = matrix(X, ncol = 2L)
	}
	if (is.null(colnames(X))) {
		colnames(X) = c(
			"(Intercept)",
			"treatment",
			if (ncol(X) > 2L) covariate_names else NULL
		)
	}
	X
}

reference_select_covariate_to_drop = function(X_curr, coef_hat){
	covariate_cols = seq.int(3L, ncol(X_curr))
	if (length(covariate_cols) == 0L) return(NA_integer_)
	coef_mags = abs(coef_hat[covariate_cols])
	if (length(coef_mags) == 0L || all(!is.finite(coef_mags))){
		return(tail(covariate_cols, 1L))
	}
	covariate_cols[which.max(replace(coef_mags, !is.finite(coef_mags), -Inf))]
}

reference_logistic_potential_outcomes = function(X_fit, coef_hat, j_treat){
	X1 = X_fit
	X0 = X_fit
	X1[, j_treat] = 1
	X0[, j_treat] = 0
	risk1_i = stats::plogis(as.numeric(X1 %*% coef_hat))
	risk0_i = stats::plogis(as.numeric(X0 %*% coef_hat))
	list(
		X1 = X1,
		X0 = X0,
		risk1_i = risk1_i,
		risk0_i = risk0_i,
		risk1 = mean(risk1_i),
		risk0 = mean(risk0_i)
	)
}

reference_cache_standardized_effects = function(cached_values, effects, inference_ready = FALSE){
	cached_values$summary_table = effects$summary_table
	cached_values$full_coefficients = effects$full_coefficients
	cached_values$full_vcov = effects$full_vcov
	cached_values$risk1 = effects$risk1
	cached_values$risk0 = effects$risk0
	cached_values$rd = effects$rd
	cached_values$se_rd = effects$se_rd
	cached_values$log_rr = effects$log_rr
	cached_values$rr = effects$rr
	cached_values$se_log_rr = effects$se_log_rr
	cached_values$gcomp_standardized_effects_inference_ready = isTRUE(inference_ready)
	cached_values
}

test_that("shared design-matrix normalization matches former inline logic by g-computation family", {
	family_cases = list(
		incidence = list(
			X = matrix(seq_len(20L), ncol = 4L),
			covariates = c("age", "risk")
		),
		kk_incidence = list(
			X = matrix(seq_len(18L), ncol = 3L),
			covariates = "baseline"
		),
		proportion = list(
			X = c(1, 1, 1, 0),
			covariates = character(0)
		)
	)

	for (family in names(family_cases)) {
		case = family_cases[[family]]
		expect_equal(
			EDI:::gcomp_normalize_treatment_design_matrix(case$X, case$covariates),
			reference_normalize_treatment_design_matrix(case$X, case$covariates),
			info = family
		)
	}
})

test_that("shared covariate-drop helper matches former inline logic by g-computation family", {
	X = matrix(seq_len(30L), ncol = 5L)
	coef_cases = list(
		incidence = c(0.1, 0.4, -0.8, 1.6, 0.3),
		kk_incidence = c(0.1, 0.4, NA_real_, Inf, -0.5),
		proportion = c(0.1, 0.4, NA_real_, NA_real_, NA_real_)
	)

	for (family in names(coef_cases)) {
		coef_hat = coef_cases[[family]]
		expect_identical(
			EDI:::gcomp_select_covariate_to_drop(X, coef_hat),
			reference_select_covariate_to_drop(X, coef_hat),
			info = family
		)
	}
})

test_that("incidence and KK incidence logistic potential-outcome helper matches former inline math", {
	X = cbind(
		`(Intercept)` = 1,
		age = c(-1.5, -0.5, 0.25, 1.25),
		treatment = c(0, 1, 0, 1),
		risk = c(0.2, -0.1, 0.4, 0.7)
	)
	coef_hat = c(`(Intercept)` = -0.25, age = 0.15, treatment = 0.6, risk = -0.35)

	for (family in c("incidence", "kk_incidence")) {
		actual = EDI:::gcomp_logistic_potential_outcomes(X, coef_hat, j_treat = 3L)
		expected = reference_logistic_potential_outcomes(X, coef_hat, j_treat = 3L)
		expect_equal(actual, expected, tolerance = 1e-15, info = family)
	}
})

test_that("incidence and KK incidence cache-map helper matches former inline assignments", {
	effects = list(
		summary_table = cbind(Value = c(-0.2, 0.5), `Std. Error` = c(0.1, 0.2)),
		full_coefficients = c(`(Intercept)` = -0.2, treatment = 0.5),
		full_vcov = diag(c(0.01, 0.04)),
		risk1 = 0.65,
		risk0 = 0.45,
		rd = 0.2,
		se_rd = 0.04,
		log_rr = log(0.65 / 0.45),
		rr = 0.65 / 0.45,
		se_log_rr = 0.08
	)
	initial_cache = list(
		rand_distr_cache = list(draws = 1:3),
		m_cache = list(existing = TRUE)
	)

	for (family in c("incidence", "kk_incidence")) {
		actual = EDI:::gcomp_cache_standardized_effects(initial_cache, effects, inference_ready = TRUE)
		expected = reference_cache_standardized_effects(initial_cache, effects, inference_ready = TRUE)
		expect_equal(actual, expected, info = family)
	}
})
