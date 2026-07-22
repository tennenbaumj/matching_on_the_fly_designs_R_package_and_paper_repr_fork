#' Internal G-Computation Design-Matrix Helpers
#'
#' Normalizes a g-computation treatment design matrix without accessing an
#' inference object. G-computation families use the same two leading columns:
#' an intercept followed by treatment. Existing column names are preserved.
#'
#' @keywords internal
#' @noRd
gcomp_normalize_treatment_design_matrix = function(X, covariate_names = NULL){
	if (is.null(dim(X))) {
		X = matrix(X, ncol = 2L)
	}
	if (is.null(colnames(X))) {
		if (ncol(X) > 2L && is.function(covariate_names)) {
			covariate_names = covariate_names()
		}
		colnames(X) = c(
			"(Intercept)",
			"treatment",
			if (ncol(X) > 2L) covariate_names else NULL
		)
	}
	X
}

#' Selects the covariate column to drop from a g-computation fit retry.
#'
#' The intercept and treatment columns are always retained in positions one and
#' two. Among covariates, the largest finite coefficient magnitude is dropped;
#' if no coefficient is usable, the last covariate is dropped deterministically.
#'
#' @keywords internal
#' @noRd
gcomp_select_covariate_to_drop = function(X_curr, coef_hat){
	if (ncol(X_curr) <= 2L) return(NA_integer_)
	covariate_cols = seq.int(3L, ncol(X_curr))
	coef_mags = abs(coef_hat[covariate_cols])
	if (length(coef_mags) == 0L || all(!is.finite(coef_mags))){
		return(tail(covariate_cols, 1L))
	}
	covariate_cols[which.max(replace(coef_mags, !is.finite(coef_mags), -Inf))]
}

#' Computes logistic potential outcomes for a treatment design matrix.
#'
#' Returns treated and control counterfactual design matrices, subject-level
#' risks, and their empirical means. Variance and estimand-specific logic stay
#' with the calling inference family.
#'
#' @keywords internal
#' @noRd
gcomp_logistic_potential_outcomes = function(X_fit, coef_hat, j_treat){
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

#' Checks whether standardized g-computation effects are cache-ready.
#'
#' Estimate-only readiness requires a cached point estimate. Full-inference
#' readiness is tracked explicitly because estimate-only paths intentionally
#' cache \code{NA} standard errors.
#'
#' @keywords internal
#' @noRd
gcomp_standardized_effect_cache_is_ready = function(cached_values, estimate_only = FALSE){
	if (isTRUE(estimate_only)){
		return(!is.null(cached_values$rd) || !is.null(cached_values$rr))
	}
	isTRUE(cached_values$gcomp_standardized_effects_inference_ready)
}

#' Maps standardized g-computation effects into the inference cache.
#'
#' This helper is deliberately scoped to incidence g-computation classes, which
#' share the same standardized risk, RD, RR, and delta-method cache fields.
#'
#' @keywords internal
#' @noRd
gcomp_cache_standardized_effects = function(cached_values, effects, inference_ready = FALSE){
	for (field in c(
		"summary_table",
		"full_coefficients",
		"full_vcov",
		"risk1",
		"risk0",
		"rd",
		"se_rd",
		"log_rr",
		"rr",
		"se_log_rr"
	)){
		cached_values[[field]] = effects[[field]]
	}
	cached_values$gcomp_standardized_effects_inference_ready = isTRUE(inference_ready)
	cached_values
}

#' Maps a failed standardized g-computation fit into the inference cache.
#'
#' @keywords internal
#' @noRd
gcomp_cache_failed_standardized_effects = function(cached_values, inference_ready = TRUE){
	gcomp_cache_standardized_effects(
		cached_values = cached_values,
		effects = list(
			summary_table = NULL,
			full_coefficients = NULL,
			full_vcov = NULL,
			risk1 = NULL,
			risk0 = NULL,
			rd = NA_real_,
			se_rd = NA_real_,
			log_rr = NA_real_,
			rr = NA_real_,
			se_log_rr = NA_real_
		),
		inference_ready = inference_ready
	)
}
