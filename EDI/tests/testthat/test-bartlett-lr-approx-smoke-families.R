library(testthat)
library(EDI)

# Every InferenceParamBootstrap family automatically inherits approx Bartlett
# support (supports_bartlett_likelihood_ratio_approx() delegates to
# supports_lik_ratio_param_bootstrap()), since the Monte-Carlo factor reuses the
# same simulate_under_lik_null()/refit machinery already validated for
# compute_lik_ratio_bootstrap_two_sided_pval(). This file smoke-tests that
# delegation across a representative cross-section of families (mirroring
# test-parametric-bootstrap-lr-smoke-families.R), plus the one known carve-out
# (Zero-Inflated Poisson / Hurdle Poisson, whose raw LR is itself miscalibrated).

assert_bartlett_approx_smoke <- function(inf, B = 15L, alpha = 0.25, seed = 31001L){
	priv <- inf$.__enclos_env__$private
	expect_true(isTRUE(priv$supports_lik_ratio_param_bootstrap()))
	expect_true(isTRUE(priv$supports_bartlett_likelihood_ratio_approx()))
	expect_true("lik_ratio_bartlett_approx" %in% inf$get_supported_testing_types())

	inf$num_cores <- 1L
	est <- inf$compute_estimate()
	expect_true(is.finite(est))

	inf$set_seed(seed)
	pval <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = B)
	expect_true(is.finite(pval))
	expect_gte(pval, 0)
	expect_lte(pval, 1)

	inf$set_seed(seed)
	ci <- inf$compute_lik_ratio_bartlett_approx_confidence_interval(alpha = alpha, B = B)
	expect_length(ci, 2)
	expect_true(all(is.finite(ci)))
	expect_true(ci[[1]] <= est && est <= ci[[2]])
}

make_bartlett_smoke_poisson_design <- function(seed = 41001L, n = 90L){
	set.seed(seed)
	x1 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	mu <- exp(0.3 + 0.4 * w01 + 0.2 * x1)
	y <- rpois(n, mu)

	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_bartlett_smoke_negbin_design <- function(seed = 41002L, n = 100L){
	set.seed(seed)
	x1 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	mu <- exp(0.35 + 0.35 * w01 + 0.2 * x1)
	y <- rnbinom(n, mu = mu, size = 2.5)

	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_bartlett_smoke_ols_kk_design <- function(seed = 41003L, n = 60L){
	set.seed(seed)
	x1 <- rnorm(n)
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i <- des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i]))
		y_i <- 0.5 + 0.4 * w_i + 0.3 * x1[i] + rnorm(1, sd = 0.6)
		des$add_one_subject_response(i, y_i, 1L)
	}
	des
}

make_bartlett_smoke_weibull_design <- function(seed = 41004L, n = 90L){
	set.seed(seed)
	x1 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	linpred <- 0.3 + 0.35 * w01 + 0.2 * x1
	event_time <- rweibull(n, shape = 1.3, scale = exp(linpred))
	censor_time <- rexp(n, rate = 0.3)
	y <- pmin(event_time, censor_time)
	dead <- as.integer(event_time <= censor_time)

	des <- DesignFixedBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y, dead)
	des
}

make_bartlett_smoke_zinb_design <- function(seed = 41005L, n = 110L){
	set.seed(seed)
	x1 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	p_zero <- plogis(-1.0 + 0.5 * w01 - 0.3 * x1)
	mu <- exp(0.35 + 0.25 * w01 + 0.2 * x1)
	is_zero <- rbinom(n, 1L, p_zero)
	y_count <- rnbinom(n, mu = mu, size = 2.5)
	y <- ifelse(is_zero == 1L, 0L, y_count)

	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_bartlett_smoke_zero_augmented_poisson_design <- function(seed = 41006L, n = 80L){
	set.seed(seed)
	x1 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	p_zero <- plogis(-0.85 + 0.35 * w01 - 0.25 * x1)
	mu <- exp(0.25 + 0.15 * w01 + 0.2 * x1)
	is_zero <- rbinom(n, 1L, p_zero)
	y_count <- rpois(n, mu)
	y <- ifelse(is_zero == 1L, 0L, y_count)

	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_bartlett_smoke_ordinal_design <- function(seed = 41007L, n = 110L){
	set.seed(seed)
	x1 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	eta <- 0.4 * w01 + 0.3 * x1
	cut_1 <- plogis(-1.0 - eta)
	cut_2 <- plogis(-0.05 - eta)
	cut_3 <- plogis(0.9 - eta)
	u <- runif(n)
	y <- ifelse(u <= cut_1, 1L, ifelse(u <= cut_2, 2L, ifelse(u <= cut_3, 3L, 4L)))

	des <- DesignFixedBernoulli$new(n = n, response_type = "ordinal", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

test_that("Bartlett approx smoke: InferenceCountPoisson", {
	inf <- InferenceCountPoisson$new(make_bartlett_smoke_poisson_design(), model_formula = ~ x1, verbose = FALSE)
	assert_bartlett_approx_smoke(inf, B = 15L, seed = 51001L)
})

test_that("Bartlett approx smoke: InferenceCountNegBin", {
	inf <- InferenceCountNegBin$new(make_bartlett_smoke_negbin_design(), model_formula = ~ x1, verbose = FALSE)
	assert_bartlett_approx_smoke(inf, B = 21L, seed = 51002L)
})

test_that("Bartlett approx smoke: InferenceContinKKOLSOneLik", {
	inf <- InferenceContinKKOLSOneLik$new(make_bartlett_smoke_ols_kk_design(), verbose = FALSE)
	assert_bartlett_approx_smoke(inf, B = 15L, seed = 51003L)
})

test_that("Bartlett approx smoke: InferenceSurvivalWeibullRegr", {
	inf <- InferenceSurvivalWeibullRegr$new(make_bartlett_smoke_weibull_design(), model_formula = ~ x1, verbose = FALSE)
	assert_bartlett_approx_smoke(inf, B = 9L, seed = 51004L)
})

test_that("Bartlett approx smoke: InferenceOrdinalPropOddsRegr", {
	inf <- InferenceOrdinalPropOddsRegr$new(make_bartlett_smoke_ordinal_design(), model_formula = ~ x1, verbose = FALSE)
	assert_bartlett_approx_smoke(inf, B = 9L, seed = 51005L)
})

test_that("Bartlett approx smoke: InferenceCountZeroInflatedNegBin (the ZA family that is NOT carved out)", {
	inf <- InferenceCountZeroInflatedNegBin$new(
		make_bartlett_smoke_zinb_design(),
		model_formula = ~ x1,
		use_rcpp = TRUE,
		verbose = FALSE
	)
	assert_bartlett_approx_smoke(inf, B = 21L, seed = 51006L)
})

test_that("Bartlett approx is carved out for Zero-Inflated Poisson and Hurdle Poisson (raw LR miscalibration)", {
	des <- make_bartlett_smoke_zero_augmented_poisson_design()

	for (class_gen in list(InferenceCountZeroInflatedPoisson, InferenceCountHurdlePoisson)) {
		inf <- class_gen$new(des, model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE)
		priv <- inf$.__enclos_env__$private

		# Structurally still on the parametric-bootstrap branch...
		expect_true(isTRUE(priv$supports_lik_ratio_param_bootstrap()))
		# ...but Bartlett (built on the same miscalibrated raw LR) must stay off.
		expect_false(isTRUE(priv$supports_bartlett_likelihood_ratio_approx()))
		expect_false("lik_ratio_bartlett_approx" %in% inf$get_supported_testing_types())

		inf$set_seed(51007L)
		pval <- inf$compute_lik_ratio_bartlett_approx_two_sided_pval(delta = 0, B = 9L)
		ci <- inf$compute_lik_ratio_bartlett_approx_confidence_interval(alpha = 0.25, B = 9L)

		expect_true(is.na(pval))
		expect_true(all(is.na(ci)))
	}
})
