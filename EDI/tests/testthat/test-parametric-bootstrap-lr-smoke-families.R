library(testthat)
library(EDI)

make_param_boot_weibull_censored_design <- function(seed = 20260521L, n = 100L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	linpred <- 0.25 + 0.35 * w01 + 0.25 * x1 - 0.15 * x2
	event_time <- rweibull(n, shape = 1.4, scale = exp(linpred))
	censor_time <- rexp(n, rate = 0.35)
	y <- pmin(event_time, censor_time)
	dead <- as.integer(event_time <= censor_time)

	des <- DesignFixedBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y, dead)
	des
}

make_param_boot_zinb_design <- function(seed = 20260522L, n = 120L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	p_zero <- plogis(-1.0 + 0.55 * w01 - 0.35 * x1 + 0.15 * x2)
	mu <- exp(0.35 + 0.25 * w01 + 0.25 * x1 - 0.10 * x2)
	is_zero <- rbinom(n, 1L, p_zero)
	y_count <- rnbinom(n, mu = mu, size = 2.5)
	y <- ifelse(is_zero == 1L, 0L, y_count)

	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_param_boot_zero_augmented_poisson_design <- function(seed = 20260525L, n = 80L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	p_zero <- plogis(-0.85 + 0.35 * w01 - 0.25 * x1)
	mu <- exp(0.25 + 0.15 * w01 + 0.20 * x1 - 0.10 * x2)
	is_zero <- rbinom(n, 1L, p_zero)
	y_count <- rpois(n, mu)
	y <- ifelse(is_zero == 1L, 0L, y_count)

	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_param_boot_ordinal_design <- function(seed = 20260523L, n = 120L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	eta <- 0.45 * w01 + 0.30 * x1 - 0.20 * x2
	cut_1 <- plogis(-1.10 - eta)
	cut_2 <- plogis(-0.10 - eta)
	cut_3 <- plogis(0.90 - eta)
	u <- runif(n)
	y <- ifelse(u <= cut_1, 1L, ifelse(u <= cut_2, 2L, ifelse(u <= cut_3, 3L, 4L)))

	des <- DesignFixedBernoulli$new(n = n, response_type = "ordinal", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_param_boot_count_kk_glmm_design <- function(seed = 20260524L, n = 72L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "count", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i <- des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
		mu_i <- exp(0.20 + 0.35 * w_i + 0.20 * x1[i] - 0.10 * x2[i] + rnorm(1, sd = 0.15))
		y_i <- rpois(1L, lambda = mu_i)
		des$add_one_subject_response(i, y_i, 1L)
	}
	des
}

assert_param_bootstrap_lr_smoke <- function(inf, B = 9L, min_success = 3L, seed = 9001L){
	priv <- inf$.__enclos_env__$private
	expect_true(isTRUE(priv$supports_lik_ratio_param_bootstrap()))

	inf$set_seed(seed)
	inf$num_cores <- 1L
	p_boot <- inf$compute_lik_ratio_bootstrap_two_sided_pval(
		delta = 0,
		B = B,
		show_progress = FALSE,
		min_number_usable_samples = min_success,
		max_attempts_per_replicate = 2L
	)
	diag <- inf$get_last_param_bootstrap_diagnostics()

	expect_true(is.finite(p_boot))
	expect_gte(p_boot, 0)
	expect_lte(p_boot, 1)
	expect_true(is.list(diag))
	expect_equal(diag$B, B)
	expect_equal(diag$n_success + diag$n_failure, diag$B)
	expect_equal(sum(unlist(diag$reason_counts, use.names = FALSE)), diag$B)
	expect_gte(diag$n_success, min_success)
	expect_gt(diag$success_fraction, 0)
}

test_that("raw LR parametric bootstrap is disabled for zero-augmented Poisson models", {
	des <- make_param_boot_zero_augmented_poisson_design()

	for (class_gen in list(InferenceCountZeroInflatedPoisson, InferenceCountHurdlePoisson)) {
		inf <- class_gen$new(des, model_formula = ~ x1 + x2, use_rcpp = TRUE, verbose = FALSE)
		p_boot <- inf$compute_lik_ratio_bootstrap_two_sided_pval(
			delta = 0,
			B = 9L,
			show_progress = FALSE
		)

		expect_true(is.na(p_boot))
		expect_true(inf$is_nonestimable("se"))
		expect_equal(
			inf$get_nonestimable_reason(),
			"zero_augmented_poisson_parametric_lrt_bootstrap_disabled_due_raw_lrt_miscalibration"
		)
	}
})

test_that("parametric bootstrap LR smoke tests cover censoring, mixture, ordinal, and GLMM-style families", {
	assert_param_bootstrap_lr_smoke(
		InferenceSurvivalWeibullRegr$new(
			make_param_boot_weibull_censored_design(),
			model_formula = ~ x1 + x2,
			verbose = FALSE
		),
		B = 9L,
		min_success = 3L,
		seed = 9201L
	)

	assert_param_bootstrap_lr_smoke(
		InferenceCountZeroInflatedNegBin$new(
			make_param_boot_zinb_design(),
			model_formula = ~ x1 + x2,
			use_rcpp = TRUE,
			verbose = FALSE
		),
		B = 21L,
		min_success = 2L,
		seed = 9202L
	)

	assert_param_bootstrap_lr_smoke(
		InferenceOrdinalPropOddsRegr$new(
			make_param_boot_ordinal_design(),
			model_formula = ~ x1 + x2,
			verbose = FALSE
		),
		B = 9L,
		min_success = 3L,
		seed = 9203L
	)

	assert_param_bootstrap_lr_smoke(
		InferenceCountKKGLMM$new(
			make_param_boot_count_kk_glmm_design(),
			model_formula = ~ x1 + x2,
			use_rcpp = TRUE,
			verbose = FALSE
		),
		B = 9L,
		min_success = 3L,
		seed = 9204L
	)
})
