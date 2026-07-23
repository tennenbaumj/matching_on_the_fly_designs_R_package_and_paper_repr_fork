library(testthat)
library(EDI)

make_param_boot_logit_design <- function(seed = 20260723L, n = 150L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	p <- plogis(-0.4 + 0.8 * w01 + 0.35 * x1 - 0.25 * x2)
	y <- rbinom(n, 1, p)

	des <- EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

test_that("compute_param_bootstrap_estimate is finite and deterministic under a fixed seed", {
	des <- make_param_boot_logit_design()

	inf1 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf1$set_seed(555)
	inf1$num_cores <- 1L

	inf2 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf2$set_seed(555)
	inf2$num_cores <- 1L

	inf3 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf3$set_seed(556)
	inf3$num_cores <- 1L

	est1 <- inf1$compute_param_bootstrap_estimate(B = 25L, show_progress = FALSE)
	est2 <- inf2$compute_param_bootstrap_estimate(B = 25L, show_progress = FALSE)
	est3 <- inf3$compute_param_bootstrap_estimate(B = 25L, show_progress = FALSE)

	expect_true(is.finite(est1))
	expect_equal(est1, est2, tolerance = 0)
	expect_false(isTRUE(all.equal(est1, est3)))
})

test_that("compute_param_bootstrap_estimate satisfies the bias-correction reconciliation identity", {
	des <- make_param_boot_logit_design(seed = 20260724L, n = 130L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(909)
	inf$num_cores <- 1L

	est <- inf$compute_param_bootstrap_estimate(B = 31L, show_progress = FALSE)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.list(diag))
	expect_equal(diag$B, 31L)
	expect_true(is.finite(diag$raw_estimate))
	expect_equal(diag$n_success + diag$n_failure, diag$B)
	expect_true(diag$n_success >= 1L)

	finite_reps <- diag$replicate_estimates[is.finite(diag$replicate_estimates)]
	expected <- 2 * diag$raw_estimate - mean(finite_reps)
	expect_equal(est, expected, tolerance = 1e-12)
})

test_that("compute_param_bootstrap_estimate errors when the family does not support it", {
	des <- make_param_boot_logit_design(seed = 20260725L, n = 90L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv_env <- inf$.__enclos_env__$private
	unlockBinding("supports_lik_ratio_param_bootstrap", priv_env)
	assign("supports_lik_ratio_param_bootstrap", function() FALSE, envir = priv_env)
	lockBinding("supports_lik_ratio_param_bootstrap", priv_env)

	expect_error(
		inf$compute_param_bootstrap_estimate(B = 10L, show_progress = FALSE),
		"does not support parametric-bootstrap"
	)
})

test_that("compute_param_bootstrap_confidence_interval returns a finite, ordered, reproducible interval", {
	des <- make_param_boot_logit_design(seed = 20260726L, n = 140L)

	inf1 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf1$set_seed(321)
	inf1$num_cores <- 1L

	inf2 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf2$set_seed(321)
	inf2$num_cores <- 1L

	ci1 <- inf1$compute_param_bootstrap_confidence_interval(alpha = 0.2, B = 41L, show_progress = FALSE)
	ci2 <- inf2$compute_param_bootstrap_confidence_interval(alpha = 0.2, B = 41L, show_progress = FALSE)

	expect_true(all(is.finite(ci1)))
	expect_length(ci1, 2)
	expect_true(ci1[[1]] <= ci1[[2]])
	expect_equal(ci1, ci2, tolerance = 0)
})

test_that("compute_param_bootstrap_confidence_interval satisfies the reflected-quantile reconciliation identity", {
	des <- make_param_boot_logit_design(seed = 20260727L, n = 130L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(707)
	inf$num_cores <- 1L

	ci <- inf$compute_param_bootstrap_confidence_interval(alpha = 0.1, B = 51L, show_progress = FALSE)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.list(diag))
	expect_equal(diag$B, 51L)

	finite_reps <- diag$replicate_estimates[is.finite(diag$replicate_estimates)]
	expected_lo <- 2 * diag$raw_estimate - stats::quantile(finite_reps, 0.95, names = FALSE)
	expected_hi <- 2 * diag$raw_estimate - stats::quantile(finite_reps, 0.05, names = FALSE)

	expect_equal(unname(ci[[1]]), expected_lo, tolerance = 1e-12)
	expect_equal(unname(ci[[2]]), expected_hi, tolerance = 1e-12)
})

test_that("compute_param_bootstrap_confidence_interval errors when the family does not support it", {
	des <- make_param_boot_logit_design(seed = 20260728L, n = 90L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv_env <- inf$.__enclos_env__$private
	unlockBinding("supports_lik_ratio_param_bootstrap", priv_env)
	assign("supports_lik_ratio_param_bootstrap", function() FALSE, envir = priv_env)
	lockBinding("supports_lik_ratio_param_bootstrap", priv_env)

	expect_error(
		inf$compute_param_bootstrap_confidence_interval(alpha = 0.05, B = 10L, show_progress = FALSE),
		"does not support parametric-bootstrap"
	)
})

test_that("compute_param_bootstrap_pval is a finite probability, deterministic under a fixed seed", {
	des <- make_param_boot_logit_design(seed = 20260729L, n = 140L)

	inf1 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf1$set_seed(111)
	inf1$num_cores <- 1L

	inf2 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf2$set_seed(111)
	inf2$num_cores <- 1L

	p1 <- inf1$compute_param_bootstrap_pval(delta = 0, B = 41L, show_progress = FALSE)
	p2 <- inf2$compute_param_bootstrap_pval(delta = 0, B = 41L, show_progress = FALSE)

	expect_true(is.finite(p1))
	expect_true(p1 >= 0 && p1 <= 1)
	expect_equal(p1, p2, tolerance = 0)
})

test_that("compute_param_bootstrap_pval satisfies the reflected-empirical-CDF reconciliation identity", {
	des <- make_param_boot_logit_design(seed = 20260730L, n = 130L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(222)
	inf$num_cores <- 1L

	delta <- 0.1
	p <- inf$compute_param_bootstrap_pval(delta = delta, B = 51L, show_progress = FALSE)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.list(diag))
	expect_equal(diag$B, 51L)

	finite_reps <- diag$replicate_estimates[is.finite(diag$replicate_estimates)]
	t_delta <- 2 * diag$raw_estimate - delta
	n <- length(finite_reps)
	left_tail <- (1 + sum(finite_reps <= t_delta)) / (1 + n)
	right_tail <- (1 + sum(finite_reps >= t_delta)) / (1 + n)
	expected <- min(1, 2 * min(left_tail, right_tail))

	expect_equal(p, expected, tolerance = 1e-12)
})

test_that("compute_param_bootstrap_pval errors when the family does not support it", {
	des <- make_param_boot_logit_design(seed = 20260731L, n = 90L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv_env <- inf$.__enclos_env__$private
	unlockBinding("supports_lik_ratio_param_bootstrap", priv_env)
	assign("supports_lik_ratio_param_bootstrap", function() FALSE, envir = priv_env)
	lockBinding("supports_lik_ratio_param_bootstrap", priv_env)

	expect_error(
		inf$compute_param_bootstrap_pval(delta = 0, B = 10L, show_progress = FALSE),
		"does not support parametric-bootstrap"
	)
})
