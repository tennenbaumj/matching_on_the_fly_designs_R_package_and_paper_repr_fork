library(testthat)
library(EDI)

make_param_boot_logit_design <- function(seed = 20260518L, n = 120L){
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

test_that("parametric bootstrap LR reusable worker path matches generic path", {
	des <- make_param_boot_logit_design()

	inf_worker <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf_worker$set_seed(4242)
	inf_worker$num_cores <- 1L

	inf_generic <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf_generic$set_seed(4242)
	inf_generic$num_cores <- 1L
	inf_generic$.__enclos_env__$private$reusable_bootstrap_worker_enabled <- FALSE

	p_worker <- inf_worker$compute_lik_ratio_bootstrap_two_sided_pval(B = 21L, show_progress = FALSE)
	p_generic <- inf_generic$compute_lik_ratio_bootstrap_two_sided_pval(B = 21L, show_progress = FALSE)
	ci_worker <- inf_worker$compute_lik_ratio_bootstrap_confidence_interval(alpha = 0.2, B = 9L, show_progress = FALSE)
	ci_generic <- inf_generic$compute_lik_ratio_bootstrap_confidence_interval(alpha = 0.2, B = 9L, show_progress = FALSE)

	expect_equal(p_worker, p_generic, tolerance = 1e-12)
	expect_equal(ci_worker, ci_generic, tolerance = 1e-12)
	expect_true(is.finite(p_worker))
	expect_true(all(is.finite(ci_worker)))

	diag_worker <- inf_worker$get_last_param_bootstrap_diagnostics()
	expect_true(is.list(diag_worker))
	expect_equal(diag_worker$B, 9L)
	expect_true(isTRUE(diag_worker$used_reusable_worker))
	expect_true(isTRUE(diag_worker$used_deterministic_mode))
	expect_equal(diag_worker$n_success + diag_worker$n_failure, diag_worker$B)
	expect_equal(sum(unlist(diag_worker$reason_counts, use.names = FALSE)), diag_worker$B)
	expect_equal(length(diag_worker$replicate_results), diag_worker$B)
	expect_true(all(vapply(diag_worker$replicate_results, function(x) is.list(x) && !is.null(x$reason), logical(1))))
})

test_that("parametric bootstrap LR enforces a minimum usable replicate threshold", {
	des <- make_param_boot_logit_design(seed = 20260519L, n = 80L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(99)

	expect_error(
		inf$compute_lik_ratio_bootstrap_two_sided_pval(B = 4L, min_number_usable_samples = 5L, show_progress = FALSE),
		"B must be at least min_number_usable_samples"
	)
})

test_that("parametric bootstrap worker loader stores a simulated null draw", {
	des <- make_param_boot_logit_design(seed = 20260520L, n = 90L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv <- inf$.__enclos_env__$private

	spec <- priv$get_likelihood_test_spec()
	eval_obs <- priv$get_memoized_likelihood_test_eval(
		delta = 0,
		testing_type = "lik_ratio",
		spec = spec,
		include_full_negloglik = TRUE,
		include_null_negloglik = TRUE
	)
	worker_state <- priv$create_param_bootstrap_worker_state(spec, 0, eval_obs$null_fit)
	boot_spec <- priv$simulate_under_lik_null(worker_state$spec, 0, worker_state$null_fit)

	expect_true(priv$load_param_bootstrap_draw_into_worker(worker_state, boot_spec))
	expect_true(is.list(worker_state$state_env$current_param_bootstrap_draw))
	expect_true(is.function(worker_state$state_env$current_param_bootstrap_draw$fit_null))
	expect_true(is.function(worker_state$state_env$current_param_bootstrap_draw$neg_loglik))
	expect_false(is.null(worker_state$state_env$current_param_bootstrap_draw$full_fit))
})

test_that("generic parametric bootstrap LR is serial-parallel deterministic under a fixed seed", {
	des <- make_param_boot_logit_design(seed = 20260521L, n = 100L)
	if (isTRUE(EDI:::edi_env$mirai_has_been_used)) {
		set_num_cores(2L, force_mirai = TRUE)
		on.exit(unset_num_cores(), add = TRUE)
	}

	inf_serial <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf_serial$set_seed(777)
	inf_serial$num_cores <- 1L
	inf_serial$.__enclos_env__$private$reusable_bootstrap_worker_enabled <- FALSE

	inf_parallel <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf_parallel$set_seed(777)
	inf_parallel$num_cores <- 2L
	inf_parallel$.__enclos_env__$private$reusable_bootstrap_worker_enabled <- FALSE

	p_serial <- inf_serial$compute_lik_ratio_bootstrap_two_sided_pval(B = 21L, show_progress = FALSE)
	p_parallel <- inf_parallel$compute_lik_ratio_bootstrap_two_sided_pval(B = 21L, show_progress = FALSE)
	ci_serial <- inf_serial$compute_lik_ratio_bootstrap_confidence_interval(alpha = 0.2, B = 9L, show_progress = FALSE)
	ci_parallel <- inf_parallel$compute_lik_ratio_bootstrap_confidence_interval(alpha = 0.2, B = 9L, show_progress = FALSE)

	expect_equal(p_serial, p_parallel, tolerance = 0)
	expect_equal(ci_serial, ci_parallel, tolerance = 0)
	expect_true(is.finite(p_serial))
	expect_true(all(is.finite(ci_serial)))

	diag_serial <- inf_serial$get_last_param_bootstrap_diagnostics()
	diag_parallel <- inf_parallel$get_last_param_bootstrap_diagnostics()
	expect_true(isTRUE(diag_serial$used_deterministic_mode))
	expect_true(isTRUE(diag_parallel$used_deterministic_mode))
})

test_that("parametric bootstrap LR keeps deterministic mode off when no seed is set", {
	des <- make_param_boot_logit_design(seed = 20260522L, n = 90L)
	if (isTRUE(EDI:::edi_env$mirai_has_been_used)) {
		set_num_cores(2L, force_mirai = TRUE)
		on.exit(unset_num_cores(), add = TRUE)
	}
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$num_cores <- 2L
	inf$.__enclos_env__$private$reusable_bootstrap_worker_enabled <- FALSE

	p <- inf$compute_lik_ratio_bootstrap_two_sided_pval(B = 9L, show_progress = FALSE)
	diag <- inf$get_last_param_bootstrap_diagnostics()

	expect_true(is.finite(p))
	expect_true(is.list(diag))
	expect_false(isTRUE(diag$used_deterministic_mode))
})
