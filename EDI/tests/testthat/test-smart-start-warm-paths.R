test_that("likelihood inversion paths populate and reuse null-fit warm-start caches", {
	set.seed(202)
	n <- 60
	x <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	des <- EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	linpred <- -0.2 + 0.7 * w01 + 0.3 * x
	des$add_all_subject_responses(rbinom(n, 1, plogis(linpred)))

	inf <- InferenceIncidLogRegr$new(des, verbose = FALSE, smart_cold_start_default = TRUE)
	inf$set_testing_type("lik_ratio")
	p1 <- inf$compute_asymp_two_sided_pval(0.05)
	p2 <- inf$compute_asymp_two_sided_pval(0.08)
	cache_lr <- inf$.__enclos_env__$private$likelihood_null_warm_cache[["likelihood_test:lik_ratio"]]
	expect_true(is.finite(p1) && is.finite(p2))
	expect_true(!is.null(cache_lr$start))
	expect_true(length(cache_lr$start) > 0L)

	inf$set_testing_type("score")
	ci <- inf$compute_asymp_confidence_interval(alpha = 0.2)
	cache_score <- inf$.__enclos_env__$private$likelihood_null_warm_cache[["likelihood_test:score"]]
	expect_true(all(is.finite(ci)))
	expect_true(!is.null(cache_score$start))

	inf$set_testing_type("gradient")
	p3 <- inf$compute_asymp_two_sided_pval(0.06)
	p4 <- inf$compute_asymp_two_sided_pval(0.09)
	cache_gradient <- inf$.__enclos_env__$private$likelihood_null_warm_cache[["likelihood_test:gradient"]]
	expect_true(is.finite(p3) && is.finite(p4))
	expect_true(!is.null(cache_gradient$start))
	expect_true(length(cache_gradient$start) > 0L)

	ci_gradient <- inf$compute_asymp_confidence_interval(alpha = 0.2)
	cache_gradient_ci <- inf$.__enclos_env__$private$likelihood_null_warm_cache[["gradient_ci"]]
	expect_true(all(is.finite(ci_gradient)))
	expect_true(!is.null(cache_gradient_ci$start))
	expect_true(length(cache_gradient_ci$start) > 0L)
})

test_that("gradient testing type is available on representative likelihood families", {
	set.seed(505)
	n <- 80
	x <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2

	des_logit <- EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des_logit$add_all_subjects_to_experiment(data.frame(x = x))
	des_logit$overwrite_all_subject_assignments(w)
	des_logit$add_all_subject_responses(rbinom(n, 1, plogis(-0.1 + 0.5 * w01 + 0.2 * x)))
	inf_logit <- InferenceIncidLogRegr$new(des_logit, verbose = FALSE, smart_cold_start_default = TRUE)
	expect_true("gradient" %in% inf_logit$get_supported_testing_types())
	inf_logit$set_testing_type("gradient")
	p_logit <- inf_logit$compute_asymp_two_sided_pval(0)
	expect_true(is.finite(p_logit))
	expect_gte(p_logit, 0)
	expect_lte(p_logit, 1)

	des_pois <- EDI:::DesignFixed$new(n = n, response_type = "count", verbose = FALSE)
	des_pois$add_all_subjects_to_experiment(data.frame(x = x))
	des_pois$overwrite_all_subject_assignments(w)
	des_pois$add_all_subject_responses(rpois(n, lambda = exp(0.2 + 0.35 * w01 + 0.15 * x)))
	inf_pois <- InferenceCountPoisson$new(des_pois, verbose = FALSE, smart_cold_start_default = TRUE)
	expect_true("gradient" %in% inf_pois$get_supported_testing_types())
	inf_pois$set_testing_type("gradient")
	p_pois <- inf_pois$compute_asymp_two_sided_pval(0)
	expect_true(is.finite(p_pois))
	expect_gte(p_pois, 0)
	expect_lte(p_pois, 1)

	des_nb <- EDI:::DesignFixed$new(n = n, response_type = "count", verbose = FALSE)
	des_nb$add_all_subjects_to_experiment(data.frame(x = x))
	des_nb$overwrite_all_subject_assignments(w)
	mu_nb <- exp(0.15 + 0.3 * w01 + 0.1 * x)
	des_nb$add_all_subject_responses(rnbinom(n, mu = mu_nb, size = 1.5))
	inf_nb <- InferenceCountNegBin$new(des_nb, verbose = FALSE)
	expect_true("gradient" %in% inf_nb$get_supported_testing_types())
	inf_nb$set_testing_type("gradient")
	p_nb <- inf_nb$compute_asymp_two_sided_pval(0)
	expect_true(is.finite(p_nb))
	expect_gte(p_nb, 0)
	expect_lte(p_nb, 1)

	des_weib <- EDI:::DesignFixed$new(n = n, response_type = "survival", verbose = FALSE)
	des_weib$add_all_subjects_to_experiment(data.frame(x = x))
	des_weib$overwrite_all_subject_assignments(w)
	y_surv <- exp(0.4 + 0.3 * w01 + 0.15 * x + rnorm(n, sd = 0.2))
	dead <- rbinom(n, 1, 0.85)
	des_weib$add_all_subject_responses(y_surv, dead)
	inf_weib <- InferenceSurvivalWeibullRegr$new(des_weib, verbose = FALSE, smart_cold_start_default = TRUE)
	expect_true("gradient" %in% inf_weib$get_supported_testing_types())
	inf_weib$set_testing_type("gradient")
	p_weib <- inf_weib$compute_asymp_two_sided_pval(0)
	ci_weib <- inf_weib$compute_asymp_confidence_interval(alpha = 0.2)
	cache_weib_ci <- inf_weib$.__enclos_env__$private$likelihood_null_warm_cache[["gradient_ci"]]
	expect_true(is.finite(p_weib))
	expect_gte(p_weib, 0)
	expect_lte(p_weib, 1)
	expect_true(all(is.finite(ci_weib)))
	expect_true(is.null(cache_weib_ci) || !is.null(cache_weib_ci$start))

	inf_cox <- InferenceSurvivalCoxPHRegr$new(des_weib, verbose = FALSE)
	expect_true("gradient" %in% inf_cox$get_supported_testing_types())
	inf_cox$set_testing_type("gradient")
	p_cox <- inf_cox$compute_asymp_two_sided_pval(0)
	ci_cox <- inf_cox$compute_asymp_confidence_interval(alpha = 0.2)
	cache_cox_ci <- inf_cox$.__enclos_env__$private$likelihood_null_warm_cache[["gradient_ci"]]
	expect_true(is.finite(p_cox))
	expect_gte(p_cox, 0)
	expect_lte(p_cox, 1)
	expect_true(all(is.finite(ci_cox)))
	expect_true(!is.null(cache_cox_ci$start))
})

test_that("bootstrap reusable workers retain warm-start state across refits", {
	set.seed(303)
	n <- 50
	x <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2
	des <- EDI:::DesignFixed$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	mu <- exp(0.25 + 0.4 * w01 + 0.2 * x)
	des$add_all_subject_responses(rpois(n, mu))

	inf <- InferenceCountPoisson$new(des, verbose = FALSE, smart_cold_start_default = TRUE)
	worker_state <- inf$.__enclos_env__$private$create_bootstrap_worker_state()
	idx1 <- sample.int(n, n, replace = TRUE)
	idx2 <- sample.int(n, n, replace = TRUE)

	inf$.__enclos_env__$private$load_bootstrap_sample_into_worker(worker_state, idx1)
	val1 <- inf$.__enclos_env__$private$compute_bootstrap_worker_estimate(worker_state)
	start1 <- worker_state$worker_priv$fit_warm_start
	inf$.__enclos_env__$private$load_bootstrap_sample_into_worker(worker_state, idx2)
	val2 <- inf$.__enclos_env__$private$compute_bootstrap_worker_estimate(worker_state)
	start2 <- worker_state$worker_priv$fit_warm_start

	expect_true(is.finite(val1))
	expect_true(is.finite(val2))
	expect_true(!is.null(start1) && length(start1) > 0L)
	expect_true(!is.null(start2) && length(start2) > 0L)
})

test_that("randomization CI p-value search reuses object-level warm starts across nearby deltas", {
	set.seed(404)
	n <- 40
	x <- rnorm(n)
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x = x[i]))
	}
	add_all_subject_responses_seq(des, rpois(n, lambda = exp(0.2 + 0.3 * des$.__enclos_env__$private$w)))

	inf <- InferenceCountPoisson$new(des, verbose = FALSE, smart_cold_start_default = TRUE)
	priv <- inf$.__enclos_env__$private
	perms <- priv$generate_permutations(25L)
	cache_env <- new.env(parent = emptyenv())
	ctrl <- priv$normalize_randomization_ci_search_control(NULL, r = 25L, pval_epsilon = 0.05)

	p1 <- priv$compute_randomization_ci_pval_cached(inf, 25L, 0.05, "log", perms, ctrl, cache_env)
	start1 <- priv$fit_warm_start
	p2 <- priv$compute_randomization_ci_pval_cached(inf, 25L, 0.08, "log", perms, ctrl, cache_env)
	start2 <- priv$fit_warm_start

	expect_true(is.finite(p1) || is.na(p1))
	expect_true(is.finite(p2) || is.na(p2))
	expect_true(!is.null(start1) && length(start1) > 0L)
	expect_true(!is.null(start2) && length(start2) > 0L)
})
