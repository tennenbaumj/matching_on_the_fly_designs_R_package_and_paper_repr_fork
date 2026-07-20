.libPaths(c("Rlib", .libPaths()))
library(testthat)
library(EDI)
library(data.table)

make_exact_fisher_2x2 <- function(w, y) {
	matrix(
		c(
			sum(w == 1L & y == 1L),
			sum(w == 1L & y == 0L),
			sum(w == 0L & y == 1L),
			sum(w == 0L & y == 0L)
		),
		nrow = 2,
		byrow = TRUE
	)
}

format_exact_fisher_tables <- function(table_list) {
	table_list <- Filter(function(tab) sum(tab[1, ]) > 0L && sum(tab[2, ]) > 0L, table_list)
	if (length(table_list) == 1L) {
		return(table_list[[1]])
	}
	table_array <- array(0, dim = c(2L, 2L, length(table_list)))
	for (k in seq_along(table_list)) {
		table_array[, , k] <- table_list[[k]]
	}
	table_array
}

build_blocking_exact_fisher_tables <- function(des) {
	des_priv <- des$.__enclos_env__$private
	w <- des_priv$w
	y <- des_priv$y
	strata_cols <- des_priv$strata_cols
	if (is.null(strata_cols) || length(strata_cols) == 0L) {
		return(make_exact_fisher_2x2(w, y))
	}
	Xraw <- des_priv$Xraw
	strata_keys <- vapply(seq_along(w), function(i) {
		paste(vapply(strata_cols, function(col) {
			val <- Xraw[i, ][[col]]
			if (is.na(val)) "NA" else as.character(val)
		}, character(1)), collapse = "|")
	}, character(1))
	format_exact_fisher_tables(lapply(split(seq_along(w), strata_keys), function(idx) {
		make_exact_fisher_2x2(w[idx], y[idx])
	}))
}

build_kk_exact_fisher_tables <- function(des) {
	des_priv <- des$.__enclos_env__$private
	w <- des_priv$w
	y <- des_priv$y
	m <- as.integer(des_priv$m)
	m[is.na(m)] <- 0L
	table_list <- lapply(sort(unique(m[m > 0L])), function(match_id) {
		idx <- which(m == match_id)
		make_exact_fisher_2x2(w[idx], y[idx])
	})
	reservoir_idx <- which(m == 0L)
	if (length(reservoir_idx) > 0L) {
		table_list[[length(table_list) + 1L]] <- make_exact_fisher_2x2(w[reservoir_idx], y[reservoir_idx])
	}
	format_exact_fisher_tables(table_list)
}

test_that("compute_rand_confidence_interval works for continuous response", {
	set.seed(123)
	n <- 40
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", verbose = TRUE)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.table(x1 = rnorm(1), x2 = rnorm(1)))
	}
	# Add treatment effect of 1.0
	treatment <- des$.__enclos_env__$private$w
	y <- rnorm(n) + treatment * 1.0
	add_all_subject_responses_seq(des, y)

	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = TRUE)

	# Compute randomization CI
	# Using small nsim for speed in tests
	ci <- inf$compute_rand_confidence_interval(alpha = 0.05, r = 100, pval_epsilon = 0.05)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	# The estimate should be within the CI
	est <- inf$compute_estimate()
	expect_true(est >= ci[1] && est <= ci[2])

	message("Continuous Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("compute_rand_confidence_interval works for proportion response", {
	set.seed(123)
	n <- 100
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "proportion", verbose = FALSE)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.table(x1 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	# Simulate proportions
	mu <- plogis(-0.5 + treatment * 1.0)
	y <- rbeta(n, shape1 = mu * 10, shape2 = (1 - mu) * 10)
	add_all_subject_responses_seq(des, y)

	inf <- InferencePropBetaRegr$new(des, verbose = FALSE)

	# Compute randomization CI
	ci <- inf$compute_rand_confidence_interval(alpha = 0.05, r = 100, pval_epsilon = 0.05)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	expect_true(all(is.finite(ci)))

	est <- inf$compute_estimate()
	expect_true(est >= ci[1] && est <= ci[2])

	message("Proportion Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("proportion quantile randomization CI shifts on the logit scale once", {
	set.seed(124)
	n <- 48
	des <- DesignSeqOneByOneiBCRD$new(n = n, response_type = "proportion", verbose = FALSE)
	x <- rnorm(n)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.table(x1 = x[i]))
	}
	w <- des$.__enclos_env__$private$w
	mu <- plogis(-0.25 + 0.7 * w + 0.2 * x)
	y <- rbeta(n, shape1 = mu * 12, shape2 = (1 - mu) * 12)
	add_all_subject_responses_seq(des, y)

	inf <- InferencePropQuantileRegr$new(des, model_formula = ~ x1, verbose = FALSE)
	ci <- inf$compute_rand_confidence_interval(
		alpha = 0.10,
		r = 31,
		pval_epsilon = 0.10,
		show_progress = FALSE,
		ci_search_control = list(max_expansions = 2L, mc_enable = FALSE)
	)

	expect_equal(length(ci), 2)
	expect_true(all(is.finite(ci)))
	expect_true(ci[1] < ci[2])
})

test_that("compute_rand_confidence_interval works for survival response (uncensored)", {
	set.seed(123)
	n <- 50
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.table(x1 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	# Simulate survival times (log-normal)
	y <- exp(1.0 + treatment * 0.8 + rnorm(n, 0, 0.5))
	add_all_subject_responses_seq(des, y, deads = rep(1, n)) # All events, no censoring

	inf <- InferenceSurvivalWeibullRegr$new(des, verbose = FALSE)

	# Compute randomization CI
	ci <- inf$compute_rand_confidence_interval(alpha = 0.05, r = 100, pval_epsilon = 0.05)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	# Survival times (ratios) should be positive
	expect_true(all(ci > 0))

	est <- inf$compute_estimate()
	expect_true(est >= ci[1] && est <= ci[2])

	message("Survival Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("compute_rand_confidence_interval works for ordinal response (cumulative logit)", {
	set.seed(456)
	n <- 40
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "ordinal", verbose = FALSE)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.table(x1 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	score <- rnorm(n, mean = treatment * 0.5)
	y <- 1L + (score > 0) + (score > 1)
	add_all_subject_responses_seq(des, y)

	inf <- InferenceOrdinalPropOddsRegr$new(des, verbose = FALSE)

	ci <- inf$compute_rand_confidence_interval(alpha = 0.05, r = 100, pval_epsilon = 0.05, show_progress = FALSE)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	est <- inf$compute_estimate()
	expect_true(est >= ci[1] && est <= ci[2])
	expect_true(all(is.finite(ci)))

	message("Ordinal Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("compute_rand_confidence_interval throws error for unsupported types", {
	n <- 20
	des_incid <- DesignSeqOneByOneEfron$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in 1:n) des_incid$add_one_subject_to_experiment_and_assign(data.table(x=1))
	add_all_subject_responses_seq(des_incid, rbinom(n, 1, 0.5))
	inf_incid <- InferenceIncidLogRegr$new(des_incid)
	expect_error(inf_incid$compute_rand_confidence_interval(), "Randomization confidence intervals are not supported for incidence")

	des_count <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	for (i in 1:n) des_count$add_one_subject_to_experiment_and_assign(data.table(x=1))
	add_all_subject_responses_seq(des_count, rpois(n, 5))
	inf_count <- InferenceCountNegBin$new(des_count)
	ci_count <- inf_count$compute_rand_confidence_interval(alpha = 0.05, r = 100, pval_epsilon = 0.05)
	expect_equal(length(ci_count), 2)
	expect_true(ci_count[1] < ci_count[2])
	expect_true(all(is.finite(ci_count)))
})

test_that("FixedRerandomization incidence randomization uses design draws, not Zhang", {
	set.seed(123)
	n <- 20
	X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des <- DesignFixedRerandomization$new(
		n = n,
		response_type = "incidence",
		obj_val_cutoff = 100,
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(rbinom(n, 1, 0.5))

	old_asserts <- EDI:::should_run_asserts()
	on.exit(toggle_asserts(old_asserts), add = TRUE)
	toggle_asserts(TRUE)

	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	expect_no_error(p <- inf$compute_rand_two_sided_pval(r = 11, show_progress = FALSE))
	expect_true(is.finite(p))
	expect_error(
		inf$compute_rand_two_sided_pval(r = 11, type = "Zhang", show_progress = FALSE),
		"Randomization type dispatch"
	)
})

test_that("Zhang incidence inference is available through randomization and exact APIs", {
	set.seed(321)
	n <- 24
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.table(x1 = rnorm(1), x2 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	prob <- plogis(-0.2 + 0.8 * treatment)
	add_all_subject_responses_seq(des, rbinom(n, 1, prob))

	inf_rand_serial <- InferenceIncidLogRegr$new(des, verbose = FALSE)
	inf_rand_parallel <- InferenceIncidLogRegr$new(des, verbose = FALSE)
	inf_serial <- InferenceIncidenceExactZhang$new(des, verbose = FALSE)
	inf_parallel <- InferenceIncidenceExactZhang$new(des, verbose = FALSE)

	ci_rand_serial <- inf_rand_serial$compute_rand_confidence_interval(alpha = 0.10, pval_epsilon = 0.01, show_progress = FALSE)
	ci_rand_parallel <- inf_rand_parallel$compute_rand_confidence_interval(alpha = 0.10, pval_epsilon = 0.01, show_progress = FALSE)
	ci_exact_serial <- inf_serial$compute_exact_confidence_interval(alpha = 0.10, pval_epsilon = 0.01)
	ci_exact_parallel <- inf_parallel$compute_exact_confidence_interval(alpha = 0.10, pval_epsilon = 0.01)
	p_rand <- inf_rand_serial$compute_rand_two_sided_pval(delta = 0)
	p_exact <- inf_serial$compute_exact_two_sided_pval_for_treatment_effect(delta = 0)

	expect_length(ci_rand_serial, 2)
	expect_true(all(is.finite(ci_rand_serial)))
	expect_equal(ci_rand_parallel, ci_rand_serial, tolerance = 1e-8)
	expect_equal(ci_exact_serial, ci_rand_serial, tolerance = 1e-8)
	expect_equal(ci_exact_parallel, ci_rand_serial, tolerance = 1e-8)
	expect_true(is.finite(p_rand))
	expect_equal(p_exact, p_rand, tolerance = 1e-12)
})

test_that("Fisher exact inference matches fisher.test for iBCRD incidence", {
	set.seed(2026)
	n <- 20
	des <- DesignSeqOneByOneiBCRD$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.table(x1 = rnorm(1)))
	}
	w <- des$.__enclos_env__$private$w
	y <- rbinom(n, 1, plogis(-0.3 + 0.7 * w))
	add_all_subject_responses_seq(des, y)

	ref <- stats::fisher.test(make_exact_fisher_2x2(w, y), conf.level = 0.90)
	inf_exact <- InferenceIncidExactFisher$new(des, verbose = FALSE)

	ci_exact <- inf_exact$compute_exact_confidence_interval(alpha = 0.10)
	p_exact <- inf_exact$compute_exact_two_sided_pval_for_treatment_effect(delta = 0)

	expect_equal(inf_exact$compute_estimate(), log(as.numeric(ref$estimate)), tolerance = 1e-12)
	expect_equal(unname(ci_exact), log(as.numeric(ref$conf.int)), tolerance = 1e-12)
	expect_equal(p_exact, ref$p.value, tolerance = 1e-12)
})

test_that("Fisher exact inference matches mantelhaen.test for blocking incidence", {
	set.seed(2027)
	n <- 16
	x_dat <- data.table(site = rep(c("A", "B"), each = n / 2), x1 = rnorm(n))
	des <- DesignSeqOneByOneSPBR$new(
		strata_cols = "site",
		block_size = 4,
		n = n,
		response_type = "incidence",
		verbose = FALSE
	)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(x_dat[i, ])
	}
	w <- des$.__enclos_env__$private$w
	y <- rbinom(n, 1, plogis(-0.5 + 0.9 * w + 0.4 * (x_dat$site == "B")))
	add_all_subject_responses_seq(des, y)

	ref_tables <- build_blocking_exact_fisher_tables(des)
	ref <- stats::mantelhaen.test(ref_tables, exact = TRUE, conf.level = 0.95)
	inf_exact <- InferenceIncidExactFisher$new(des, verbose = FALSE)

	ci_exact <- inf_exact$compute_exact_confidence_interval(alpha = 0.05)
	p_exact <- inf_exact$compute_exact_two_sided_pval_for_treatment_effect(delta = 0)

	expect_equal(inf_exact$compute_estimate(), log(as.numeric(ref$estimate)), tolerance = 1e-12)
	expect_equal(unname(ci_exact), log(as.numeric(ref$conf.int)), tolerance = 1e-12)
	expect_equal(p_exact, ref$p.value, tolerance = 1e-12)
	expect_error(
		inf_exact$compute_exact_two_sided_pval_for_treatment_effect(delta = 0.2),
		"Stratified Fisher exact inference only supports delta = 0"
	)
})

test_that("Fisher exact inference matches mantelhaen.test for KK incidence", {
	set.seed(2028)
	x_dat <- data.table(x1 = c(-3, 3, -3.01, 3.01, 0, 0.01))
	des <- DesignSeqOneByOneKK14$new(
		n = nrow(x_dat),
		response_type = "incidence",
		t_0_pct = 0.34,
		lambda = 0.99,
		verbose = FALSE
	)
	for (i in seq_len(nrow(x_dat))) {
		des$add_one_subject_to_experiment_and_assign(x_dat[i, ])
	}
	m <- des$.__enclos_env__$private$m
	expect_true(any(m > 0L))

	w <- des$.__enclos_env__$private$w
	y <- c(1L, 0L, 0L, 1L, 1L, 0L)
	add_all_subject_responses_seq(des, y)

	ref_tables <- build_kk_exact_fisher_tables(des)
	ref <- stats::mantelhaen.test(ref_tables, exact = TRUE, conf.level = 0.95)
	inf_exact <- InferenceIncidExactFisher$new(des, verbose = FALSE)

	ci_exact <- inf_exact$compute_exact_confidence_interval(alpha = 0.05)
	p_exact <- inf_exact$compute_exact_two_sided_pval_for_treatment_effect(delta = 0)

	expect_equal(inf_exact$compute_estimate(), log(as.numeric(ref$estimate)), tolerance = 1e-12)
	expect_equal(unname(ci_exact), log(as.numeric(ref$conf.int)), tolerance = 1e-12)
	expect_equal(p_exact, ref$p.value, tolerance = 1e-12)
})

test_that("Fisher exact inference is rejected for unsupported incidence designs", {
	n <- 12
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.table(x1 = i))
	}
	add_all_subject_responses_seq(des, rep(c(0L, 1L), length.out = n))

	inf_exact <- InferenceIncidExactFisher$new(des, verbose = FALSE)

	expect_error(
		inf_exact$compute_exact_two_sided_pval_for_treatment_effect(),
		"Fisher exact inference requires iBCRD, blocking, or matching designs"
	)
})
