test_that("Inference works for continuous", {
	n <- 10
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rnorm(n))

	# Simple Mean Diff
	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	expect_true(is.numeric(est))

	# OLS
	inf_ols <- InferenceContinOLS$new(des, verbose = FALSE)
	est_ols <- inf_ols$compute_estimate()
	expect_true(is.numeric(est_ols))
})

test_that("Inference works for incidence", {
	n <- 30
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rbinom(n, 1, 0.5))

	inf <- InferenceIncidLogRegr$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	expect_true(is.numeric(est))
})

test_that("Simple incidence proportion difference uses pooled-variance t inference", {
	des <- EDI:::DesignFixed$new(n = 10, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = 1:10))
	des$overwrite_all_subject_assignments(c(1, 1, 1, 1, 1, -1, -1, -1, -1, -1))
	des$add_all_subject_responses(c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0))

	inf <- InferenceAllSimpleMeanDiffPooledVar$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	ci <- inf$compute_asymp_confidence_interval(alpha = 0.05)
	pval <- inf$compute_asymp_two_sided_pval()

	y_t <- c(1, 1, 0, 1, 0)
	y_c <- c(0, 1, 0, 0, 0)
	p_t <- mean(y_t)
	p_c <- mean(y_c)
	s2_t <- stats::var(y_t)
	s2_c <- stats::var(y_c)
	df <- length(y_t) + length(y_c) - 2
	s2_pooled <- ((length(y_t) - 1) * s2_t + (length(y_c) - 1) * s2_c) / df
	se <- sqrt(s2_pooled * (1 / length(y_t) + 1 / length(y_c)))
	crit <- stats::qt(0.975, df = df)
	expected_est <- p_t - p_c
	expected_ci <- c(expected_est - crit * se, expected_est + crit * se)
	expected_pval <- 2 * stats::pt(-abs(expected_est / se), df = df)

	expect_equal(est, expected_est, tolerance = 1e-12)
	expect_equal(unname(ci), expected_ci, tolerance = 1e-12)
	expect_equal(pval, expected_pval, tolerance = 1e-12)

	# Test with continuous
	des_cont <- DesignSeqOneByOneBernoulli$new(n = 10, response_type = "continuous", verbose = FALSE)
	set.seed(1)
	for (i in 1:10) {
		des_cont$add_one_subject_to_experiment_and_assign(data.frame(x = i))
	}
	y_cont = rnorm(10)
	add_all_subject_responses_seq(des_cont, y_cont)
	
	inf_cont <- InferenceAllSimpleMeanDiffPooledVar$new(des_cont, verbose = FALSE)
	est_cont <- inf_cont$compute_estimate()
	expect_true(is.numeric(est_cont))
	expect_true(is.finite(inf_cont$compute_asymp_two_sided_pval()))
})

test_that("Simple mean difference pooled-variance inference matches pooled t.test", {
	des <- EDI:::DesignFixed$new(n = 9, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = 1:9))
	w <- c(1, 1, 1, 1, -1, -1, -1, -1, -1)
	y <- c(2.2, 2.6, 1.8, 3.1, 0.9, 1.0, 1.4, 0.8, 1.2)
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)

	inf <- InferenceAllSimpleMeanDiffPooledVar$new(des, verbose = FALSE)
	tt <- stats::t.test(y[w == 1], y[w == -1], var.equal = TRUE)

	expect_equal(inf$compute_asymp_two_sided_pval(), tt$p.value, tolerance = 1e-12)
	expect_equal(as.numeric(inf$compute_asymp_confidence_interval()), as.numeric(tt$conf.int), tolerance = 1e-12)
})

test_that("CMH inference is gated to blocked incidence designs", {
	des <- DesignFixedBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c("A", "A", "A", "A", "B", "B", "B", "B")))
	des$overwrite_all_subject_assignments(c(1, -1, 1, -1, 1, -1, 1, -1))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	inf <- InferenceIncidCMH$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	ci <- inf$compute_asymp_confidence_interval()
	pval <- inf$compute_asymp_two_sided_pval()

	y <- c(1, 0, 1, 0, 1, 1, 0, 0)
	w <- c(1, 0, 1, 0, 1, 0, 1, 0)
	m <- c(1, 1, 1, 1, 2, 2, 2, 2)
	expected_est <- mean(y[w == 1]) - mean(y[w == 0])
	var_cmh <- 0
	for (m_j in unique(m)) {
		i_m <- m == m_j
		y_m <- y[i_m]
		n_m <- sum(i_m)
		Sigma_m <- matrix(-1 / (n_m - 1), nrow = n_m, ncol = n_m)
		diag(Sigma_m) <- 1
		var_cmh <- var_cmh + as.numeric(t(y_m) %*% Sigma_m %*% y_m)
	}
	expected_ci <- c(-0.300151946059218, 1.300151946059218)
	expected_pval <- 0.2206713619198468

	expect_equal(est, expected_est, tolerance = 1e-12)
	expect_equal(unname(ci), expected_ci, tolerance = 1e-12)
	expect_equal(as.numeric(pval), expected_pval, tolerance = 1e-8)

	des_nonblock <- DesignSeqOneByOneBernoulli$new(n = 8, response_type = "incidence", verbose = FALSE)
	for (i in 1:8) {
		des_nonblock$add_one_subject_to_experiment_and_assign(data.frame(x = i))
	}
	add_all_subject_responses_seq(des_nonblock, c(1, 0, 1, 0, 1, 0, 1, 0))
	inf_nonblock <- InferenceIncidCMH$new(des_nonblock, verbose = FALSE)
	expect_true(is.finite(inf_nonblock$compute_estimate()))
})

test_that("CMH inference requires even treatment allocation", {
	des <- DesignFixedBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		prob_T = 0.25,
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c("A", "A", "A", "A", "B", "B", "B", "B")))
	des$overwrite_all_subject_assignments(c(1, -1, 1, -1, 1, -1, 1, -1))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	expect_error(
		InferenceIncidCMH$new(des, verbose = FALSE),
		"even treatment allocation"
	)
})

test_that("CMH inference requires equal block sizes", {
	des <- DesignFixedBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c(rep("A", 3), rep("B", 5))))
	des$overwrite_all_subject_assignments(c(1, -1, 1, -1, 1, -1, 1, -1))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	expect_error(
		InferenceIncidCMH$new(des, verbose = FALSE),
		"equal_block_sizes = TRUE"
	)
})

test_that("CMH and Extended Robins standard errors match a fixed blocked simulation", {
	set.seed(20260405)
	expected <- matrix(
		c(
			0.3818813079129866, 0.3644344934278313,
			0.3535533905932738, 0.3061862178478972,
			0.2886751345948129, 0.3307189138830738,
			0.2500000000000000, 0.2651650429449553,
			0.3818813079129866, 0.3644344934278313
		),
		ncol = 2, byrow = TRUE,
		dimnames = list(
			NULL,
			c("cmh", "robins")
		)
	)

	ys <- replicate(5, rbinom(8, 1, 0.5), simplify = FALSE)
	se_pairs <- t(vapply(ys, function(y) {
		des <- DesignFixedBlocking$new(
			strata_cols = "stratum",
			n = 8,
			response_type = "incidence",
			verbose = FALSE
		)
		des$add_all_subjects_to_experiment(data.frame(stratum = c(rep("A", 4), rep("B", 4))))
		des$overwrite_all_subject_assignments(c(1, -1, 1, -1, 1, -1, 1, -1))
		des$add_all_subject_responses(y)

		inf_cmh <- InferenceIncidCMH$new(des, verbose = FALSE)
		inf_robins <- InferenceIncidExtendedRobins$new(des, verbose = FALSE)
		c(
			cmh = inf_cmh$.__enclos_env__$private$get_standard_error(),
			robins = inf_robins$.__enclos_env__$private$get_standard_error()
		)
	}, numeric(2)))

	expect_equal(se_pairs, expected, tolerance = 1e-12)
})

test_that("CMH get_standard_error block and non-block paths agree when D = 3 and R = n - 1", {
	# Mathematical basis: when y has exactly 3 ones (3 discordant pairs, rest all-zero),
	# E[(y^T w_r)^2] = 3 for BOTH paired and Bernoulli randomization, and with
	# se_est_num_vectors = n - 1 the formula reduces to 2/n * sqrt(3) = SE_block.
	set.seed(42)
	n <- 50L
	n_pairs <- n / 2L
	y <- c(1, 0, 1, 0, 1, 0, rep(0L, n - 6L))

	# block path: paired blocking design, 3 discordant + 22 concordant (0,0) pairs
	des_block <- DesignFixedBlocking$new(
		strata_cols = "s",
		n = n,
		B_target = n_pairs,
		response_type = "incidence",
		verbose = FALSE
	)
	des_block$add_all_subjects_to_experiment(
		data.frame(s = rep(as.character(seq_len(n_pairs)), each = 2L))
	)
	des_block$overwrite_all_subject_assignments(rep(c(1L, -1L), n_pairs))
	des_block$add_all_subject_responses(y)
	inf_block <- InferenceIncidCMH$new(des_block, verbose = FALSE)
	se_block <- inf_block$.__enclos_env__$private$get_standard_error()
	expect_equal(se_block, 2 / n * sqrt(3), tolerance = 1e-12)

	# non-block path: Bernoulli design, same y, se_est_num_vectors = n - 1
	des_nonblock <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) {
		des_nonblock$add_one_subject_to_experiment_and_assign(data.frame(x = i))
	}
	add_all_subject_responses_seq(des_nonblock, y)
	inf_nonblock <- InferenceIncidCMH$new(des_nonblock, se_est_num_vectors = n - 1L, verbose = FALSE)
	se_nonblock <- inf_nonblock$.__enclos_env__$private$get_standard_error()
	expect_equal(se_block, se_nonblock, tolerance = 0.15)
})

test_that("CMH get_standard_error block and non-block paths agree for many D with n = 100", {
	# For n = 100, n_pairs = 50 pairs, y with exactly D ones (D discordant pairs),
	# SE_block = 2/100 * sqrt(D) exactly.
	# SE_nonblock with R = 4*(n-1)/(D+1) vectors satisfies E[SE_nonblock^2] = SE_block^2,
	# so both paths agree in expectation for each case below.
	n <- 100L
	n_pairs <- n / 2L
	cases <- list(
		list(D = 1L, R = 198L),
		list(D = 2L, R = 132L),
		list(D = 3L, R = 99L),
		list(D = 5L, R = 66L),
		list(D = 8L, R = 44L)
	)
	for (k in seq_along(cases)) {
		D <- cases[[k]]$D
		R <- cases[[k]]$R
		y <- c(rep(c(1L, 0L), D), rep(0L, n - 2L * D))

		# block path
		des_block <- DesignFixedBlocking$new(
			strata_cols = "s",
			n = n,
			B_target = n_pairs,
			response_type = "incidence",
			verbose = FALSE
		)
		des_block$add_all_subjects_to_experiment(
			data.frame(s = rep(as.character(seq_len(n_pairs)), each = 2L))
		)
		des_block$overwrite_all_subject_assignments(rep(c(1L, -1L), n_pairs))
		des_block$add_all_subject_responses(y)
		inf_block <- InferenceIncidCMH$new(des_block, verbose = FALSE)
		se_block <- inf_block$.__enclos_env__$private$get_standard_error()
		expect_equal(se_block, 2 / n * sqrt(D), tolerance = 1e-12,
			label = sprintf("block SE exact for D=%d", D))

		# non-block path
		set.seed(k * 7L)
		des_nonblock <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
		for (i in seq_len(n)) {
			des_nonblock$add_one_subject_to_experiment_and_assign(data.frame(x = i))
		}
		add_all_subject_responses_seq(des_nonblock, y)
		inf_nonblock <- InferenceIncidCMH$new(des_nonblock, se_est_num_vectors = R, verbose = FALSE)
		se_nonblock <- inf_nonblock$.__enclos_env__$private$get_standard_error()
		expect_equal(se_block, se_nonblock, tolerance = 0.25,
			label = sprintf("non-block SE ≈ block SE for D=%d, R=%d", D, R))
	}
})

test_that("CMH and Extended Robins confidence intervals use normal critical values", {
	des <- DesignFixedBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c(rep("A", 4), rep("B", 4))))
	des$overwrite_all_subject_assignments(c(1, -1, 1, -1, 1, -1, 1, -1))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	for (inf in list(
		InferenceIncidCMH$new(des, verbose = FALSE),
		InferenceIncidExtendedRobins$new(des, verbose = FALSE)
	)) {
		est <- inf$compute_estimate()
		se <- inf$.__enclos_env__$private$get_standard_error()
		ci <- inf$compute_asymp_confidence_interval(alpha = 0.05)
		expect_equal(unname(ci), est + c(-1, 1) * qnorm(0.975) * se, tolerance = 1e-12)
	}
})

test_that("Extended Robins standard error matches the blockwise formula", {
	des <- DesignFixedBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c(rep("A", 4), rep("B", 4))))
	des$overwrite_all_subject_assignments(c(1, -1, 1, -1, 1, -1, 1, -1))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	inf <- InferenceIncidExtendedRobins$new(des, verbose = FALSE)
	se_cpp <- inf$.__enclos_env__$private$get_standard_error()

	m <- des$get_block_ids()
	B <- length(unique(m))
	n_B <- sum(m == 1L)
	n_B_over_two <- n_B / 2
	variance_tot <- 0
	for (b in 1:B) {
		y_b <- des$get_y()[m == b]
		w_b <- des$get_w()[m == b]
		p_hat_T_b <- sum(y_b[w_b == 1]) / n_B_over_two
		p_hat_C_b <- sum(y_b[w_b == -1]) / n_B_over_two
		m_1_b <- max(p_hat_T_b, p_hat_C_b)
		m_0_b <- min(p_hat_T_b, p_hat_C_b)
		variance_tot <- variance_tot +
			m_1_b * (1 - m_1_b) / n_B_over_two +
			m_0_b * (1 - m_0_b) / n_B_over_two +
			((2 * m_0_b - m_1_b) * (1 - m_1_b) - m_0_b * (1 - m_0_b)) / n_B
	}
	p_hat_T <- mean(des$get_y()[des$get_w() == 1])
	p_hat_C <- mean(des$get_y()[des$get_w() == -1])
	var_robbins_ext <- 1 / des$get_n() * (
		p_hat_T * (1 - p_hat_T) + p_hat_C * (1 - p_hat_C)
	)
	se_r <- sqrt(variance_tot / (B * B) + var_robbins_ext)

	expect_equal(se_cpp, se_r, tolerance = 1e-12)
})

test_that("G-computation risk-ratio intervals error when log-scale bounds overflow", {
	des <- DesignFixediBCRD$new(n = 4, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = 1:4))
	des$overwrite_all_subject_assignments(c(1, -1, 1, -1))
	des$add_all_subject_responses(c(1, 0, 1, 0))

	inf <- InferenceIncidGCompRiskRatio$new(des, verbose = FALSE)
	priv <- inf$.__enclos_env__$private
	priv$cached_values$rr <- 1
	priv$cached_values$s_beta_hat_T <- 1
	priv$cached_values$log_rr <- 1000
	priv$cached_values$se_log_rr <- 100

	expect_error(
		inf$compute_asymp_confidence_interval(),
		"could not compute a finite delta-method confidence interval"
	)
})

test_that("Inference works for count", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rpois(n, 5))

	inf <- InferenceCountNegBin$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	expect_true(is.numeric(est))

	inf_robust_pois <- InferenceCountRobustPoisson$new(des, verbose = FALSE)
	est_robust_pois <- inf_robust_pois$compute_estimate()
	expect_true(is.numeric(est_robust_pois))

	inf_quasi_pois <- InferenceCountQuasiPoisson$new(des, verbose = FALSE)
	est_quasi_pois <- inf_quasi_pois$compute_estimate()
	expect_true(is.numeric(est_quasi_pois))

	if (requireNamespace("glmmTMB", quietly = TRUE)) {
		inf_zinb <- InferenceCountZeroInflatedNegBin$new(des, verbose = FALSE)
		est_zinb <- inf_zinb$compute_estimate()
		expect_true(is.numeric(est_zinb))
	}

	inf_hurdle_nb <- InferenceCountHurdleNegBin$new(des, verbose = FALSE)
	est_hurdle_nb <- inf_hurdle_nb$compute_estimate()
	expect_true(is.numeric(est_hurdle_nb))

	est_hurdle_nb_fast <- inf_hurdle_nb$compute_estimate(estimate_only = TRUE)
	expect_true(is.numeric(est_hurdle_nb_fast))
	expect_equal(est_hurdle_nb_fast, est_hurdle_nb, tolerance = 1e-8)
})

test_that("KK count combined-likelihood multi inference handles full-width covariates", {
	skip_if_not_installed("MASS")
	skip_if_not_installed("dplyr")

	set.seed(1)
	cars_subset <- MASS::Cars93 |>
		stats::na.omit() |>
		dplyr::slice_sample(n = 48, replace = TRUE) |>
		dplyr::select(-Make, -Model)
	X <- model.matrix(Price ~ 0 + ., data = cars_subset)
	X <- apply(X, 2, scale)
	X <- X / ncol(X)
	X <- as.data.frame(X)
	X <- dplyr::select(X, dplyr::where(~ !any(is.na(.))))
	y <- round(cars_subset$Price - min(cars_subset$Price))

	des <- DesignSeqOneByOneKK14$new(n = nrow(X), response_type = "count", verbose = FALSE)
	for (i in seq_len(nrow(X))) {
		w_i <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		y_i <- round(y[i] * exp(rnorm(1, mean = 0, sd = 0.1) * w_i))
		des$add_one_subject_response(i, y_i, 1)
	}

	inf <- InferenceCountKKCondPoissonOneLik$new(des, verbose = FALSE)
	est <- inf$compute_estimate()

	expect_true(is.numeric(est))
	expect_length(est, 1L)
	expect_true(is.finite(est))
})

test_that("Inference works for proportion", {
	n <- 10
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "proportion", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, runif(n))

	inf <- InferencePropBetaRegr$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	expect_true(is.numeric(est))
})

test_that("Inference works for survival", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "survival")
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rexp(n), deads = rbinom(n, 1, 0.8))

	# KM Diff
	inf <- InferenceSurvivalKMDiff$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	expect_true(is.numeric(est))

	# Log-rank
	inf_logrank <- InferenceSurvivalLogRank$new(des, verbose = FALSE)
	est_logrank <- inf_logrank$compute_estimate()
	expect_true(is.numeric(est_logrank))

	# Cox PH
	inf_cox <- InferenceSurvivalCoxPHRegr$new(des, verbose = FALSE)
	est_cox <- inf_cox$compute_estimate()
	expect_true(is.numeric(est_cox))

	# Stratified Cox PH
	inf_strat_cox <- InferenceSurvivalStratCoxPHRegr$new(des, verbose = FALSE)
	est_strat_cox <- inf_strat_cox$compute_estimate()
	expect_true(is.numeric(est_strat_cox))
})

test_that("Inference works for ordinal partial proportional odds", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "ordinal")
	set.seed(10)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	y_levels <- sample(1:4, n, replace = TRUE)
	add_all_subject_responses_seq(des, y_levels)

	inf_ppod <- InferenceOrdinalPartialProportionalOddsRegr$new(des, nonparallel = c("x"), verbose = FALSE)
	est_ppod <- inf_ppod$compute_estimate()
	expect_true(is.numeric(est_ppod))
	pval_ppod <- inf_ppod$compute_asymp_two_sided_pval()
	expect_true(is.numeric(pval_ppod))
	expect_false(inherits(inf_ppod, "InferenceOrdinalPartialProportionalOddsAbstract"))
})

test_that("Inference works for incidence KK Newcombe IVWC", {
	n <- 20
	set.seed(1)
	x_dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.5))
	y <- rbinom(n, 1, 0.5)

	seq_des <- DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in 1:n) {
		seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(seq_des, y)

	inf <- InferenceIncidKKNewcombeRiskDiff$new(seq_des, verbose = FALSE)
	est <- inf$compute_estimate()
	expect_true(is.numeric(est))
	
	ci <- inf$compute_asymp_confidence_interval()
	expect_true(is.numeric(ci))
	expect_length(ci, 2)
	
	pval <- inf$compute_asymp_two_sided_pval()
	expect_true(is.numeric(pval))
})

test_that("ordinal hardening drops QR-ranked covariates only when enabled", {
	build_cars93_ordinal_design <- function(design_ctor) {
		dat <- stats::na.omit(MASS::Cars93)
		y_num <- dat$Price
		qs <- stats::quantile(y_num, probs = c(0.25, 0.5, 0.75))
		y_ord <- as.integer(cut(y_num, breaks = c(-Inf, qs, Inf), labels = FALSE))
		x_dat <- subset(dat, select = -Price)
		for (j in seq_len(ncol(x_dat))) {
			if (is.numeric(x_dat[[j]])) {
				x_dat[[j]] <- ifelse(
					x_dat[[j]] <= stats::median(x_dat[[j]], na.rm = TRUE),
					"low",
					"high"
				)
			}
		}

		des <- design_ctor(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
		for (i in seq_len(nrow(x_dat))) {
			des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		}
		add_all_subject_responses_seq(des, y_ord)
		des
	}

	kk_des <- build_cars93_ordinal_design(DesignSeqOneByOneKK14$new)

	adj_hardened <- InferenceOrdinalAdjCatLogitRegr$new(kk_des, verbose = FALSE, harden = TRUE)
	adj_raw <- InferenceOrdinalAdjCatLogitRegr$new(kk_des, verbose = FALSE, harden = FALSE)
	expect_true(is.finite(adj_hardened$compute_estimate()))
	expect_true(is.finite(adj_hardened$compute_asymp_two_sided_pval()))
	expect_true(all(is.finite(adj_hardened$compute_asymp_confidence_interval())))
	expect_true(is.finite(adj_raw$compute_asymp_two_sided_pval()))
	expect_true(all(is.finite(adj_raw$compute_asymp_confidence_interval())))

	kk_adj_hardened <- InferenceOrdinalKKCondAdjCatLogitRegr$new(kk_des, verbose = FALSE, harden = TRUE)
	kk_adj_raw <- InferenceOrdinalKKCondAdjCatLogitRegr$new(kk_des, verbose = FALSE, harden = FALSE)
	expect_true(is.finite(kk_adj_hardened$compute_estimate()))
	expect_true(is.finite(kk_adj_hardened$compute_asymp_two_sided_pval()))
	expect_true(all(is.finite(kk_adj_hardened$compute_asymp_confidence_interval())))
	expect_false(is.finite(kk_adj_raw$compute_estimate()))
	expect_false(is.finite(kk_adj_raw$compute_asymp_two_sided_pval()))
	expect_false(all(is.finite(kk_adj_raw$compute_asymp_confidence_interval())))
})

test_that("bootstrap debug preserves per-iteration records for sequential count designs", {
	dat <- stats::na.omit(MASS::Cars93)
	x_dat <- as.data.frame(subset(dat, select = -Price))
	y_count <- as.integer(round(dat$Price))

	strata_cols_to_use <- names(x_dat)[1:min(2, ncol(x_dat))]
	x_spbr <- x_dat
	for (col in strata_cols_to_use) {
		if (is.numeric(x_spbr[[col]])) {
			med <- stats::median(x_spbr[[col]], na.rm = TRUE)
			x_spbr[[col]] <- factor(ifelse(x_spbr[[col]] <= med, "low", "high"))
		}
	}

	des_spbr <- DesignSeqOneByOneSPBR$new(
		strata_cols = strata_cols_to_use,
		block_size = 4,
		response_type = "count",
		n = nrow(x_spbr)
	)
	for (i in seq_len(nrow(x_spbr))) {
		des_spbr$add_one_subject_to_experiment_and_assign(x_spbr[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(des_spbr, y_count)

	inf_spbr <- InferenceAllSimpleMeanDiff$new(des_spbr, verbose = FALSE)
	debug_spbr <- inf_spbr$approximate_bootstrap_distribution_beta_hat_T(
		B = 12,
		show_progress = FALSE,
		debug = TRUE
	)
	expect_length(debug_spbr$values, 12)
	expect_length(debug_spbr$errors, 12)
	expect_length(debug_spbr$warnings, 12)
	expect_true(all(vapply(debug_spbr$errors, is.character, logical(1))))
	expect_true(all(vapply(debug_spbr$warnings, is.character, logical(1))))

	des_kk21 <- DesignSeqOneByOneKK21$new(response_type = "count", n = nrow(x_dat))
	for (i in seq_len(nrow(x_dat))) {
		des_kk21$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(des_kk21, y_count)

	inf_kk21 <- InferenceCountNegBin$new(des_kk21, verbose = FALSE)
	debug_kk21 <- inf_kk21$approximate_bootstrap_distribution_beta_hat_T(
		B = 8,
		show_progress = FALSE,
		debug = TRUE
	)
	expect_length(debug_kk21$values, 8)
	expect_length(debug_kk21$errors, 8)
	expect_length(debug_kk21$warnings, 8)
})

test_that("proportion g-computation bootstrap worker keeps mutable screening state", {
	dat <- stats::na.omit(MASS::Cars93)
	x_dat <- as.data.frame(subset(dat, select = -Price))
	y_prop <- pmin(0.99, pmax(0.01, dat$Price / max(dat$Price)))

	strata_cols_to_use <- names(x_dat)[1:min(2, ncol(x_dat))]
	x_pocock <- x_dat
	for (col in strata_cols_to_use) {
		if (is.numeric(x_pocock[[col]])) {
			med <- stats::median(x_pocock[[col]], na.rm = TRUE)
			x_pocock[[col]] <- factor(ifelse(x_pocock[[col]] <= med, "low", "high"))
		}
	}

	des_pocock <- DesignSeqOneByOnePocockSimon$new(
		strata_cols = strata_cols_to_use,
		response_type = "proportion",
		n = nrow(x_pocock)
	)
	for (i in seq_len(nrow(x_pocock))) {
		des_pocock$add_one_subject_to_experiment_and_assign(x_pocock[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(des_pocock, y_prop)

	inf_pocock <- InferencePropGCompMeanDiff$new(des_pocock, verbose = FALSE)
	debug_pocock <- inf_pocock$approximate_bootstrap_distribution_beta_hat_T(
		B = 20,
		show_progress = FALSE,
		debug = TRUE
	)
	expect_lt(debug_pocock$prop_illegal_values, 1)
	expect_true(any(is.finite(debug_pocock$values)))
	expect_true(all(is.finite(inf_pocock$compute_bootstrap_confidence_interval(B = 20, show_progress = FALSE))))
	expect_true(is.finite(inf_pocock$compute_bootstrap_two_sided_pval(B = 20)))

	des_urn <- DesignSeqOneByOneUrn$new(
		response_type = "proportion",
		n = nrow(x_dat)
	)
	for (i in seq_len(nrow(x_dat))) {
		des_urn$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(des_urn, y_prop)

	inf_urn <- InferencePropGCompMeanDiff$new(des_urn, verbose = FALSE)
	debug_urn <- inf_urn$approximate_bootstrap_distribution_beta_hat_T(
		B = 12,
		show_progress = FALSE,
		debug = TRUE
	)
	expect_lt(debug_urn$prop_illegal_values, 1)
	expect_true(any(is.finite(debug_urn$values)))
	expect_true(all(is.finite(inf_urn$compute_bootstrap_confidence_interval(B = 12, show_progress = FALSE))))
	expect_true(is.finite(inf_urn$compute_bootstrap_two_sided_pval(B = 12)))
})
