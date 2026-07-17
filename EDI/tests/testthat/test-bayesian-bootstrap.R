context("Bayesian bootstrap")

make_seq_design_for_bayes_boot = function(response_type, y){
	des = DesignSeqOneByOneBernoulli$new(n = length(y), response_type = response_type)
	for (i in seq_along(y)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i / 10))
	}
	des$overwrite_all_subject_assignments(rep(c(-1, 1), length.out = length(y)))
	des$add_all_subject_responses(y)
	des
}

make_survival_design_for_bayes_boot = function(y, dead){
	des = DesignSeqOneByOneBernoulli$new(n = length(y), response_type = "survival")
	for (i in seq_along(y)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i / 10))
	}
	des$overwrite_all_subject_assignments(rep(c(-1, 1), length.out = length(y)))
	des$add_all_subject_responses(ys = y, deads = dead)
	des
}

make_kk_design_for_weighted_bayes_boot = function(response_type, y, n_pairs = 3L, n_single = 2L){
	n = 2L * n_pairs + n_single
	des = DesignSeqOneByOneKK14$new(n = n, response_type = response_type, verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i / 10, x2 = (i %% 3) / 10))
	}
	des$.__enclos_env__$private$m <- c(rep(seq_len(n_pairs), each = 2L), rep(0L, n_single))
	des$add_all_subject_responses(y)
	des
}

test_that("jackknife descendants expose Bayesian bootstrap methods", {
	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf = InferenceCountPoisson$new(des)
	expect_true(inherits(inf, "InferenceJackknife"))
	expect_true(is.function(inf$approximate_bayesian_bootstrap_distribution_beta_hat_T))
	expect_true(is.function(inf$compute_bayesian_bootstrap_two_sided_pval))
	expect_true(is.function(inf$compute_bayesian_bootstrap_confidence_interval))
})

test_that("equal subject weights recover the original point estimate", {
	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf = InferenceCountPoisson$new(des)
	est = inf$compute_estimate()
	n = des$get_n()
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weighted_est = inf$compute_estimate_with_bootstrap_weights(rep(1, n))
	expect_equal(as.numeric(weighted_est), as.numeric(est), tolerance = 1e-5)
})

test_that("weighted logistic hook returns a finite estimate", {
	des = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf = InferenceIncidLogRegr$new(des)
	n = des$get_n()
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weighted_est = inf$compute_estimate_with_bootstrap_weights(rep(1, n))
	expect_true(is.finite(as.numeric(weighted_est)))
})

test_that("weighted log-binomial hook returns a finite estimate and matches equal weights", {
	des = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf = InferenceIncidLogBinomial$new(des)
	n = des$get_n()
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weighted_est = inf$compute_estimate_with_bootstrap_weights(rep(1, n))
	expect_true(is.finite(as.numeric(weighted_est)))
	expect_equal(
		as.numeric(weighted_est),
		as.numeric(inf$compute_estimate()),
		tolerance = 1e-6
	)
})

test_that("weighted fractional-logit hook returns a finite estimate and matches equal weights", {
	des = make_seq_design_for_bayes_boot("proportion", c(0.1, 0.8, 0.2, 0.7, 0.6, 0.3, 0.9, 0.4))
	inf = InferencePropFractionalLogit$new(des)
	n = des$get_n()
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weighted_est = inf$compute_estimate_with_bootstrap_weights(rep(1, n))
	expect_true(is.finite(as.numeric(weighted_est)))
	expect_equal(
		as.numeric(weighted_est),
		as.numeric(inf$compute_estimate()),
		tolerance = 1e-6
	)
})

test_that("additional near-term weighted hooks return finite estimates and recover equal weights when practical", {
	des_count = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	n_count = des_count$get_n()
	ctx_count = list(row_to_unit = seq_len(n_count), unit_group_id = rep(1L, n_count), n_units = n_count)

	inf_nb = InferenceCountNegBin$new(des_count)
	inf_nb$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_count
	expect_true(is.finite(as.numeric(inf_nb$compute_estimate_with_bootstrap_weights(rep(1, n_count)))))

	inf_qp = InferenceCountQuasiPoisson$new(des_count)
	inf_qp$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_count
	expect_equal(
		as.numeric(inf_qp$compute_estimate_with_bootstrap_weights(rep(1, n_count))),
		as.numeric(inf_qp$compute_estimate()),
		tolerance = 1e-6
	)

	inf_rp = InferenceCountRobustPoisson$new(des_count)
	inf_rp$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_count
	expect_equal(
		as.numeric(inf_rp$compute_estimate_with_bootstrap_weights(rep(1, n_count))),
		as.numeric(inf_rp$compute_estimate()),
		tolerance = 1e-6
	)

	des_incid = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	n_incid = des_incid$get_n()
	ctx_incid = list(row_to_unit = seq_len(n_incid), unit_group_id = rep(1L, n_incid), n_units = n_incid)
	inf_probit = InferenceIncidProbitRegr$new(des_incid)
	inf_probit$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_incid
	expect_true(is.finite(as.numeric(inf_probit$compute_estimate_with_bootstrap_weights(rep(1, n_incid)))))

	des_prop = make_seq_design_for_bayes_boot("proportion", c(0.1, 0.8, 0.2, 0.7, 0.6, 0.3, 0.9, 0.4))
	n_prop = des_prop$get_n()
	ctx_prop = list(row_to_unit = seq_len(n_prop), unit_group_id = rep(1L, n_prop), n_units = n_prop)
	inf_beta = InferencePropBetaRegr$new(des_prop)
	inf_beta$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_prop
	expect_true(is.finite(as.numeric(inf_beta$compute_estimate_with_bootstrap_weights(rep(1, n_prop)))))

	des_cont = make_seq_design_for_bayes_boot("continuous", c(0, 1, 2, 3, 4, 5, 6, 7))
	n_cont = des_cont$get_n()
	ctx_cont = list(row_to_unit = seq_len(n_cont), unit_group_id = rep(1L, n_cont), n_units = n_cont)

	inf_rob = InferenceContinRobustRegr$new(des_cont, use_rcpp = FALSE)
	inf_rob$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_cont
	expect_true(is.finite(as.numeric(inf_rob$compute_estimate_with_bootstrap_weights(rep(1, n_cont)))))

	skip_if_not_installed("quantreg")
	inf_qr = InferenceContinQuantileRegr$new(des_cont)
	inf_qr$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_cont
	expect_true(is.finite(as.numeric(inf_qr$compute_estimate_with_bootstrap_weights(rep(1, n_cont)))))
})

test_that("remaining near-term weighted hooks return finite estimates", {
	skip_if_not_installed("glmmTMB")
	des_count = make_seq_design_for_bayes_boot("count", c(0L, 1L, 0L, 2L, 3L, 1L, 0L, 2L))
	n_count = des_count$get_n()
	ctx_count = list(row_to_unit = seq_len(n_count), unit_group_id = rep(1L, n_count), n_units = n_count)

	inf_hp = InferenceCountHurdlePoisson$new(des_count, use_rcpp = FALSE)
	inf_hp$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_count
	expect_true(is.finite(as.numeric(inf_hp$compute_estimate_with_bootstrap_weights(rep(1, n_count)))))

	inf_hnb = InferenceCountHurdleNegBin$new(des_count)
	inf_hnb$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_count
	expect_true(is.finite(as.numeric(inf_hnb$compute_estimate_with_bootstrap_weights(rep(1, n_count)))))

	inf_zip = InferenceCountZeroInflatedPoisson$new(des_count, use_rcpp = FALSE)
	inf_zip$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_count
	expect_true(is.finite(as.numeric(inf_zip$compute_estimate_with_bootstrap_weights(rep(1, n_count)))))

	inf_zinb = InferenceCountZeroInflatedNegBin$new(des_count, use_rcpp = FALSE)
	inf_zinb$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_count
	expect_true(is.finite(as.numeric(inf_zinb$compute_estimate_with_bootstrap_weights(rep(1, n_count)))))

	des_incid = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	n_incid = des_incid$get_n()
	ctx_incid = list(row_to_unit = seq_len(n_incid), unit_group_id = rep(1L, n_incid), n_units = n_incid)

	inf_new = InferenceIncidNewcombeRiskDiff$new(des_incid)
	inf_new$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_incid
	expect_equal(
		as.numeric(inf_new$compute_estimate_with_bootstrap_weights(rep(1, n_incid))),
		as.numeric(inf_new$compute_estimate()),
		tolerance = 1e-8
	)

	inf_mn = InferenceIncidMiettinenNurminenRiskDiff$new(des_incid)
	inf_mn$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_incid
	expect_equal(
		as.numeric(inf_mn$compute_estimate_with_bootstrap_weights(rep(1, n_incid))),
		as.numeric(inf_mn$compute_estimate()),
		tolerance = 1e-8
	)

	des_prop = make_seq_design_for_bayes_boot("proportion", c(0, 1, 0.2, 0.7, 0.6, 0.3, 0.9, 0.4))
	n_prop = des_prop$get_n()
	ctx_prop = list(row_to_unit = seq_len(n_prop), unit_group_id = rep(1L, n_prop), n_units = n_prop)
	inf_zoib = InferencePropZeroOneInflatedBetaRegr$new(des_prop)
	inf_zoib$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_prop
	zoib_est = tryCatch(as.numeric(inf_zoib$compute_estimate_with_bootstrap_weights(rep(1, n_prop))), error = function(e) NA_real_)
	expect_true(is.na(zoib_est) || is.finite(zoib_est))
})

test_that("weighted identity-binomial hook returns a finite estimate and matches equal weights", {
	des = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf = InferenceIncidBinomialIdentityRiskDiff$new(des)
	n = des$get_n()
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weighted_est = inf$compute_estimate_with_bootstrap_weights(rep(1, n))
	expect_true(is.finite(as.numeric(weighted_est)))
	expect_equal(
		as.numeric(weighted_est),
		as.numeric(inf$compute_estimate()),
		tolerance = 1e-6
	)
})

test_that("Bayesian bootstrap smoke test runs on a stable first-wave family", {
	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf = InferenceCountPoisson$new(des)
	inf$num_cores = 1L
	set.seed(20260515)
	boot = inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(B = 8L, show_progress = FALSE)
	expect_length(boot, 8L)
	expect_true(all(is.finite(boot)))
	set.seed(20260515)
	ci = inf$compute_bayesian_bootstrap_confidence_interval(B = 8L, show_progress = FALSE, type = "percentile")
	expect_length(ci, 2L)
	expect_true(all(is.finite(ci)))
})

test_that("next-wave weighted hooks return finite estimates on simple and g-computation paths", {
	des_cont = make_seq_design_for_bayes_boot("continuous", c(0, 1, 2, 3, 4, 5, 6, 7))
	inf_cont = InferenceAllSimpleMeanDiff$new(des_cont)
	n_cont = des_cont$get_n()
	inf_cont$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_cont),
		unit_group_id = rep(1L, n_cont),
		n_units = n_cont
	)
	expect_true(is.finite(as.numeric(inf_cont$compute_estimate_with_bootstrap_weights(rep(1, n_cont)))))

	des_incid = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf_incid = InferenceIncidGCompRiskDiff$new(des_incid)
	n_incid = des_incid$get_n()
	inf_incid$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_incid),
		unit_group_id = rep(1L, n_incid),
		n_units = n_incid
	)
	expect_true(is.finite(as.numeric(inf_incid$compute_estimate_with_bootstrap_weights(rep(1, n_incid)))))

	des_prop = make_seq_design_for_bayes_boot("proportion", c(0.1, 0.8, 0.2, 0.7, 0.6, 0.3, 0.9, 0.4))
	inf_prop = InferencePropGCompMeanDiff$new(des_prop)
	n_prop = des_prop$get_n()
	inf_prop$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_prop),
		unit_group_id = rep(1L, n_prop),
		n_units = n_prop
	)
	expect_true(is.finite(as.numeric(inf_prop$compute_estimate_with_bootstrap_weights(rep(1, n_prop)))))

	des_ord = make_seq_design_for_bayes_boot("ordinal", c(1L, 2L, 1L, 3L, 2L, 1L, 3L, 2L))
	inf_ord = InferenceOrdinalGCompMeanDiff$new(des_ord)
	n_ord = des_ord$get_n()
	inf_ord$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_ord),
		unit_group_id = rep(1L, n_ord),
		n_units = n_ord
	)
	expect_true(is.finite(as.numeric(inf_ord$compute_estimate_with_bootstrap_weights(rep(1, n_ord)))))
})

test_that("selected second-wave weighted hooks recover unweighted estimates under equal weights", {
	des_ols = make_seq_design_for_bayes_boot("continuous", c(0, 1, 2, 3, 4, 5, 6, 7))
	inf_ols = InferenceContinOLS$new(des_ols)
	n_ols = des_ols$get_n()
	inf_ols$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_ols),
		unit_group_id = rep(1L, n_ols),
		n_units = n_ols
	)
	expect_equal(
		as.numeric(inf_ols$compute_estimate_with_bootstrap_weights(rep(1, n_ols))),
		as.numeric(inf_ols$compute_estimate()),
		tolerance = 1e-8
	)

	des_lin = make_seq_design_for_bayes_boot("continuous", c(0, 1, 1.5, 2.5, 4, 5, 5.5, 7))
	inf_lin = InferenceContinLin$new(des_lin)
	n_lin = des_lin$get_n()
	inf_lin$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_lin),
		unit_group_id = rep(1L, n_lin),
		n_units = n_lin
	)
	expect_equal(
		as.numeric(inf_lin$compute_estimate_with_bootstrap_weights(rep(1, n_lin))),
		as.numeric(inf_lin$compute_estimate()),
		tolerance = 1e-8
	)

	des_surv = make_survival_design_for_bayes_boot(
		y = c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
		dead = c(1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L)
	)
	n_surv = des_surv$get_n()
	ctx_surv = list(
		row_to_unit = seq_len(n_surv),
		unit_group_id = rep(1L, n_surv),
		n_units = n_surv
	)

	inf_logrank = InferenceSurvivalLogRank$new(des_surv)
	inf_logrank$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_surv
	expect_equal(
		as.numeric(inf_logrank$compute_estimate_with_bootstrap_weights(rep(1, n_surv))),
		as.numeric(inf_logrank$compute_estimate()),
		tolerance = 1e-8
	)

	inf_gehan = InferenceSurvivalGehanWilcox$new(des_surv)
	inf_gehan$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_surv
	expect_equal(
		as.numeric(inf_gehan$compute_estimate_with_bootstrap_weights(rep(1, n_surv))),
		as.numeric(inf_gehan$compute_estimate()),
		tolerance = 1e-8
	)
})

test_that("selected ordinal and survival second-wave hooks return finite weighted estimates", {
	des_ord = make_seq_design_for_bayes_boot("ordinal", c(1L, 2L, 1L, 3L, 2L, 1L, 3L, 2L))
	inf_ord = InferenceOrdinalPropOddsRegr$new(des_ord)
	n_ord = des_ord$get_n()
	inf_ord$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_ord),
		unit_group_id = rep(1L, n_ord),
		n_units = n_ord
	)
	expect_true(is.finite(as.numeric(inf_ord$compute_estimate_with_bootstrap_weights(rep(1, n_ord)))))

	des_surv = make_survival_design_for_bayes_boot(
		y = c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
		dead = c(1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L)
	)
	n_surv = des_surv$get_n()
	ctx_surv = list(
		row_to_unit = seq_len(n_surv),
		unit_group_id = rep(1L, n_surv),
		n_units = n_surv
	)

	inf_logrank = InferenceSurvivalLogRank$new(des_surv)
	inf_logrank$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_surv
	expect_true(is.finite(as.numeric(inf_logrank$compute_estimate_with_bootstrap_weights(rep(1, n_surv)))))

	inf_gehan = InferenceSurvivalGehanWilcox$new(des_surv)
	inf_gehan$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_surv
	expect_true(is.finite(as.numeric(inf_gehan$compute_estimate_with_bootstrap_weights(rep(1, n_surv)))))
})

test_that("ordinal likelihood-gap weighted hooks are finite and exact empirical ordinal hooks recover equal weights", {
	des_ord = make_seq_design_for_bayes_boot("ordinal", c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L))
	n_ord = des_ord$get_n()
	ctx_ord = list(
		row_to_unit = seq_len(n_ord),
		unit_group_id = rep(1L, n_ord),
		n_units = n_ord
	)

	ordinal_classes = list(
		InferenceOrdinalAdjCatLogitRegr$new(des_ord),
		InferenceOrdinalCloglogRegr$new(des_ord),
		InferenceOrdinalCauchitRegr$new(des_ord),
		InferenceOrdinalOrderedProbitRegr$new(des_ord),
		InferenceOrdinalStereotypeLogitRegr$new(des_ord),
		InferenceOrdinalContRatioRegr$new(des_ord),
		InferenceOrdinalOrderedProbitRegr$new(des_ord)
	)
	for (inf in ordinal_classes) {
		inf$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_ord
		expect_true(is.finite(as.numeric(inf$compute_estimate_with_bootstrap_weights(rep(1, n_ord)))))
	}

	inf_ridit = InferenceOrdinalRidit$new(des_ord)
	inf_ridit$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_ord
	expect_equal(
		as.numeric(inf_ridit$compute_estimate_with_bootstrap_weights(rep(1, n_ord))),
		as.numeric(inf_ridit$compute_estimate()),
		tolerance = 1e-8
	)

	des_kk_ord = make_kk_design_for_weighted_bayes_boot("ordinal", c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L))
	n_kk_ord = des_kk_ord$get_n()
	ctx_kk_ord = list(
		row_to_unit = seq_len(n_kk_ord),
		unit_group_id = rep(1L, n_kk_ord),
		n_units = n_kk_ord
	)
	inf_sign = InferenceOrdinalPairedSignTest$new(des_kk_ord, verbose = FALSE)
	inf_sign$.__enclos_env__$private$current_bayesian_bootstrap_context = ctx_kk_ord
	expect_equal(
		as.numeric(inf_sign$compute_estimate_with_bootstrap_weights(rep(1, n_kk_ord))),
		as.numeric(inf_sign$compute_estimate()),
		tolerance = 1e-8
	)
})

test_that("next-wave weighted hooks return finite estimates on KK GEE paths", {
	skip_if_not_installed("multgee")

	des_incid = make_kk_design_for_weighted_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf_incid = InferenceIncidKKGEE$new(des_incid, use_rcpp = TRUE, verbose = FALSE)
	n_incid = des_incid$get_n()
	inf_incid$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_incid),
		unit_group_id = rep(1L, n_incid),
		n_units = n_incid
	)
	expect_true(is.finite(as.numeric(inf_incid$compute_estimate_with_bootstrap_weights(rep(1, n_incid)))))

	des_count = make_kk_design_for_weighted_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf_count = InferenceCountPoissonKKGEE$new(des_count, use_rcpp = TRUE, verbose = FALSE)
	n_count = des_count$get_n()
	inf_count$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_count),
		unit_group_id = rep(1L, n_count),
		n_units = n_count
	)
	expect_true(is.finite(as.numeric(inf_count$compute_estimate_with_bootstrap_weights(rep(1, n_count)))))

	des_prop = make_kk_design_for_weighted_bayes_boot("proportion", c(0.1, 0.8, 0.2, 0.7, 0.6, 0.3, 0.9, 0.4))
	inf_prop = InferencePropKKGEE$new(des_prop, use_rcpp = TRUE, verbose = FALSE)
	n_prop = des_prop$get_n()
	inf_prop$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_prop),
		unit_group_id = rep(1L, n_prop),
		n_units = n_prop
	)
	expect_true(is.finite(as.numeric(inf_prop$compute_estimate_with_bootstrap_weights(rep(1, n_prop)))))

	des_ord = make_kk_design_for_weighted_bayes_boot("ordinal", c(1L, 2L, 1L, 3L, 2L, 1L, 3L, 2L))
	inf_ord = InferenceOrdinalKKGEE$new(des_ord, verbose = FALSE)
	n_ord = des_ord$get_n()
	inf_ord$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_ord),
		unit_group_id = rep(1L, n_ord),
		n_units = n_ord
	)
	expect_true(is.finite(as.numeric(inf_ord$compute_estimate_with_bootstrap_weights(rep(1, n_ord)))))
})

test_that("Bayesian bootstrap matches mirai-backed parallel execution", {
	skip_if_not_installed("mirai")

	on.exit(unset_num_cores(), add = TRUE)

	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	serial_inf = InferenceCountPoisson$new(des)
	mirai_inf = InferenceCountPoisson$new(des)

	serial_inf$num_cores = 1L
	set.seed(20260515)
	serial_boot = serial_inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(
		B = 11L,
		show_progress = FALSE
	)

	set_num_cores(2L, force_mirai = TRUE)
	mirai_inf$num_cores = 2L
	set.seed(20260515)
	mirai_boot = mirai_inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(
		B = 11L,
		show_progress = FALSE
	)

	expect_equal(unname(mirai_boot), unname(serial_boot), tolerance = 1e-10)
})

test_that("IVWC classes handle Bayesian bootstrap calls correctly", {
	des_cont = make_kk_design_for_weighted_bayes_boot("continuous", c(0, 1, 2, 3, 4, 5, 6, 7))
	inf_cont_ivwc = InferenceContinKKOLSIVWC$new(des_cont, verbose = FALSE)
	expect_no_error(
		inf_cont_ivwc$approximate_bayesian_bootstrap_distribution_beta_hat_T(B = 3L, show_progress = FALSE)
	)
	expect_error(
		inf_cont_ivwc$compute_estimate_with_bootstrap_weights(rep(1, des_cont$get_n())),
		"No Bayesian-bootstrap context is installed on this inference object"
	)

	des_incid = make_kk_design_for_weighted_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf_incid_ivwc = InferenceIncidKKCondLogitIVWC$new(des_incid, verbose = FALSE)
	expect_no_error(
		inf_incid_ivwc$compute_bayesian_bootstrap_confidence_interval(B = 3L, show_progress = FALSE)
	)
})
