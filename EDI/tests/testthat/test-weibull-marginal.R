test_that("InferenceSurvivalKKWeibullMarginal matches survreg cluster-robust fit in R", {
	n <- 60
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "survival", verbose = FALSE)
	set.seed(5)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	w01 <- as.numeric(des$get_w() == 1)
	y_lat <- rexp(n) * exp(0.7 * w01 + rnorm(n, sd = 0.3))
	cens <- as.numeric(quantile(y_lat, 0.85))
	y <- pmin(y_lat, cens)
	dead <- as.numeric(y_lat <= cens)
	add_all_subject_responses_seq(des, y, deads = dead)

	inf <- InferenceSurvivalKKWeibullMarginal$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	ci <- inf$compute_asymp_confidence_interval()
	se <- as.numeric(ci[2] - est) / qnorm(0.975)
	pv <- inf$compute_asymp_two_sided_pval()

	# reference: identical model fit in R by survreg with the same pair-clusters
	ipriv <- inf$.__enclos_env__$private
	X_cov <- as.matrix(ipriv$get_X())
	dat <- data.frame(.y__ = ipriv$y, .dead__ = ipriv$dead, treatment = ipriv$w, X_cov, check.names = FALSE)
	dat$.cl__ <- ipriv$get_cluster_ids()
	rhs <- paste(c("treatment", colnames(X_cov)), collapse = " + ")
	fit_r <- survival::survreg(
		stats::as.formula(paste0("survival::Surv(.y__, .dead__) ~ ", rhs)),
		data = dat, dist = "weibull", cluster = .cl__, robust = TRUE
	)
	est_r <- as.numeric(stats::coef(fit_r)[["treatment"]])
	se_r <- sqrt(stats::vcov(fit_r)["treatment", "treatment"])

	expect_equal(est, est_r, tolerance = 1e-3)
	expect_equal(se, se_r, tolerance = 1e-3)
	expect_equal(pv, 2 * stats::pnorm(-abs(est / se)), tolerance = 1e-8)
})

test_that("InferenceSurvivalKKWeibullMarginal matches survreg with no covariate adjustment", {
	n <- 40
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "survival", verbose = FALSE)
	set.seed(6)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
	}
	y <- rexp(n) * exp(0.5 * as.numeric(des$get_w() == 1))
	add_all_subject_responses_seq(des, y, deads = rep(1L, n))

	inf <- InferenceSurvivalKKWeibullMarginal$new(des, model_formula = ~ 1, verbose = FALSE)
	est <- inf$compute_estimate()
	ci <- inf$compute_asymp_confidence_interval()
	se <- as.numeric(ci[2] - est) / qnorm(0.975)

	ipriv <- inf$.__enclos_env__$private
	dat <- data.frame(.y__ = ipriv$y, .dead__ = ipriv$dead, treatment = ipriv$w)
	dat$.cl__ <- ipriv$get_cluster_ids()
	fit_r <- survival::survreg(
		survival::Surv(.y__, .dead__) ~ treatment,
		data = dat, dist = "weibull", cluster = .cl__, robust = TRUE
	)
	expect_equal(est, as.numeric(stats::coef(fit_r)[["treatment"]]), tolerance = 1e-3)
	expect_equal(se, sqrt(stats::vcov(fit_r)["treatment", "treatment"]), tolerance = 1e-3)
})

test_that("InferenceSurvivalKKWeibullMarginal randomization inference path returns valid pval", {
	n <- 40
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "survival", verbose = FALSE)
	set.seed(7)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
	}
	y <- rexp(n) * exp(0.8 * as.numeric(des$get_w() == 1))
	add_all_subject_responses_seq(des, y, deads = rep(1L, n))

	inf <- InferenceSurvivalKKWeibullMarginal$new(des, verbose = FALSE)
	inf$set_seed(1)
	pv <- inf$compute_rand_two_sided_pval(r = 101, show_progress = FALSE)
	expect_true(is.finite(pv))
	expect_true(pv >= 0 && pv <= 1)
})
