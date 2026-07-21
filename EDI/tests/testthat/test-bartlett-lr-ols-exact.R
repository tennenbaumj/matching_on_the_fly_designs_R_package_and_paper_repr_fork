library(testthat)
library(EDI)

# InferenceContinKKOLSOneLik holds sigma2 fixed at the full model's unbiased
# estimate (RSS_full/(n-p)) rather than re-profiling it under each null fit.
# Under that convention, the package's LR statistic is algebraically identical
# to the classical F(1, n-p) statistic for testing a single coefficient, an
# EXACT finite-sample pivot (not merely a higher-order asymptotic approximation).
# These tests validate our get_bartlett_factor_exact() against base R's lm().

make_kk_ols_inference <- function(seed = 1, n = 40){
	set.seed(seed)
	x1 <- rnorm(n)
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i <- des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i]))
		y_i <- 0.5 + 0.8 * w_i + 0.4 * x1[i] + rnorm(1, sd = 1.3)
		des$add_one_subject_response(i, y_i, 1L)
	}
	InferenceContinKKOLSOneLik$new(des, verbose = FALSE)
}

lm_reference <- function(inf, delta = 0){
	priv <- inf$.__enclos_env__$private
	spec <- priv$get_likelihood_test_spec()
	X <- spec$X
	colnames(X) <- paste0("v", seq_len(ncol(X)))
	dat_full <- data.frame(X)
	dat_full$y <- as.numeric(spec$y)
	form_full <- as.formula(paste("y ~ 0 +", paste(colnames(X), collapse = " + ")))
	lm_full <- lm(form_full, data = dat_full)
	vname <- colnames(X)[spec$j]
	est <- coef(lm_full)[[vname]]
	se <- summary(lm_full)$coefficients[vname, "Std. Error"]
	df_resid <- lm_full$df.residual
	tstat <- (est - delta) / se
	pval <- 2 * pt(-abs(tstat), df_resid)
	list(pval = pval, est = est, se = se, df_resid = df_resid)
}

test_that("InferenceContinKKOLSOneLik opts in to exact Bartlett and advertises the testing type", {
	inf <- make_kk_ols_inference(seed = 101)
	priv <- inf$.__enclos_env__$private
	expect_true(isTRUE(priv$supports_bartlett_likelihood_ratio_exact()))
	expect_true("lik_ratio_bartlett_exact" %in% inf$get_supported_testing_types())
})

test_that("Exact Bartlett p-value matches base R lm()'s classical t-test p-value exactly", {
	for (seed in c(101, 202, 303)) {
		inf <- make_kk_ols_inference(seed = seed)
		ref <- lm_reference(inf, delta = 0)
		bartlett_pval <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0)
		expect_true(is.finite(bartlett_pval))
		expect_equal(bartlett_pval, ref$pval, tolerance = 1e-6)
	}
})

test_that("Exact Bartlett p-value differs meaningfully from the raw chi-square-approximated LR p-value", {
	inf <- make_kk_ols_inference(seed = 404)
	raw_pval <- inf$compute_lik_ratio_two_sided_pval(delta = 0)
	bartlett_pval <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0)
	expect_true(is.finite(raw_pval) && is.finite(bartlett_pval))
	expect_gt(bartlett_pval, raw_pval)
})

test_that("Exact Bartlett p-value at a nonzero delta also matches lm()'s classical t-test", {
	inf <- make_kk_ols_inference(seed = 505)
	est <- inf$compute_estimate()
	delta <- est - 0.3
	ref <- lm_reference(inf, delta = delta)
	bartlett_pval <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = delta)
	expect_equal(bartlett_pval, ref$pval, tolerance = 1e-6)
})

test_that("Exact Bartlett confidence interval matches the classical t-based OLS interval", {
	inf <- make_kk_ols_inference(seed = 606)
	priv <- inf$.__enclos_env__$private
	spec <- priv$get_likelihood_test_spec()
	X <- spec$X
	rss_full <- sum((spec$y - X %*% as.numeric(spec$full_fit$b))^2)
	df_resid <- nrow(X) - ncol(X)
	sig2 <- spec$full_fit$sigma2_hat
	XtX_inv_jj <- solve(t(X) %*% X)[spec$j, spec$j]
	se_classical <- sqrt(sig2 * XtX_inv_jj)
	est <- inf$compute_estimate()
	alpha <- 0.1
	expected_ci <- est + c(-1, 1) * qt(1 - alpha / 2, df_resid) * se_classical

	bartlett_ci <- inf$compute_lik_ratio_bartlett_exact_confidence_interval(alpha = alpha)
	expect_length(bartlett_ci, 2)
	expect_true(all(is.finite(bartlett_ci)))
	expect_equal(as.numeric(bartlett_ci), expected_ci, tolerance = 1e-4)
})

test_that("Exact Bartlett factor is well-defined and positive across a range of delta, including near the estimate", {
	inf <- make_kk_ols_inference(seed = 707)
	priv <- inf$.__enclos_env__$private
	spec <- priv$get_likelihood_test_spec()
	est <- inf$compute_estimate()

	for (delta in c(est, est - 2, est + 2, est - 0.001)) {
		null_fit <- spec$fit_null(delta)
		factor <- priv$get_bartlett_factor_exact(spec = spec, delta = delta, full_fit = spec$full_fit, null_fit = null_fit)
		expect_true(is.finite(factor))
		expect_gt(factor, 0)
	}
})

test_that("Smart wrapper prefers exact over approx for InferenceContinKKOLSOneLik", {
	inf <- make_kk_ols_inference(seed = 808)
	priv <- inf$.__enclos_env__$private
	expect_true(priv$supports_bartlett_likelihood_ratio_exact())
	expect_true(priv$supports_bartlett_likelihood_ratio_approx())

	exact_pval <- inf$compute_lik_ratio_bartlett_exact_two_sided_pval(delta = 0)
	smart_pval <- expect_no_warning(inf$compute_lik_ratio_bartlett_two_sided_pval(delta = 0))
	expect_equal(smart_pval, exact_pval, tolerance = 1e-8)

	exact_ci <- inf$compute_lik_ratio_bartlett_exact_confidence_interval(alpha = 0.1)
	smart_ci <- expect_no_warning(inf$compute_lik_ratio_bartlett_confidence_interval(alpha = 0.1))
	expect_equal(as.numeric(smart_ci), as.numeric(exact_ci), tolerance = 1e-8)
})
