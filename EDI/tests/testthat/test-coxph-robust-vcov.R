library(testthat)
library(EDI)

test_that("clustered Cox robust covariance matches survival::coxph", {
	skip_if_not_installed("survival")
	set.seed(30)
	n <- 240L
	p <- 3L
	X <- matrix(rnorm(n * p), nrow = n, ncol = p)
	colnames(X) <- paste0("x", seq_len(p))
	beta <- c(0.35, -0.2, 0.1)
	y <- round(rexp(n, rate = 0.1 * exp(drop(X %*% beta))), digits = 1L)
	y[y == 0] <- 0.1
	dead <- rbinom(n, 1L, 0.8)
	cluster <- sample(rep(seq_len(60L), each = 4L))

	res_cpp <- fast_coxph_regression_cpp(X, y, dead, cluster = cluster)
	res_r <- survival::coxph(
		survival::Surv(y, dead) ~ x1 + x2 + x3 + survival::cluster(cluster),
		data = data.frame(y, dead, cluster, X),
		ties = "breslow",
		robust = TRUE
	)

	expect_equal(
		as.numeric(res_cpp$coefficients),
		as.numeric(stats::coef(res_r)),
		tolerance = 1e-7
	)
	expect_equal(
		unname(res_cpp$vcov),
		unname(stats::vcov(res_r)),
		tolerance = 1e-7
	)
})
