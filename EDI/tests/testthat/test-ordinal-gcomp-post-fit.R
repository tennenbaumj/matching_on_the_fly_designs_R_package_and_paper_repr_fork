library(testthat)
library(EDI)

test_that("gcomp_ordinal_proportional_odds_post_fit_cpp mean1/mean0/md match R reference", {
	set.seed(77)
	n <- 300L
	X <- cbind(treatment = rep(0:1, length.out = n), matrix(rnorm(n * 3L), ncol = 3L))
	latent <- drop(X %*% c(0.5, -0.2, 0.3, -0.1)) + stats::rlogis(n)
	y <- as.numeric(cut(latent, breaks = c(-Inf, -1, 0.5, Inf)))
	fit <- fast_ordinal_regression_cpp(X, y, estimate_only = TRUE)
	coef_hat <- as.numeric(fit$b)
	alpha_hat <- as.numeric(fit$alpha)

	result <- gcomp_ordinal_proportional_odds_post_fit_cpp(
		X_fit_sexp = X, coef_hat_sexp = coef_hat, alpha_hat_sexp = alpha_hat, j_treat = 1L
	)

	mean_ref <- function(design) {
		eta <- drop(design %*% coef_hat)
		mean(1 + rowSums(vapply(alpha_hat, function(a) stats::plogis(eta - a), numeric(n))))
	}
	X1 <- X; X1[, 1L] <- 1
	X0 <- X; X0[, 1L] <- 0
	expect_equal(result$mean1, mean_ref(X1), tolerance = 1e-12)
	expect_equal(result$mean0, mean_ref(X0), tolerance = 1e-12)
	expect_equal(result$md, mean_ref(X1) - mean_ref(X0), tolerance = 1e-12)
})

test_that("ordinal G-computation post-fit results match an R reference", {
	set.seed(32)
	n <- 240L
	X <- cbind(treatment = rep(0:1, length.out = n), matrix(rnorm(n * 2L), ncol = 2L))
	latent <- drop(X %*% c(0.45, -0.3, 0.2)) + stats::rlogis(n)
	y <- as.numeric(cut(latent, breaks = c(-Inf, -1, 0, 1, Inf)))
	fit <- fast_ordinal_regression_with_var_cpp(X, y)

	result <- ordinal_gcomp_post_fit_cpp(
		X,
		y,
		coef_hat = as.numeric(fit$b),
		alpha_hat = as.numeric(fit$alpha),
		j_treat = 1L
	)

	n_alpha <- length(fit$alpha)
	theta <- c(fit$alpha, fit$b)
	X1 <- X
	X0 <- X
	X1[, 1L] <- 1
	X0[, 1L] <- 0
	mean_from_theta <- function(par, design) {
		alpha <- par[seq_len(n_alpha)]
		beta <- par[n_alpha + seq_len(ncol(X))]
		eta <- drop(design %*% beta)
		mean(1 + rowSums(vapply(alpha, function(a) stats::plogis(eta - a), numeric(n))))
	}
	md_from_theta <- function(par) mean_from_theta(par, X1) - mean_from_theta(par, X0)
	h <- 1e-5
	grad <- vapply(seq_along(theta), function(j) {
		p_plus <- theta
		p_minus <- theta
		p_plus[j] <- p_plus[j] + h
		p_minus[j] <- p_minus[j] - h
		(md_from_theta(p_plus) - md_from_theta(p_minus)) / (2 * h)
	}, numeric(1))
	se_md <- sqrt(drop(t(grad) %*% fit$vcov %*% grad))
	beta_idx <- n_alpha + seq_along(fit$b)

	expect_equal(result$mean1, mean_from_theta(theta, X1), tolerance = 1e-12)
	expect_equal(result$mean0, mean_from_theta(theta, X0), tolerance = 1e-12)
	expect_equal(result$md, md_from_theta(theta), tolerance = 1e-12)
	expect_equal(result$se_md, se_md, tolerance = 1e-10)
	expect_equal(result$vcov, fit$vcov[beta_idx, beta_idx], tolerance = 1e-12)
})
