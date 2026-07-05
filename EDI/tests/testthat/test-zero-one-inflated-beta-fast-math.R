library(testthat)
library(EDI)

test_that("ZOIB fast math and reusable buffers preserve score and Hessian", {
	skip_if_not_installed("numDeriv")
	set.seed(42)
	n <- 24L
	X <- cbind(1, rnorm(n))
	X_zi <- cbind(1, rep(c(-1, 1), length.out = n))
	y <- c(0, 1, 0, 1, seq(0.08, 0.92, length.out = n - 4L))
	params <- c(0.15, -0.2, log(8), -1.1, 0.25, -1.3, -0.15)

	nll <- function(par) {
		p <- ncol(X)
		p_zi <- ncol(X_zi)
		mu <- plogis(drop(X %*% par[seq_len(p)]))
		phi <- exp(par[p + 1L])
		eta0 <- drop(X_zi %*% par[p + 1L + seq_len(p_zi)])
		eta1 <- drop(X_zi %*% par[p + 1L + p_zi + seq_len(p_zi)])
		denom <- 1 + exp(eta0) + exp(eta1)
		pi0 <- exp(eta0) / denom
		pi1 <- exp(eta1) / denom
		pib <- 1 / denom
		is_zero <- y <= 0
		is_one <- y >= 1
		is_beta <- !(is_zero | is_one)
		loglik <- sum(log(pi0[is_zero])) + sum(log(pi1[is_one])) + sum(log(pib[is_beta]))
		loglik <- loglik + sum(dbeta(y[is_beta], mu[is_beta] * phi,
			(1 - mu[is_beta]) * phi, log = TRUE))
		-loglik
	}

	score <- EDI:::get_zero_one_inflated_beta_score_cpp(X, X_zi, y, params)
	hessian <- EDI:::get_zero_one_inflated_beta_hessian_cpp(X, X_zi, y, params)
	expect_equal(as.numeric(score), as.numeric(-numDeriv::grad(nll, params)), tolerance = 1e-6)
	expect_equal(unname(hessian), unname(-numDeriv::hessian(nll, params)), tolerance = 1e-4)
})
