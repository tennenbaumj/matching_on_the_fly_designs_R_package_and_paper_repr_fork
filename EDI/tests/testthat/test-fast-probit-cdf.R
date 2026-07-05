library(testthat)
library(EDI)

test_that("fast probit CDF matches the normal CDF across the fitted range", {
	grid <- seq(-7.9, 7.9, by = 0.01)
	recovered <- numeric(length(grid))
	block_size <- 160L
	for (start in seq.int(1L, length(grid), by = block_size)) {
		idx <- start:min(start + block_size - 1L, length(grid))
		eta <- grid[idx]
		X <- diag(length(idx))
		score <- EDI:::get_probit_regression_score_cpp(X, rep(1, length(idx)), eta)
		recovered[idx] <- stats::dnorm(eta) / as.numeric(score)
	}

	reference <- stats::pnorm(grid)
	expect_lt(max(abs(recovered - reference)), 2e-15)
	tail_idx <- which(grid < -5)
	expect_lt(
		max(abs(log(recovered[tail_idx]) - stats::pnorm(grid[tail_idx], log.p = TRUE))),
		2e-12
	)
})

test_that("ordinal probit fast-erfc score and Hessian match numerical derivatives", {
	skip_if_not_installed("numDeriv")
	set.seed(5101)
	n <- 90L; p <- 3L
	X <- matrix(rnorm(n * p), n, p)
	y <- rep(1:3, length.out = n)
	params <- c(-0.8, 1.1, 0.25, -0.4, 0.15)

	nll <- function(par) {
		K <- 3L; n_alpha <- K - 1L
		alpha <- par[seq_len(n_alpha)]
		beta  <- par[seq(n_alpha + 1L, length(par))]
		eta   <- as.numeric(X %*% beta)
		total <- 0
		for (i in seq_len(n)) {
			yi <- y[i]
			z_u <- if (yi == K) Inf else alpha[yi] - eta[i]
			z_l <- if (yi == 1L) -Inf else alpha[yi - 1L] - eta[i]
			prob <- max(pnorm(z_u) - pnorm(z_l), 1e-12)
			total <- total - log(prob)
		}
		total
	}

	score   <- EDI:::get_ordinal_probit_regression_score_cpp(X, y, params)
	hessian <- EDI:::get_ordinal_probit_regression_hessian_cpp(X, y, params)

	expect_equal(as.numeric(score),   as.numeric(-numDeriv::grad(nll, params)),    tolerance = 1e-7)
	expect_equal(unname(hessian),     unname(-numDeriv::hessian(nll, params)),      tolerance = 1e-5)
})
