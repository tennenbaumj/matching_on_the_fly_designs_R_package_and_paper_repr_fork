library(testthat)
library(EDI)

test_that("hurdle-Poisson GLMM large-lambda fast path matches exact log1p", {
	skip_if_not_installed("numDeriv")
	X <- cbind(1, c(-1, 0.5, 1, -0.5, 0.25, 0.75))
	y <- c(16, 19, 18, 21, 17, 20)
	group_id <- rep(1:3, each = 2L)
	params <- c(log(18), 0.03, -4)
	n_gh <- 7L

	J <- matrix(0, n_gh, n_gh)
	for (i in seq_len(n_gh - 1L)) {
		J[i, i + 1L] <- sqrt(i / 2)
		J[i + 1L, i] <- J[i, i + 1L]
	}
	eig <- eigen(J, symmetric = TRUE)
	nodes <- eig$values
	log_weights <- log(eig$vectors[1L, ]^2)

	exact_nll <- function(par) {
		eta0 <- drop(X %*% par[1:2])
		b <- sqrt(2) * exp(par[3]) * nodes
		total <- 0
		for (g in unique(group_id)) {
			idx <- which(group_id == g)
			log_terms <- vapply(seq_along(nodes), function(k) {
				eta <- eta0[idx] + b[k]
				lambda <- exp(eta)
				log_weights[k] + sum(
					y[idx] * eta - lambda - lfactorial(y[idx]) - log1p(-exp(-lambda))
				)
			}, numeric(1))
			m <- max(log_terms)
			total <- total - (m + log(sum(exp(log_terms - m))))
		}
		total
	}

	all_lambda <- exp(drop(X %*% params[1:2]) +
		rep(sqrt(2) * exp(params[3]) * nodes, each = nrow(X)))
	expect_true(all(all_lambda > 16))

	nll_cpp <- EDI:::get_hurdle_poisson_glmm_neg_loglik_cpp(X, y, group_id, params, n_gh)
	score_cpp <- EDI:::get_hurdle_poisson_glmm_score_cpp(X, y, group_id, params, n_gh)
	hessian_cpp <- EDI:::get_hurdle_poisson_glmm_hessian_cpp(X, y, group_id, params, n_gh)

	expect_equal(nll_cpp, exact_nll(params), tolerance = 1e-12)
	expect_equal(
		as.numeric(score_cpp),
		as.numeric(-numDeriv::grad(exact_nll, params)),
		tolerance = 1e-7
	)
	expect_equal(
		unname(hessian_cpp),
		unname(-numDeriv::hessian(exact_nll, params)),
		tolerance = 1e-5
	)
})
