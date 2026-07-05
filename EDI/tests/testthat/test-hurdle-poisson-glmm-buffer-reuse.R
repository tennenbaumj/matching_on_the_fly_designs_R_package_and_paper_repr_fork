library(testthat)
library(EDI)

test_that("hurdle-Poisson GLMM scratch buffer handles unequal groups", {
	skip_if_not_installed("numDeriv")
	set.seed(35)
	group_id <- c(rep(1L, 4L), 2L, rep(3L, 3L), rep(4L, 2L))
	X <- cbind(1, rep(0:1, length.out = length(group_id)), rnorm(length(group_id)))
	y <- c(2, 3, 1, 4, 2, 5, 2, 3, 1, 4)
	params <- c(0.5, -0.15, 0.2, -1)
	n_gh <- 7L

	nll <- function(par) {
		EDI:::get_hurdle_poisson_glmm_neg_loglik_cpp(X, y, group_id, par, n_gh)
	}
	score <- EDI:::get_hurdle_poisson_glmm_score_cpp(X, y, group_id, params, n_gh)
	hessian <- EDI:::get_hurdle_poisson_glmm_hessian_cpp(X, y, group_id, params, n_gh)

	expect_equal(as.numeric(score), as.numeric(-numDeriv::grad(nll, params)), tolerance = 1e-7)
	expect_equal(unname(hessian), unname(-numDeriv::hessian(nll, params)), tolerance = 1e-5)

	perm <- sample(seq_along(group_id))
	expect_equal(
		EDI:::get_hurdle_poisson_glmm_neg_loglik_cpp(X[perm, ], y[perm], group_id[perm], params, n_gh),
		nll(params),
		tolerance = 1e-12
	)
	expect_equal(
		as.numeric(EDI:::get_hurdle_poisson_glmm_score_cpp(
			X[perm, ], y[perm], group_id[perm], params, n_gh
		)),
		as.numeric(score),
		tolerance = 1e-12
	)
})
