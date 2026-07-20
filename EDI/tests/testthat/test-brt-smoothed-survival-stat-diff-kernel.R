library(EDI)

test_that("compute_survival_stat_diff_rand_bootstrap_parallel_cpp adds noise before the multiplicative sharp-null shift (both do_rmst values)", {
	set.seed(20260730)
	n = 60; B = 4
	y0 = rexp(n, rate = 0.2) + 0.5
	dead = rbinom(n, 1, 0.8)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.1), nrow = n, ncol = B)
	delta = 0.2

	for (do_rmst in c(TRUE, FALSE)) {
		res_noisy = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
			as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, do_rmst, noise_mat, 1L
		)
		res_ref = numeric(B)
		for (b in seq_len(B)) {
			res_ref[b] = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
				as.numeric(y0 + noise_mat[, b]), as.integer(dead), i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
				delta, do_rmst, NULL, 1L
			)
		}
		expect_equal(res_noisy, res_ref, tolerance = 1e-8, info = paste("do_rmst =", do_rmst))

		res_null = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, do_rmst, NULL, 1L)
		res_zero = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, do_rmst, matrix(0, n, B), 1L)
		expect_equal(res_null, res_zero, tolerance = 1e-10, info = paste("do_rmst =", do_rmst))
	}
})
