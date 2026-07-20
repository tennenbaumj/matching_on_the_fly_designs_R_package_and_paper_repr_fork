library(EDI)

test_that("compute_rand_bootstrap_mean_diff_parallel_cpp adds noise before the sharp-null shift", {
	set.seed(20260724)
	n = 30; B = 5
	y0 = rnorm(n, sd = 2)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.5), nrow = n, ncol = B)
	delta = -0.2
	transform_code = 0L
	clamp = .Machine$double.eps

	res_noisy = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, noise_mat, 1L
	)
	res_ref = numeric(B)
	for (b in seq_len(B)) {
		res_ref[b] = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
			as.numeric(y0 + noise_mat[, b]), i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, transform_code, clamp, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-10)

	res_null = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, NULL, 1L
	)
	res_zero = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, matrix(0, n, B), 1L
	)
	expect_equal(res_null, res_zero, tolerance = 1e-12)
})
