library(EDI)

test_that("rand_bootstrap_draw_matrices packs a noise_mat only when draws carry smooth_noise", {
	set.seed(20260721)
	n = 10
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(rnorm(n))
	inf = InferenceAllSimpleWilcox$new(des)
	priv = inf$.__enclos_env__$private

	draws = priv$generate_rand_bootstrap_draws(B = 4L, materialize_w = TRUE)
	mats_no_noise = priv$rand_bootstrap_draw_matrices(draws)
	expect_null(mats_no_noise$noise_mat)

	for (b in seq_along(draws)) draws[[b]][["smooth_noise"]] = rep(b, n)
	mats_noise = priv$rand_bootstrap_draw_matrices(draws)
	expect_true(is.matrix(mats_noise$noise_mat))
	expect_equal(dim(mats_noise$noise_mat), c(n, 4L))
	expect_equal(mats_noise$noise_mat[, 3], rep(3, n))

	# all-or-nothing contract: one draw missing smooth_noise -> NULL, like w_b today
	draws2 = draws
	draws2[[2]][["smooth_noise"]] = NULL
	expect_null(priv$rand_bootstrap_draw_matrices(draws2)$noise_mat)
})

test_that("Ridit and Jonckheere-Terpstra fast kernels decline smoothed (noise) draws", {
	set.seed(20260722)
	n = 12
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "ordinal")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(sample(1:4, n, replace = TRUE))

	inf_ridit = InferenceOrdinalRidit$new(des)
	priv_ridit = inf_ridit$.__enclos_env__$private
	priv_ridit$shared()
	draws = priv_ridit$generate_rand_bootstrap_draws(B = 3L, materialize_w = TRUE)
	for (b in seq_along(draws)) draws[[b]][["smooth_noise"]] = rnorm(n)
	expect_null(priv_ridit$compute_fast_rand_bootstrap_distr(as.numeric(priv_ridit$y), draws, 0, "none"))

	inf_jt = InferenceOrdinalJonckheereTerpstraTest$new(des)
	priv_jt = inf_jt$.__enclos_env__$private
	draws_jt = priv_jt$generate_rand_bootstrap_draws(B = 3L, materialize_w = TRUE)
	for (b in seq_along(draws_jt)) draws_jt[[b]][["smooth_noise"]] = rnorm(n)
	expect_null(priv_jt$compute_fast_rand_bootstrap_distr(as.numeric(priv_jt$y), draws_jt, 0, "none"))
})
