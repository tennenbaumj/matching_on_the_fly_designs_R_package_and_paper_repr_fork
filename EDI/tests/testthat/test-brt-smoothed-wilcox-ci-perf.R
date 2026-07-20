library(EDI)

SlowInferenceAllSimpleWilcox = R6::R6Class(
	"SlowInferenceAllSimpleWilcox",
	inherit = InferenceAllSimpleWilcox,
	lock_objects = FALSE,
	private = list(
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			NULL
		}
	)
)

test_that("smoothed CI: fast-kernel result matches the forced-slow-fallback result", {
	set.seed(20260731)
	n = 30
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(rnorm(n))

	fast_inf = InferenceAllSimpleWilcox$new(des)
	slow_inf = SlowInferenceAllSimpleWilcox$new(des)

	set.seed(99)
	fast_ci = fast_inf$compute_rand_bootstrap_confidence_interval(B = 99, type = "smoothed", show_progress = FALSE)
	set.seed(99)
	slow_ci = slow_inf$compute_rand_bootstrap_confidence_interval(B = 99, type = "smoothed", show_progress = FALSE)

	expect_equal(as.numeric(fast_ci), as.numeric(slow_ci), tolerance = 1e-6)
})

test_that("smoothed CI: fast kernel is dramatically faster than the forced-slow fallback", {
	# This is the actual motivating case from the original investigation: CI inversion
	# pre-materializes fresh assignments (materialize_w = TRUE) once and reuses them (common
	# random numbers) across every delta evaluated during root-finding, so the fast kernel can
	# engage on every one of those evaluations. A *standalone* compute_rand_bootstrap_two_sided_pval
	# call (no CI inversion) deliberately draws w lazily per-iteration (materialize_w = FALSE) for
	# CRN reasons unrelated to smoothing, so rand_bootstrap_draw_matrices() can't build i_mat/w_mat
	# and the fast kernel never engages there regardless of this fix — that path's cost is
	# untouched by design, not a regression, so it is not asserted on here.
	set.seed(20260732)
	n = 30
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(rnorm(n))

	fast_inf = InferenceAllSimpleWilcox$new(des)
	slow_inf = SlowInferenceAllSimpleWilcox$new(des)

	set.seed(99)
	t_fast = system.time(
		fast_ci <- fast_inf$compute_rand_bootstrap_confidence_interval(B = 99, type = "smoothed", show_progress = FALSE)
	)[["elapsed"]]
	set.seed(99)
	t_slow = system.time(
		slow_ci <- slow_inf$compute_rand_bootstrap_confidence_interval(B = 99, type = "smoothed", show_progress = FALSE)
	)[["elapsed"]]

	expect_equal(as.numeric(fast_ci), as.numeric(slow_ci), tolerance = 1e-6)
	# Measured ~50x on this workload during development; require at least 5x here to allow for
	# machine variance while still catching a regression back to the fully-slow path.
	expect_true(t_fast < t_slow / 5, info = sprintf("t_fast=%.3fs t_slow=%.3fs", t_fast, t_slow))
})
