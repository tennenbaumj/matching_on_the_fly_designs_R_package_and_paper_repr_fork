library(testthat)
library(EDI)

# Verify that warm-start chaining in compute_weighted_gcomp_estimate / weighted_gcomp_fit
# produces results consistent with independent (fresh-object) fits.
# Logistic regression is strictly convex → both paths converge to the same MLE within tolerance.

make_incid_kk_design <- function(n, seed) {
	set.seed(seed)
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	w = des$get_w()
	p = plogis(-0.3 + 0.7 * w + 0.4 * X$x1)
	des$add_all_subject_responses(rbinom(n, 1, p))
	des
}

test_that("KK gcomp warm-start chaining matches independent fits (RD)", {
	des = make_incid_kk_design(60L, seed = 20260705)
	set.seed(1L)
	n = des$get_t()
	weight_sets = lapply(1:15, function(b) { wts = rexp(n); wts / mean(wts) })

	# Chained: reuse one object so gcomp_boot_beta accumulates across calls
	inf_chain = InferenceIncidKKGCompRiskDiff$new(des, verbose = FALSE)
	priv_chain = inf_chain$.__enclos_env__$private
	priv_chain$shared(estimate_only = TRUE)  # run primary fit so warm start is set
	chain_vals = vapply(weight_sets, function(wts) {
		rd = priv_chain$compute_weighted_gcomp_estimate(wts)
		if (is.null(rd)) NA_real_ else rd
	}, numeric(1))

	# Reference: fresh object per replicate (no chaining, cold from primary warm start)
	ref_vals = vapply(weight_sets, function(wts) {
		inf_fresh = InferenceIncidKKGCompRiskDiff$new(des, verbose = FALSE)
		priv_fresh = inf_fresh$.__enclos_env__$private
		priv_fresh$shared(estimate_only = TRUE)
		rd = priv_fresh$compute_weighted_gcomp_estimate(wts)
		if (is.null(rd)) NA_real_ else rd
	}, numeric(1))

	expect_equal(chain_vals, ref_vals, tolerance = 1e-6)
	expect_true(all(is.finite(chain_vals)))
})

test_that("KK gcomp warm-start chaining matches independent fits (RR)", {
	des = make_incid_kk_design(60L, seed = 20260706)
	set.seed(2L)
	n = des$get_t()
	weight_sets = lapply(1:15, function(b) { wts = rexp(n); wts / mean(wts) })

	inf_chain = InferenceIncidKKGCompRiskRatio$new(des, verbose = FALSE)
	priv_chain = inf_chain$.__enclos_env__$private
	priv_chain$shared(estimate_only = TRUE)
	chain_vals = vapply(weight_sets, function(wts) {
		rr = priv_chain$compute_weighted_gcomp_estimate(wts)
		if (is.null(rr)) NA_real_ else rr
	}, numeric(1))

	ref_vals = vapply(weight_sets, function(wts) {
		inf_fresh = InferenceIncidKKGCompRiskRatio$new(des, verbose = FALSE)
		priv_fresh = inf_fresh$.__enclos_env__$private
		priv_fresh$shared(estimate_only = TRUE)
		rr = priv_fresh$compute_weighted_gcomp_estimate(wts)
		if (is.null(rr)) NA_real_ else rr
	}, numeric(1))

	expect_equal(chain_vals, ref_vals, tolerance = 1e-6)
	expect_true(all(is.finite(chain_vals)))
})

test_that("non-KK gcomp warm-start chaining matches independent fits via weighted_gcomp_fit", {
	set.seed(20260707)
	n = 80L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	des$add_all_subject_responses(rbinom(n, 1, plogis(-0.2 + 0.6 * w + 0.3 * X$x1)))

	set.seed(3L)
	weight_sets = lapply(1:15, function(b) { wts = rexp(n); wts / mean(wts) })

	# Chained: reuse one object
	inf_chain = InferenceIncidGCompRiskDiff$new(des, verbose = FALSE)
	priv_chain = inf_chain$.__enclos_env__$private
	priv_chain$shared(estimate_only = TRUE)
	chain_vals = vapply(weight_sets, function(wts) {
		effects = priv_chain$weighted_gcomp_effects_from_row_weights(wts)
		if (is.null(effects)) NA_real_ else effects$rd
	}, numeric(1))

	# Reference: fresh object per replicate
	ref_vals = vapply(weight_sets, function(wts) {
		inf_fresh = InferenceIncidGCompRiskDiff$new(des, verbose = FALSE)
		priv_fresh = inf_fresh$.__enclos_env__$private
		priv_fresh$shared(estimate_only = TRUE)
		effects = priv_fresh$weighted_gcomp_effects_from_row_weights(wts)
		if (is.null(effects)) NA_real_ else effects$rd
	}, numeric(1))

	expect_equal(chain_vals, ref_vals, tolerance = 1e-6)
	expect_true(all(is.finite(chain_vals)))
})
