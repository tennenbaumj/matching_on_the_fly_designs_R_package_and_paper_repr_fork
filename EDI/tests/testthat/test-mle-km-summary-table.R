library(testthat)
library(EDI)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_design = function(response_type, n = 60L, n_cov = 2L, seed = NULL){
	if (!is.null(seed)) set.seed(seed)
	X = as.data.frame(replicate(n_cov, rnorm(n), simplify = FALSE))
	names(X) = paste0("x", seq_len(n_cov))
	des = DesignFixedBernoulli$new(n = n, response_type = response_type, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des
}

make_kk_design = function(response_type, n = 32L, n_cov = 2L, seed = NULL){
	if (!is.null(seed)) set.seed(seed)
	X_all = as.data.frame(replicate(n_cov, rnorm(n), simplify = FALSE))
	names(X_all) = paste0("x", seq_len(n_cov))
	des = DesignSeqOneByOneKK14$new(n = n, response_type = response_type, verbose = FALSE)
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X_all[i, , drop = FALSE])
	list(des = des, X = X_all)
}

# Check that a CI is a valid length-2 finite vector with lower <= upper
expect_valid_ci = function(ci){
	expect_true(is.numeric(ci))
	expect_length(ci, 2L)
	expect_true(all(is.finite(ci)))
	expect_lte(ci[1L], ci[2L])
}

# ---------------------------------------------------------------------------
# Part 1: Direct tests of InferenceMLEorKMSummaryTable's own shared() path
#
# No concrete production class currently reaches this path (they all override
# it via InferenceAsympLikStdModCache).  We exercise it with two minimal
# local subclasses:
#   1. generate_mod() returns an R model object  → vcov() S3 dispatch
#   2. generate_mod() returns a plain list        → $vcov fallback
# ---------------------------------------------------------------------------

local({
	# --- Mock 1: wraps lm() — uses stats::vcov() dispatch ---
	MockLMInf = R6::R6Class("MockLMInf",
		lock_objects = FALSE,
		inherit = EDI:::InferenceMLEorKMSummaryTable,
		public = list(
			initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
				super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			}
		),
		private = list(
			generate_mod = function(estimate_only = FALSE){
				X_cov = private$X                     # covariate columns only
				df    = cbind(data.frame(treatment = as.numeric(private$w)),
				              as.data.frame(X_cov))
				lm(as.numeric(private$y) ~ ., data = df)
			}
		)
	)

	# --- Mock 2: returns a plain list with $coefficients and $vcov ---
	MockListInf = R6::R6Class("MockListInf",
		lock_objects = FALSE,
		inherit = EDI:::InferenceMLEorKMSummaryTable,
		public = list(
			initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
				super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			}
		),
		private = list(
			generate_mod = function(estimate_only = FALSE){
				X_cov = private$X
				X     = cbind(`(Intercept)` = 1, treatment = as.numeric(private$w),
				              as.matrix(X_cov))
				y     = as.numeric(private$y)
				b     = as.numeric(qr.solve(X, y))
				names(b) = colnames(X)
				n    = length(y); p = length(b)
				s2   = sum((y - X %*% b)^2) / (n - p)
				vcov = s2 * solve(crossprod(X))
				list(coefficients = b, vcov = vcov)
			}
		)
	)

	setup_continuous_des = function(seed){
		des = make_design("continuous", n = 60L, n_cov = 2L, seed = seed)
		w   = des$get_w()
		X   = as.matrix(des$get_X())
		y   = 0.4 * w + 0.3 * X[, 1L] - 0.2 * X[, 2L] + rnorm(60L)
		des$add_all_subject_responses(y)
		des
	}

	test_that("MockLMInf: R model object path — lazy vcov, single fit, covariate SEs", {
		des = setup_continuous_des(seed = 101L)
		inf = MockLMInf$new(des, model_formula = ~.)

		# compute_estimate() should cache the model but NOT extract vcov yet
		est = inf$compute_estimate(estimate_only = TRUE)
		expect_true(is.finite(est))
		priv = inf$.__enclos_env__$private
		expect_false(is.null(priv$cached_mod))         # model cached
		expect_null(priv$cached_values$summary_table)  # vcov NOT yet extracted

		# compute_asymp_confidence_interval() triggers lazy vcov extraction
		ci = inf$compute_asymp_confidence_interval(alpha = 0.05)
		expect_valid_ci(ci)

		# summary_table must now be populated
		st = priv$cached_values$summary_table
		expect_false(is.null(st))
		expect_equal(colnames(st), c("Value", "Std. Error", "z value", "Pr(>|z|)"))
		expect_true("treatment" %in% rownames(st))
		expect_true("x1"        %in% rownames(st))
		expect_true("x2"        %in% rownames(st))
		expect_true(is.finite(st["treatment", "Std. Error"]))
		expect_true(is.finite(st["x1",        "Std. Error"]))
		expect_true(is.finite(st["x2",        "Std. Error"]))

		# Second CI call must hit the cache — same result, no second fit
		ci2 = inf$compute_asymp_confidence_interval(alpha = 0.05)
		expect_equal(ci, ci2)

		# get_summary() must not throw
		expect_no_error(suppressWarnings(inf$get_summary()))
	})

	test_that("MockListInf: plain-list $vcov fallback — covariate SEs", {
		des = setup_continuous_des(seed = 102L)
		inf = MockListInf$new(des, model_formula = ~.)

		ci = inf$compute_asymp_confidence_interval(alpha = 0.05)
		expect_valid_ci(ci)

		st = inf$.__enclos_env__$private$cached_values$summary_table
		expect_false(is.null(st))
		expect_true("x1" %in% rownames(st))
		expect_true("x2" %in% rownames(st))
		expect_true(is.finite(st["x1", "Std. Error"]))
		expect_true(is.finite(st["x2", "Std. Error"]))
	})

	test_that("MockLMInf: compute_estimate() then pval does not double-fit", {
		des  = setup_continuous_des(seed = 103L)
		inf  = MockLMInf$new(des, model_formula = ~.)
		priv = inf$.__enclos_env__$private

		inf$compute_estimate()
		mod_ref = priv$cached_mod             # pointer to first fit
		inf$compute_asymp_two_sided_pval()
		expect_identical(priv$cached_mod, mod_ref)  # same object, no second fit
	})
})

# ---------------------------------------------------------------------------
# Part 2: All concrete InferenceAsympLikStdModCache families with ~.
#
# These classes inherit InferenceMLEorKMSummaryTable through InferenceAsympLik.
# They use InferenceAsympLikStdModCache's own shared(), so summary_table is
# not populated via InferenceMLEorKMSummaryTable — but the inheritance chain
# still carries through. We verify:
#   • model_formula = ~. is accepted (no init error)
#   • compute_asymp_confidence_interval() returns a finite CI
#   • get_summary() does not throw
# ---------------------------------------------------------------------------

test_that("InferenceCountPoisson: model_formula = ~.", {
	set.seed(201L)
	des = make_design("count", seed = 201L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	des$add_all_subject_responses(rpois(60L, exp(0.3 * w + 0.2 * X[, 1L])))
	inf = InferenceCountPoisson$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceCountNegBin: model_formula = ~.", {
	set.seed(202L)
	des = make_design("count", seed = 202L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	des$add_all_subject_responses(rnbinom(60L, size = 3, mu = exp(0.3 * w + 0.2 * X[, 1L])))
	inf = InferenceCountNegBin$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceCountHurdlePoisson: model_formula = ~.", {
	set.seed(203L)
	des = make_design("count", seed = 203L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	y   = rpois(60L, exp(0.3 * w + 0.2 * X[, 1L]))
	des$add_all_subject_responses(y)
	inf = InferenceCountHurdlePoisson$new(des, model_formula = ~., use_rcpp = TRUE, verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceIncidLogRegr: model_formula = ~.", {
	set.seed(211L)
	des = make_design("incidence", seed = 211L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	des$add_all_subject_responses(rbinom(60L, 1L, plogis(0.4 * w + 0.3 * X[, 1L])))
	inf = InferenceIncidLogRegr$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceIncidRiskDiff: model_formula = ~.", {
	set.seed(212L)
	des = make_design("incidence", seed = 212L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	des$add_all_subject_responses(rbinom(60L, 1L, plogis(0.4 * w + 0.3 * X[, 1L])))
	inf = InferenceIncidRiskDiff$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceOrdinalPropOddsRegr: model_formula = ~.", {
	set.seed(221L)
	des = make_design("ordinal", n = 80L, seed = 221L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	y   = as.integer(cut(0.4 * w + 0.3 * X[, 1L] + rnorm(80L),
	                     breaks = c(-Inf, -0.8, 0, 0.8, Inf)))
	des$add_all_subject_responses(y)
	inf = InferenceOrdinalPropOddsRegr$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferencePropBetaRegr: model_formula = ~.", {
	set.seed(231L)
	des = make_design("proportion", seed = 231L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	mu  = plogis(0.3 * w + 0.2 * X[, 1L])
	y   = pmax(pmin(rbeta(60L, mu * 8, (1 - mu) * 8), 1 - 1e-6), 1e-6)
	des$add_all_subject_responses(y)
	inf = InferencePropBetaRegr$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferencePropZeroOneInflatedBetaRegr: model_formula = ~. (no-error, convergence not required)", {
	# ZOIB has two sub-models (beta + inflation), so ~. with 2 covariates doubles the
	# parameter count; convergence is not guaranteed on small datasets.  We only verify
	# the call doesn't throw.
	set.seed(232L)
	des = make_design("proportion", seed = 232L)
	w   = des$get_w(); X = as.matrix(des$get_X())
	mu  = plogis(0.3 * w + 0.2 * X[, 1L])
	y   = pmax(pmin(rbeta(60L, mu * 8, (1 - mu) * 8), 1 - 1e-6), 1e-6)
	des$add_all_subject_responses(y)
	inf = InferencePropZeroOneInflatedBetaRegr$new(des, model_formula = ~., verbose = FALSE)
	expect_no_error(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceSurvivalCoxPHRegr: model_formula = ~.", {
	skip_if_not_installed("survival")
	set.seed(241L)
	des  = make_design("survival", seed = 241L)
	w    = des$get_w(); X = as.matrix(des$get_X())
	time = rexp(60L, exp(-0.5 + 0.3 * w + 0.2 * X[, 1L]))
	des$add_all_subject_responses(time, rep(1L, 60L))
	inf = InferenceSurvivalCoxPHRegr$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceSurvivalWeibullRegr: model_formula = ~.", {
	set.seed(242L)
	des  = make_design("survival", seed = 242L)
	w    = des$get_w(); X = as.matrix(des$get_X())
	time = rexp(60L, exp(-0.5 + 0.3 * w + 0.2 * X[, 1L]))
	des$add_all_subject_responses(time, rep(1L, 60L))
	inf = InferenceSurvivalWeibullRegr$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceCountPoisson (KK design): model_formula = ~.", {
	set.seed(251L)
	kk  = make_kk_design("count", n = 40L, n_cov = 2L, seed = 251L)
	des = kk$des; X = as.matrix(kk$X)
	for (i in seq_len(40L)){
		w_i = des$get_w()[i]
		des$add_one_subject_response(i, rpois(1L, exp(0.25 * w_i + 0.2 * X[i, 1L])), 1L)
	}
	inf = InferenceCountPoisson$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})

test_that("InferenceIncidLogRegr (KK design): model_formula = ~.", {
	set.seed(252L)
	kk  = make_kk_design("incidence", n = 40L, n_cov = 2L, seed = 252L)
	des = kk$des; X = as.matrix(kk$X)
	for (i in seq_len(40L)){
		w_i = des$get_w()[i]
		des$add_one_subject_response(i, rbinom(1L, 1L, plogis(0.3 * w_i + 0.2 * X[i, 1L])), 1L)
	}
	inf = InferenceIncidLogRegr$new(des, model_formula = ~., verbose = FALSE)
	expect_valid_ci(inf$compute_asymp_confidence_interval())
	expect_no_error(suppressWarnings(inf$get_summary()))
})
