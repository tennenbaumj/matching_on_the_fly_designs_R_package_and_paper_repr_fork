compare_bootstrap_fast_slow_asymp <- function(fast_inf, slow_inf, B = 5L, seed = 1L, tolerance = 1e-10){
	fast_inf$num_cores = 1L
	slow_inf$num_cores = 1L
	set.seed(seed)
	fast_boot = fast_inf$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE)
	set.seed(seed)
	slow_boot = slow_inf$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE)
	expect_equal(unname(fast_boot), unname(slow_boot), tolerance = tolerance)
}

make_fixed_design <- function(response_type, X, y_fun, dead_fun = NULL){
	des = DesignFixedBernoulli$new(n = nrow(X), response_type = response_type, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = y_fun(w, X)
	if (is.null(dead_fun)) {
		des$add_all_subject_responses(y)
	} else {
		des$add_all_subject_responses(y, deads = dead_fun(w, X, y))
	}
	des
}

make_fixed_blocked_cluster_design <- function(X, y_fun, cluster_size = 2L){
	X_design = as.data.frame(X)
	strata_col = "block"
	cluster_col = ".assignment_only_cluster_id"
	cluster_ids = integer(nrow(X_design))
	next_cluster_id = 1L

	for (block in unique(X_design[[strata_col]])) {
		idx = which(X_design[[strata_col]] == block)
		cluster_ids[idx] = ((seq_along(idx) - 1L) %/% cluster_size) + next_cluster_id
		next_cluster_id = max(cluster_ids[idx]) + 1L
	}

	X_design[[cluster_col]] = factor(cluster_ids)
	des = DesignFixedBlockedCluster$new(
		n = nrow(X_design),
		strata_cols = strata_col,
		cluster_col = cluster_col,
		response_type = "continuous",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(X_design)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	des$add_all_subject_responses(y_fun(w, X_design))
	des
}

test_that("incidence and ordinal g-computation reusable workers match generic bootstrap", {
	SlowInferenceIncidGCompRiskDiff = R6::R6Class(
		"SlowInferenceIncidGCompRiskDiff",
		inherit = InferenceIncidGCompRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidGCompRiskDiff = R6::R6Class(
		"SlowInferenceIncidGCompRiskDiff",
		inherit = InferenceIncidGCompRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidGCompRiskRatio = R6::R6Class(
		"SlowInferenceIncidGCompRiskRatio",
		inherit = InferenceIncidGCompRiskRatio,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidGCompRiskRatio = R6::R6Class(
		"SlowInferenceIncidGCompRiskRatio",
		inherit = InferenceIncidGCompRiskRatio,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceOrdinalGCompMeanDiff = R6::R6Class(
		"SlowInferenceOrdinalGCompMeanDiff",
		inherit = InferenceOrdinalGCompMeanDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceOrdinalGCompMeanDiff = R6::R6Class(
		"SlowInferenceOrdinalGCompMeanDiff",
		inherit = InferenceOrdinalGCompMeanDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260415)
	n = 42L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

	incid_des = make_fixed_design(
		"incidence",
		X,
		y_fun = function(w, X){
			stats::rbinom(nrow(X), 1, stats::plogis(-0.6 + 0.9 * w + 0.4 * X$x1 - 0.2 * X$x2))
		}
	)

	ordinal_des = make_fixed_design(
		"ordinal",
		X,
		y_fun = function(w, X){
			latent = -0.3 + 0.7 * w + 0.4 * X$x1 - 0.2 * X$x2 + stats::rlogis(nrow(X))
			as.integer(cut(latent, breaks = c(-Inf, -0.6, 0.2, 1.0, Inf), labels = FALSE))
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceIncidGCompRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidGCompRiskDiff$new(incid_des, verbose = FALSE),
		seed = 201
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidGCompRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidGCompRiskDiff$new(incid_des, verbose = FALSE),
		seed = 202
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidGCompRiskRatio$new(incid_des, verbose = FALSE),
		SlowInferenceIncidGCompRiskRatio$new(incid_des, verbose = FALSE),
		seed = 203
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidGCompRiskRatio$new(incid_des, verbose = FALSE),
		SlowInferenceIncidGCompRiskRatio$new(incid_des, verbose = FALSE),
		seed = 204
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceOrdinalGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		SlowInferenceOrdinalGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		seed = 205
	)
	expect_false(inherits(InferenceOrdinalGCompMeanDiff$new(ordinal_des, verbose = FALSE), "InferenceOrdinalGCompAbstract"))
	compare_bootstrap_fast_slow_asymp(
		InferenceOrdinalGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		SlowInferenceOrdinalGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		seed = 206
	)
})

test_that("continuous lin, count negbin, and classical incidence estimators match generic bootstrap", {
	SlowInferenceContinLin = R6::R6Class(
		"SlowInferenceContinLin",
		inherit = InferenceContinLin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountNegBin = R6::R6Class(
		"SlowInferenceCountNegBin",
		inherit = InferenceCountNegBin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountNegBin = R6::R6Class(
		"SlowInferenceCountNegBin",
		inherit = InferenceCountNegBin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountHurdleNegBin = R6::R6Class(
		"SlowInferenceCountHurdleNegBin",
		inherit = InferenceCountHurdleNegBin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountHurdleNegBin = R6::R6Class(
		"SlowInferenceCountHurdleNegBin",
		inherit = InferenceCountHurdleNegBin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidRiskDiff = R6::R6Class(
		"SlowInferenceIncidRiskDiff",
		inherit = InferenceIncidRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidNewcombeRiskDiff = R6::R6Class(
		"SlowInferenceIncidNewcombeRiskDiff",
		inherit = InferenceIncidNewcombeRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidMiettinenNurminenRiskDiff = R6::R6Class(
		"SlowInferenceIncidMiettinenNurminenRiskDiff",
		inherit = InferenceIncidMiettinenNurminenRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260416)
	n = 44L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

	continuous_des = make_fixed_design(
		"continuous",
		X,
		y_fun = function(w, X){
			0.5 * w + 0.4 * X$x1 - 0.3 * X$x2 + stats::rnorm(nrow(X), sd = 0.8)
		}
	)
	count_des = make_fixed_design(
		"count",
		X,
		y_fun = function(w, X){
			stats::rnbinom(nrow(X), mu = exp(0.2 + 0.45 * w + 0.25 * X$x1), size = 2.2)
		}
	)
	hurdle_des = make_fixed_design(
		"count",
		X,
		y_fun = function(w, X){
			is_zero = stats::rbinom(nrow(X), 1, stats::plogis(0.2 - 0.8 * w + 0.2 * X$x2))
			pos = 1 + stats::rnbinom(nrow(X), mu = exp(0.1 + 0.35 * w + 0.2 * X$x1), size = 1.8)
			ifelse(is_zero == 1, 0, pos)
		}
	)
	incid_des = make_fixed_design(
		"incidence",
		X,
		y_fun = function(w, X){
			stats::rbinom(nrow(X), 1, stats::plogis(-0.5 + 1.0 * w))
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceContinLin$new(continuous_des, verbose = FALSE),
		SlowInferenceContinLin$new(continuous_des, verbose = FALSE),
		seed = 207
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountNegBin$new(count_des, verbose = FALSE),
		SlowInferenceCountNegBin$new(count_des, verbose = FALSE),
		seed = 208
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountNegBin$new(count_des, verbose = FALSE),
		SlowInferenceCountNegBin$new(count_des, verbose = FALSE),
		seed = 209
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountHurdleNegBin$new(hurdle_des, verbose = FALSE),
		SlowInferenceCountHurdleNegBin$new(hurdle_des, verbose = FALSE),
		seed = 210
	)
	expect_false(inherits(InferenceCountHurdleNegBin$new(hurdle_des, verbose = FALSE), "InferenceCountHurdleNegBinAbstract"))
	compare_bootstrap_fast_slow_asymp(
		InferenceCountHurdleNegBin$new(hurdle_des, verbose = FALSE),
		SlowInferenceCountHurdleNegBin$new(hurdle_des, verbose = FALSE),
		seed = 211
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidRiskDiff$new(incid_des, verbose = FALSE),
		seed = 212
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidNewcombeRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidNewcombeRiskDiff$new(incid_des, verbose = FALSE),
		seed = 213
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidMiettinenNurminenRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidMiettinenNurminenRiskDiff$new(incid_des, verbose = FALSE),
		seed = 214
	)
})

test_that("continuous blocked-cluster bootstrap keeps multivariate workers finite", {
	set.seed(20260417)
	block = rep(c("a", "b", "c"), each = 8L)
	X = data.frame(
		block = factor(block),
		x1 = rnorm(length(block)),
		x2 = rnorm(length(block)),
		x3 = rnorm(length(block))
	)
	des = make_fixed_blocked_cluster_design(
		X,
		y_fun = function(w, X_design){
			0.6 * w + 0.4 * X_design$x1 - 0.25 * X_design$x2 + stats::rnorm(nrow(X_design), sd = 0.5)
		}
	)

	inf_ols = InferenceContinOLS$new(des, verbose = FALSE)
	inf_lin = InferenceContinLin$new(des, verbose = FALSE)
	inf_robust = InferenceContinRobustRegr$new(des, method = "M", verbose = FALSE)

	set.seed(215)
	dbg_ols = inf_ols$approximate_bootstrap_distribution_beta_hat_T(B = 25L, show_progress = FALSE, debug = TRUE)
	set.seed(216)
	dbg_lin = inf_lin$approximate_bootstrap_distribution_beta_hat_T(B = 25L, show_progress = FALSE, debug = TRUE)
	set.seed(217)
	dbg_robust = inf_robust$approximate_bootstrap_distribution_beta_hat_T(B = 25L, show_progress = FALSE, debug = TRUE)

	expect_equal(sum(is.finite(dbg_ols$values)), 25L)
	expect_equal(sum(is.finite(dbg_lin$values)), 25L)
	expect_equal(sum(is.finite(dbg_robust$values)), 25L)
	expect_false(any(grepl("incompatible dimensions|number of rows of result", unlist(dbg_ols$errors))))
	expect_false(any(grepl("incompatible dimensions|number of rows of result", unlist(dbg_lin$errors))))
	expect_false(any(grepl("incompatible dimensions|number of rows of result", unlist(dbg_robust$errors))))

	set.seed(218)
	ci_ols = inf_ols$compute_bootstrap_confidence_interval(B = 51L, show_progress = FALSE)
	set.seed(219)
	ci_lin = inf_lin$compute_bootstrap_confidence_interval(B = 51L, show_progress = FALSE)
	set.seed(220)
	ci_robust = inf_robust$compute_bootstrap_confidence_interval(B = 51L, show_progress = FALSE)
	expect_true(all(is.finite(ci_ols)))
	expect_true(all(is.finite(ci_lin)))
	expect_true(all(is.finite(ci_robust)))
})

test_that("MLE and proportion families picked up through InferenceAsymp match generic bootstrap", {
	SlowInferenceOrdinalPropOddsRegr = R6::R6Class(
		"SlowInferenceOrdinalPropOddsRegr",
		inherit = InferenceOrdinalPropOddsRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropBetaRegr = R6::R6Class(
		"SlowInferencePropBetaRegr",
		inherit = InferencePropBetaRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropFractionalLogit = R6::R6Class(
		"SlowInferencePropFractionalLogit",
		inherit = InferencePropFractionalLogit,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropFractionalLogit = R6::R6Class(
		"SlowInferencePropFractionalLogit",
		inherit = InferencePropFractionalLogit,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropZeroOneInflatedBetaRegr = R6::R6Class(
		"SlowInferencePropZeroOneInflatedBetaRegr",
		inherit = InferencePropZeroOneInflatedBetaRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropZeroOneInflatedBetaRegr = R6::R6Class(
		"SlowInferencePropZeroOneInflatedBetaRegr",
		inherit = InferencePropZeroOneInflatedBetaRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260417)
	n = 40L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

	ordinal_des = make_fixed_design(
		"ordinal",
		X,
		y_fun = function(w, X){
			latent = 0.5 * w + 0.3 * X$x1 - 0.2 * X$x2 + stats::rlogis(nrow(X))
			as.integer(cut(latent, breaks = c(-Inf, -0.8, 0.0, 0.8, Inf), labels = FALSE))
		}
	)
	prop_des = make_fixed_design(
		"proportion",
		X,
		y_fun = function(w, X){
			stats::plogis(-0.3 + 0.8 * w + 0.2 * X$x1 + stats::rnorm(nrow(X), sd = 0.35))
		}
	)
	zoib_des = make_fixed_design(
		"proportion",
		X,
		y_fun = function(w, X){
			u = stats::runif(nrow(X))
			mu = stats::plogis(-0.4 + 0.9 * w + 0.2 * X$x1)
			out = stats::rbeta(nrow(X), shape1 = 6 * mu + 0.5, shape2 = 6 * (1 - mu) + 0.5)
			out[u < 0.15] = 0
			out[u > 0.85] = 1
			out
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceOrdinalPropOddsRegr$new(ordinal_des, verbose = FALSE),
		SlowInferenceOrdinalPropOddsRegr$new(ordinal_des, verbose = FALSE),
		seed = 215
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropBetaRegr$new(prop_des, verbose = FALSE),
		SlowInferencePropBetaRegr$new(prop_des, verbose = FALSE),
		seed = 216
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropFractionalLogit$new(prop_des, verbose = FALSE),
		SlowInferencePropFractionalLogit$new(prop_des, verbose = FALSE),
		seed = 217,
		tolerance = 1e-9
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropFractionalLogit$new(prop_des, verbose = FALSE),
		SlowInferencePropFractionalLogit$new(prop_des, verbose = FALSE),
		seed = 218,
		tolerance = 1e-9
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		SlowInferencePropZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		seed = 219
	)
	expect_false(inherits(InferencePropZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE), "InferencePropZeroOneInflatedBetaAbstract"))
	compare_bootstrap_fast_slow_asymp(
		InferencePropZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		SlowInferencePropZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		seed = 220
	)
})

test_that("survival reusable workers match generic bootstrap", {
	SlowInferenceSurvivalCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalCoxPHRegr",
		inherit = InferenceSurvivalCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalCoxPHRegr",
		inherit = InferenceSurvivalCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalStratCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalStratCoxPHRegr",
		inherit = InferenceSurvivalStratCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalStratCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalStratCoxPHRegr",
		inherit = InferenceSurvivalStratCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalKMDiff = R6::R6Class(
		"SlowInferenceSurvivalKMDiff",
		inherit = InferenceSurvivalKMDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalLogRank = R6::R6Class(
		"SlowInferenceSurvivalLogRank",
		inherit = InferenceSurvivalLogRank,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalRestrictedMeanDiff = R6::R6Class(
		"SlowInferenceSurvivalRestrictedMeanDiff",
		inherit = InferenceSurvivalRestrictedMeanDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalGehanWilcox = R6::R6Class(
		"SlowInferenceSurvivalGehanWilcox",
		inherit = InferenceSurvivalGehanWilcox,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260418)
	n = 46L
	X = data.frame(x1 = rnorm(n), x2 = sample(0:1, n, replace = TRUE), x3 = rnorm(n))
	surv_des = make_fixed_design(
		"survival",
		X,
		y_fun = function(w, X){
			base_rate = exp(-0.4 + 0.5 * w + 0.2 * X$x1 - 0.3 * X$x2)
			pmax(stats::rexp(nrow(X), rate = base_rate), 0.05)
		},
		dead_fun = function(w, X, y){
			stats::rbinom(length(y), 1, stats::plogis(1.0 - 0.12 * y + 0.2 * w))
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 221
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 222
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalStratCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalStratCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 223
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalStratCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalStratCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 224
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalKMDiff$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalKMDiff$new(surv_des, verbose = FALSE),
		seed = 225
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalLogRank$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalLogRank$new(surv_des, verbose = FALSE),
		seed = 226
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalRestrictedMeanDiff$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalRestrictedMeanDiff$new(surv_des, verbose = FALSE),
		seed = 227
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalGehanWilcox$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalGehanWilcox$new(surv_des, verbose = FALSE),
		seed = 228
	)
})

test_that("zero-augmented count reusable workers match generic bootstrap", {
	skip_if_not_installed("glmmTMB")

	SlowInferenceCountZeroInflatedPoisson = R6::R6Class(
		"SlowInferenceCountZeroInflatedPoisson",
		inherit = InferenceCountZeroInflatedPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountZeroInflatedPoisson = R6::R6Class(
		"SlowInferenceCountZeroInflatedPoisson",
		inherit = InferenceCountZeroInflatedPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountZeroInflatedNegBin = R6::R6Class(
		"SlowInferenceCountZeroInflatedNegBin",
		inherit = InferenceCountZeroInflatedNegBin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountZeroInflatedNegBin = R6::R6Class(
		"SlowInferenceCountZeroInflatedNegBin",
		inherit = InferenceCountZeroInflatedNegBin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountHurdlePoisson = R6::R6Class(
		"SlowInferenceCountHurdlePoisson",
		inherit = InferenceCountHurdlePoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountHurdlePoisson = R6::R6Class(
		"SlowInferenceCountHurdlePoisson",
		inherit = InferenceCountHurdlePoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260419)
	n = 36L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = make_fixed_design(
		"count",
		X,
		y_fun = function(w, X){
			is_zero = stats::rbinom(nrow(X), 1, stats::plogis(0.4 - 0.9 * w + 0.2 * X$x1))
			pos = stats::rpois(nrow(X), lambda = exp(0.3 + 0.4 * w + 0.2 * X$x2))
			ifelse(is_zero == 1, 0, pos)
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceCountZeroInflatedPoisson$new(des, verbose = FALSE),
		SlowInferenceCountZeroInflatedPoisson$new(des, verbose = FALSE),
		seed = 229,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountZeroInflatedPoisson$new(des, verbose = FALSE),
		SlowInferenceCountZeroInflatedPoisson$new(des, verbose = FALSE),
		seed = 230,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountZeroInflatedNegBin$new(des, verbose = FALSE),
		SlowInferenceCountZeroInflatedNegBin$new(des, verbose = FALSE),
		seed = 231,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountZeroInflatedNegBin$new(des, verbose = FALSE),
		SlowInferenceCountZeroInflatedNegBin$new(des, verbose = FALSE),
		seed = 232,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountHurdlePoisson$new(des, verbose = FALSE),
		SlowInferenceCountHurdlePoisson$new(des, verbose = FALSE),
		seed = 233,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountHurdlePoisson$new(des, verbose = FALSE),
		SlowInferenceCountHurdlePoisson$new(des, verbose = FALSE),
		seed = 234,
		tolerance = 1e-8
	)
})
