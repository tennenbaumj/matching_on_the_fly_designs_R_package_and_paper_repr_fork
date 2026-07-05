#!/usr/bin/env Rscript
# EDI Kernel Profiler Dispatcher
# ===
# Usage:  Rscript edi_kernel_profiler.R <kernel_name>
# Seed:   42 (all kernels)
# N:      1000 GLM / 500 survival+wilcox (estimate-only)  |  200 (variance)
#
# Estimate-only kernels (point-estimate benchmark paths):
#   logistic_est, ols_est, poisson_est, negbin_est, beta_est, robust_est,
#   logbin_est, probit_est, identbin_est, hurdle_p_est, zip_est, zinb_est,
#   hurdle_nb_est, prop_odds_est, adj_cat_est, cont_ratio_est, ord_probit_est,
#   ord_cloglog_est, ord_cauchit_est, coxph_est, strat_coxph_est, weibull_est,
#   logrank_est, km_diff_est, rmean_diff_est, wilcox_est,
#   gcomp_logistic_rd_est, gcomp_logistic_rr_est, gcomp_frac_logit_est, gcomp_ordinal_est
#
# Variance/full-inference kernels (Wald benchmark paths):
#   ols_var, ols_lin_var, robust_var, logistic_var, logbin_var, probit_var,
#   newcombe_var, mn_var, poisson_var, quasi_var, negbin_var, hurdle_nb_var,
#   hurdle_p_var, zip_var, zinb_var, poisson_robust_var, beta_var, prop_odds_var,
#   adj_cat_var, cont_ratio_var, ridit_var, jt_var, coxph_var, strat_coxph_var,
#   weibull_var, logrank_var, wilcox_var, gcomp_logistic_post_fit_var, gcomp_ordinal_var

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
    cat("Usage: Rscript edi_kernel_profiler.R <kernel_name>\n")
    quit(status = 1L)
}
KERNEL = args[1L]

.libPaths(c(
    file.path(Sys.getenv("HOME"), "R",
              paste0(R.version$platform, "-library"),
              paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")),
    .libPaths()
))
suppressPackageStartupMessages(library(EDI))
SEED = 42L

# ---------------------------------------------------------------------------
# Data generators
# ---------------------------------------------------------------------------
make_data = function(n, family, seed = SEED) {
    set.seed(seed)
    p = 5L
    X = matrix(rnorm(n * p), n, p); X[, 1L] = 1
    beta = rnorm(p) * 0.5
    n_treat = floor(n / 2L)
    w = sample(c(rep(1L, n_treat), rep(0L, n - n_treat)))
    eta = as.numeric(X %*% beta + 0.5 * w)

    y = switch(family,
        continuous    = as.numeric(eta + rnorm(n, 0, 0.5)),
        logistic      = as.numeric(rbinom(n, 1L, plogis(eta))),
        `log-binomial`= as.numeric(rbinom(n, 1L, pmin(0.5, exp(eta - 2)))),
        poisson       = as.numeric(rpois(n, exp(eta))),
        negbin        = as.numeric(rnbinom(n, size = 2, mu = exp(eta))),
        beta          = { mu = plogis(eta); phi = 10
                          pmax(pmin(rbeta(n, mu * phi, (1 - mu) * phi), 1 - 1e-6), 1e-6) },
        cox           = rexp(n, exp(eta)),
        ordinal       = {
            p1 = 1 / (1 + exp(eta - 1)); p2 = 1 / (1 + exp(eta + 1)) - p1; p3 = 1 - p1 - p2
            pr = pmax(cbind(p1, p2, p3), 1e-6); pr = pr / rowSums(pr)
            as.numeric(apply(pr, 1L, function(p) sample(1:3, 1L, prob = p)))
        },
        stop("Unknown family: ", family)
    )

    dead = if (family == "cox") as.integer(rbinom(n, 1L, 0.8)) else NULL
    X_cov = X[, -1L, drop = FALSE]; colnames(X_cov) = paste0("x", 1:4)
    X_bm  = cbind(`(Intercept)` = 1, treatment = w, X_cov)
    X_ord = cbind(treatment = w, X_cov)
    # IPW-style bootstrap weights: positive, mean 1, sum to n
    wts_raw = rexp(n); wts_bm = as.numeric(wts_raw / mean(wts_raw))
    list(X = X, X_bm = X_bm, X_ord = X_ord, y_bm = y,
         w_bm = as.integer(w), dead_bm = dead, wts_bm = wts_bm, n = n)
}

make_strat_cox_cache = function(n = 500L, seed = SEED) {
    set.seed(seed)
    p = 5L
    X = matrix(rnorm(n * p), n, p); X[, 1L] = 1
    beta = rnorm(p) * 0.5
    n_treat = floor(n / 2L)
    w = sample(c(rep(1L, n_treat), rep(0L, n - n_treat)))
    # inject low-cardinality strata
    strata_grid = as.matrix(expand.grid(x1 = 0:1, x2 = 0:2))
    idx = sample(rep(seq_len(nrow(strata_grid)), length.out = n))
    X[, 2L] = strata_grid[idx, 1L]; X[, 3L] = strata_grid[idx, 2L]
    beta_cov = c(0.45, -0.35, 0.20, -0.15)
    eta = 0.5 * w + drop(X[, 2:5, drop = FALSE] %*% beta_cov)
    y   = rexp(n, exp(eta))
    dead = as.integer(rbinom(n, 1L, 0.8))
    X_cov = X[, -1L, drop = FALSE]; colnames(X_cov) = paste0("x", 1:4)

    si    = EDI:::compute_survival_strata_ids_cpp(X_cov)
    si_id = as.integer(si$strata_id)
    inf   = integer(0)
    for (s in unique(si_id)) {
        i_s = which(si_id == s)
        if (length(i_s) < 2L || length(unique(w[i_s])) < 2L) next
        if (!any(dead[i_s] == 1L, na.rm = TRUE)) next
        inf = c(inf, i_s)
    }
    inf = sort(unique(inf))
    if (length(inf) >= 4L) {
        X_s = matrix(w[inf], ncol = 1L, dimnames = list(NULL, "w"))
        cache = build_stratified_cox_data_cache_cpp(X_s, y[inf], as.numeric(dead[inf]), si_id[inf])
    } else {
        X_s = matrix(w, ncol = 1L, dimnames = list(NULL, "w"))
        cache = build_cox_data_cache_cpp(X_s, y, as.numeric(dead))
    }
    list(cache = cache, w_bm = as.integer(w), y_bm = y, dead_bm = dead)
}

make_glmm_data = function(n = 400L, family = "logistic", seed = SEED, n_clusters = 80L) {
    d = make_data(n, family, seed)
    group_id     = as.integer(rep(seq_len(n_clusters), length.out = n))
    pairs_group_id = as.integer(rep(seq_len(n %/% 2L), each = 2L))
    c(d, list(group_id_bm = group_id, pairs_group_id_bm = pairs_group_id))
}

# ---------------------------------------------------------------------------
# Kernel table
# ---------------------------------------------------------------------------
# Each entry: list(desc, family, n, REPS, setup_fn, kernel_fn)
#   setup_fn()  -> sets up data, returns env vars as a list
#   kernel_fn() -> the timed call (closure over setup env)

get_kernel = function(name) {
    switch(name,

    # ================================================================
    # ESTIMATE-ONLY PATHS  (point-estimate benchmark, N=1000/500)
    # ================================================================

    logistic_est = local({
        # fast_logistic_regression_cpp(X_bm, y_bm, estimate_only=TRUE)
        # Classes: InferenceIncidLogRegr, InferencePropFractionalLogit
        d = make_data(1000L, "logistic")
        list(desc = "fast_logistic_regression_cpp(X_bm, y_bm, estimate_only=TRUE)",
             REPS = 115000L,
             fn   = function() fast_logistic_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE))
    }),

    ols_est = local({
        # fast_ols_cpp(X_bm, y_bm)
        # Classes: InferenceContinOLS, InferenceIncidRiskDiff
        d = make_data(1000L, "continuous")
        list(desc = "fast_ols_cpp(X_bm, y_bm)",
             REPS = 500000L,
             fn   = function() fast_ols_cpp(d$X_bm, d$y_bm))
    }),

    poisson_est = local({
        # fast_poisson_regression_cpp(X_bm, y_bm, estimate_only=TRUE, optimization_alg="irls")
        # Classes: InferenceCountPoisson, InferenceCountQuasiPoisson,
        #          InferenceCountRobustPoisson, InferenceIncidModifiedPoisson
        d = make_data(1000L, "poisson")
        list(desc = "fast_poisson_regression_cpp(X_bm, y_bm, estimate_only=TRUE, optimization_alg='irls')",
             REPS = 65000L,
             fn   = function() fast_poisson_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE, optimization_alg = "irls"))
    }),

    negbin_est = local({
        # fast_neg_bin_cpp(X_bm, as.integer(y_bm), estimate_only=TRUE)
        # Class: InferenceCountNegBin  (benchmark uses poisson data for point-est)
        d = make_data(1000L, "poisson")
        list(desc = "fast_neg_bin_cpp(X_bm, as.integer(y_bm), estimate_only=TRUE)",
             REPS = 33000L,
             fn   = function() fast_neg_bin_cpp(d$X_bm, as.integer(d$y_bm), estimate_only = TRUE))
    }),

    beta_est = local({
        # fast_beta_regression_cpp(X_bm, y_bm, estimate_only=TRUE)
        # Class: InferencePropBetaRegr
        d = make_data(1000L, "beta")
        list(desc = "fast_beta_regression_cpp(X_bm, y_bm, estimate_only=TRUE)",
             REPS = 9000L,
             fn   = function() fast_beta_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE))
    }),

    robust_est = local({
        # fast_robust_regression_cpp(X_bm, y_bm, estimate_only=TRUE)
        # Class: InferenceContinRobustRegr
        d = make_data(1000L, "continuous")
        list(desc = "fast_robust_regression_cpp(X_bm, y_bm, estimate_only=TRUE)",
             REPS = 44000L,
             fn   = function() fast_robust_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE))
    }),

    logbin_est = local({
        # fast_log_binomial_regression_cpp(X_bm, y_bm)
        # Class: InferenceIncidLogBinomial
        d = make_data(1000L, "log-binomial")
        list(desc = "fast_log_binomial_regression_cpp(X_bm, y_bm)",
             REPS = 10000L,
             fn   = function() fast_log_binomial_regression_cpp(d$X_bm, d$y_bm))
    }),

    probit_est = local({
        # fast_probit_regression_cpp(X_bm, y_bm, estimate_only=TRUE)
        # Class: InferenceIncidProbitRegr
        d = make_data(1000L, "logistic")
        list(desc = "fast_probit_regression_cpp(X_bm, y_bm, estimate_only=TRUE)",
             REPS = 28000L,
             fn   = function() fast_probit_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE))
    }),

    identbin_est = local({
        # fast_identity_binomial_regression_cpp(X_bm, y_bm)
        # Class: InferenceIncidBinomialIdentityRiskDiff
        d = make_data(1000L, "logistic")
        list(desc = "fast_identity_binomial_regression_cpp(X_bm, y_bm)",
             REPS = 125000L,
             fn   = function() fast_identity_binomial_regression_cpp(d$X_bm, d$y_bm))
    }),

    hurdle_p_est = local({
        # fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=TRUE, estimate_only=TRUE)
        # Class: InferenceCountHurdlePoisson
        d = make_data(1000L, "poisson")
        list(desc = "fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=TRUE, estimate_only=TRUE)",
             REPS = 8000L,
             fn   = function() fast_zero_augmented_poisson_cpp(d$X_bm, d$y_bm, d$X_bm, is_hurdle = TRUE, estimate_only = TRUE))
    }),

    zip_est = local({
        # fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=FALSE, estimate_only=TRUE)
        # Class: InferenceCountZeroInflatedPoisson
        d = make_data(1000L, "poisson")
        list(desc = "fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=FALSE, estimate_only=TRUE)",
             REPS = 3800L,
             fn   = function() fast_zero_augmented_poisson_cpp(d$X_bm, d$y_bm, d$X_bm, is_hurdle = FALSE, estimate_only = TRUE))
    }),

    zinb_est = local({
        # fast_zinb_cpp(X_bm, X_bm, y_bm, estimate_only=TRUE)
        # Class: InferenceCountZeroInflatedNegBin  (benchmark uses poisson data for point-est)
        d = make_data(1000L, "poisson")
        list(desc = "fast_zinb_cpp(X_bm, X_bm, y_bm, estimate_only=TRUE)",
             REPS = 8700L,
             fn   = function() fast_zinb_cpp(d$X_bm, d$X_bm, d$y_bm, estimate_only = TRUE))
    }),

    hurdle_nb_est = local({
        # EDI:::fast_hurdle_negbin_cpp(X_bm, as.integer(y_bm), X_bm, estimate_only=TRUE)
        # Class: InferenceCountHurdleNegBin  (benchmark uses poisson data for point-est)
        d = make_data(1000L, "poisson")
        list(desc = "EDI:::fast_hurdle_negbin_cpp(X_bm, as.integer(y_bm), X_bm, estimate_only=TRUE)",
             REPS = 8700L,
             fn   = function() EDI:::fast_hurdle_negbin_cpp(d$X_bm, as.integer(d$y_bm), d$X_bm, estimate_only = TRUE))
    }),

    prop_odds_est = local({
        # fast_ordinal_regression_cpp(X_ord, y_bm, estimate_only=TRUE)
        # Classes: InferenceOrdinalPropOddsRegr, InferenceOrdinalGCompMeanDiff (fit part)
        d = make_data(1000L, "ordinal")
        list(desc = "fast_ordinal_regression_cpp(X_ord, y_bm, estimate_only=TRUE)",
             REPS = 17000L,
             fn   = function() fast_ordinal_regression_cpp(d$X_ord, d$y_bm, estimate_only = TRUE))
    }),

    adj_cat_est = local({
        # fast_adjacent_category_logit_cpp(X_ord, y_bm)
        # Class: InferenceOrdinalAdjCatLogitRegr
        d = make_data(1000L, "ordinal")
        list(desc = "fast_adjacent_category_logit_cpp(X_ord, y_bm)",
             REPS = 29000L,
             fn   = function() fast_adjacent_category_logit_cpp(d$X_ord, d$y_bm))
    }),

    cont_ratio_est = local({
        # fast_continuation_ratio_regression_cpp(X_ord, y_bm)
        # Class: InferenceOrdinalContRatioRegr
        d = make_data(1000L, "ordinal")
        list(desc = "fast_continuation_ratio_regression_cpp(X_ord, y_bm)",
             REPS = 35000L,
             fn   = function() fast_continuation_ratio_regression_cpp(d$X_ord, d$y_bm))
    }),

    ord_probit_est = local({
        # fast_ordinal_probit_regression_cpp(X_ord, y_bm, estimate_only=TRUE)
        # Class: InferenceOrdinalOrderedProbitRegr
        d = make_data(1000L, "ordinal")
        list(desc = "fast_ordinal_probit_regression_cpp(X_ord, y_bm, estimate_only=TRUE)",
             REPS = 21000L,
             fn   = function() fast_ordinal_probit_regression_cpp(d$X_ord, d$y_bm, estimate_only = TRUE))
    }),

    ord_cloglog_est = local({
        # fast_ordinal_cloglog_regression_cpp(X_ord, y_bm, estimate_only=TRUE)
        # Class: InferenceOrdinalCloglogRegr
        d = make_data(1000L, "ordinal")
        list(desc = "fast_ordinal_cloglog_regression_cpp(X_ord, y_bm, estimate_only=TRUE)",
             REPS = 24000L,
             fn   = function() fast_ordinal_cloglog_regression_cpp(d$X_ord, d$y_bm, estimate_only = TRUE))
    }),

    ord_cauchit_est = local({
        # fast_ordinal_cauchit_regression_cpp(X_ord, y_bm, estimate_only=TRUE)
        # Class: InferenceOrdinalCauchitRegr
        d = make_data(1000L, "ordinal")
        list(desc = "fast_ordinal_cauchit_regression_cpp(X_ord, y_bm, estimate_only=TRUE)",
             REPS = 21000L,
             fn   = function() fast_ordinal_cauchit_regression_cpp(d$X_ord, d$y_bm, estimate_only = TRUE))
    }),

    coxph_est = local({
        # fast_coxph_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=TRUE)
        # Class: InferenceSurvivalCoxPHRegr
        d = make_data(500L, "cox")
        list(desc = "fast_coxph_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=TRUE)",
             REPS = 43000L,
             fn   = function() fast_coxph_regression_cpp(d$X_ord, d$y_bm, d$dead_bm, estimate_only = TRUE))
    }),

    strat_coxph_est = local({
        # fast_coxph_regression_prebuilt_cpp(cache, estimate_only=TRUE)
        # Class: InferenceSurvivalStratCoxPHRegr
        sc = make_strat_cox_cache(500L)
        list(desc = "fast_coxph_regression_prebuilt_cpp(cache, estimate_only=TRUE) [stratified]",
             REPS = 25000L,
             fn   = function() fast_coxph_regression_prebuilt_cpp(sc$cache, estimate_only = TRUE))
    }),

    weibull_est = local({
        # fast_weibull_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=TRUE)
        # Class: InferenceSurvivalWeibullRegr
        d = make_data(500L, "cox")
        list(desc = "fast_weibull_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=TRUE)",
             REPS = 150000L,
             fn   = function() fast_weibull_regression_cpp(d$X_ord, d$y_bm, d$dead_bm, estimate_only = TRUE))
    }),

    logrank_est = local({
        # EDI:::fast_logrank_stats_cpp(w_bm, y_bm, dead_bm)
        # Class: InferenceSurvivalLogRank
        d = make_data(500L, "cox")
        list(desc = "EDI:::fast_logrank_stats_cpp(w_bm, y_bm, dead_bm)",
             REPS = 750000L,
             fn   = function() EDI:::fast_logrank_stats_cpp(d$w_bm, d$y_bm, d$dead_bm))
    }),

    km_diff_est = local({
        # EDI:::get_survival_stat_diff(y_bm, dead_bm, w_bm, "median")
        # Class: InferenceSurvivalKMDiff
        d = make_data(500L, "cox")
        list(desc = "EDI:::get_survival_stat_diff(y_bm, dead_bm, w_bm, 'median')",
             REPS = 750000L,
             fn   = function() EDI:::get_survival_stat_diff(d$y_bm, d$dead_bm, d$w_bm, "median"))
    }),

    rmean_diff_est = local({
        # EDI:::get_survival_stat_diff(y_bm, dead_bm, w_bm, "restricted_mean")
        # Class: InferenceSurvivalRestrictedMeanDiff
        d = make_data(500L, "cox")
        list(desc = "EDI:::get_survival_stat_diff(y_bm, dead_bm, w_bm, 'restricted_mean')",
             REPS = 750000L,
             fn   = function() EDI:::get_survival_stat_diff(d$y_bm, d$dead_bm, d$w_bm, "restricted_mean"))
    }),

    wilcox_est = local({
        # EDI:::wilcox_hl_point_estimate_cpp(w_bm, y_bm)
        # Class: InferenceAllSimpleWilcox  (scale=0.5 -> N=500)
        d = make_data(500L, "continuous")
        list(desc = "EDI:::wilcox_hl_point_estimate_cpp(w_bm, y_bm)",
             REPS = 25000L,
             fn   = function() EDI:::wilcox_hl_point_estimate_cpp(d$w_bm, d$y_bm))
    }),

    gcomp_logistic_rd_est = local({
        # fast_logistic_regression_cpp(X_bm,y_bm,estimate_only=TRUE) + gcomp_logistic_point_estimate_cpp
        # Class: InferenceIncidGCompRiskDiff
        d = make_data(1000L, "logistic")
        list(desc = "fast_logistic_regression_cpp + EDI:::gcomp_logistic_point_estimate_cpp(rd)",
             REPS = 83000L,
             fn   = function() {
                 b = fast_logistic_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE)$b
                 EDI:::gcomp_logistic_point_estimate_cpp(d$X_bm, b, 2L)$md
             })
    }),

    gcomp_logistic_rr_est = local({
        # fast_logistic_regression_cpp + gcomp_logistic_point_estimate_cpp (rr)
        # Class: InferenceIncidGCompRiskRatio
        d = make_data(1000L, "logistic")
        list(desc = "fast_logistic_regression_cpp + EDI:::gcomp_logistic_point_estimate_cpp(rr)",
             REPS = 75000L,
             fn   = function() {
                 b = fast_logistic_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE)$b
                 res = EDI:::gcomp_logistic_point_estimate_cpp(d$X_bm, b, 2L)
                 res$mean1 / res$mean0
             })
    }),

    gcomp_frac_logit_est = local({
        # fast_logistic_regression_cpp + gcomp_fractional_logit_point_estimate_cpp
        # Class: InferencePropGCompMeanDiff (estimate-only path)
        d = make_data(1000L, "beta")
        list(desc = "fast_logistic_regression_cpp + EDI:::gcomp_fractional_logit_point_estimate_cpp",
             REPS = 100000L,
             fn   = function() {
                 b = fast_logistic_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE)$b
                 EDI:::gcomp_fractional_logit_point_estimate_cpp(d$X_bm, b, 2L)$md
             })
    }),

    gcomp_ordinal_est = local({
        # fast_ordinal_regression_cpp + gcomp_ordinal_proportional_odds_post_fit_cpp
        # Class: InferenceOrdinalGCompMeanDiff (estimate-only path)
        d = make_data(1000L, "ordinal")
        list(desc = "fast_ordinal_regression_cpp + gcomp_ordinal_proportional_odds_post_fit_cpp",
             REPS = 16000L,
             fn   = function() {
                 fit = fast_ordinal_regression_cpp(d$X_ord, d$y_bm, estimate_only = TRUE)
                 gcomp_ordinal_proportional_odds_post_fit_cpp(
                     X_fit = d$X_ord, coef_hat = as.numeric(fit$b),
                     alpha_hat = as.numeric(fit$alpha), j_treat = 1L)$md
             })
    }),

    # ================================================================
    # VARIANCE / FULL-INFERENCE PATHS  (Wald benchmark, N=200)
    # ================================================================

    ols_var = local({
        # fast_ols_with_var_cpp(X_bm, y_bm, j=2L)
        # Class: InferenceContinOLS (Wald), InferenceIncidRiskDiff (Wald)
        d = make_data(200L, "continuous")
        list(desc = "fast_ols_with_var_cpp(X_bm, y_bm, j=2L)",
             REPS = 1500000L,
             fn   = function() {
                 res = EDI:::fast_ols_with_var_cpp(d$X_bm, d$y_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    ols_lin_var = local({
        # fast_ols_cpp + ols_hc2_post_fit_cpp  (InferenceContinLin Wald)
        d = make_data(200L, "continuous")
        Xc = scale(d$X_bm[, -1L, drop = FALSE], scale = FALSE)
        X_lin = cbind(1, d$w_bm, Xc, Xc * d$w_bm)
        colnames(X_lin)[1:2] = c("(Intercept)", "treatment")
        list(desc = "fast_ols_cpp(X_lin) + EDI:::ols_hc2_post_fit_cpp(X_lin, y_bm, b, 2L)",
             REPS = 83000L,
             fn   = function() {
                 res_fit = fast_ols_cpp(X_lin, d$y_bm)
                 post_fit = EDI:::ols_hc2_post_fit_cpp(X_lin, d$y_bm, as.numeric(res_fit$b), 2L)
                 sqrt(post_fit$ssq_hat)
             })
    }),

    robust_var = local({
        # fast_robust_regression_cpp(X_bm, y_bm, j=2L)
        # Class: InferenceContinRobustRegr (Wald)
        d = make_data(200L, "continuous")
        list(desc = "fast_robust_regression_cpp(X_bm, y_bm, j=2L)",
             REPS = 188000L,
             fn   = function() {
                 res = fast_robust_regression_cpp(d$X_bm, d$y_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    logistic_var = local({
        # fast_logistic_regression_with_var_cpp(X_bm, y_bm, j=2L)
        # Class: InferenceIncidLogRegr (Wald)
        d = make_data(200L, "logistic")
        list(desc = "fast_logistic_regression_with_var_cpp(X_bm, y_bm, j=2L)",
             REPS = 375000L,
             fn   = function() {
                 res = EDI:::fast_logistic_regression_with_var_cpp(d$X_bm, d$y_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    logbin_var = local({
        # fast_log_binomial_regression_with_var_cpp(X_bm, y_bm, j=2L)
        # Class: InferenceIncidLogBinomial (Wald)
        d = make_data(200L, "log-binomial")
        list(desc = "fast_log_binomial_regression_with_var_cpp(X_bm, y_bm, j=2L)",
             REPS = 13000L,
             fn   = function() {
                 res = EDI:::fast_log_binomial_regression_with_var_cpp(d$X_bm, d$y_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    probit_var = local({
        # fast_probit_regression_with_var_cpp(X_bm, y_bm, j=2L)
        # Class: InferenceIncidProbitRegr (Wald)
        d = make_data(200L, "logistic")
        list(desc = "fast_probit_regression_with_var_cpp(X_bm, y_bm, j=2L)",
             REPS = 107000L,
             fn   = function() {
                 res = EDI:::fast_probit_regression_with_var_cpp(d$X_bm, d$y_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    newcombe_var = local({
        # newcombe_independent_ci_cpp(x_t, n_t, x_c, n_c, alpha)
        # Class: InferenceIncidNewcombeRiskDiff (Wald)
        d = make_data(200L, "logistic")
        x_t = sum(d$y_bm[d$w_bm == 1L]); n_t = sum(d$w_bm == 1L)
        x_c = sum(d$y_bm[d$w_bm == 0L]); n_c = sum(d$w_bm == 0L)
        list(desc = "newcombe_independent_ci_cpp(x_t, n_t, x_c, n_c, alpha)",
             REPS = 188000L,
             fn   = function() newcombe_independent_ci_cpp(x_t, n_t, x_c, n_c, 0.05))
    }),

    mn_var = local({
        # mn_ci_cpp(x_t, n_t, x_c, n_c, p_t, p_c, alpha, pval_epsilon)
        # Class: InferenceIncidMiettinenNurminenRiskDiff (Wald)
        d = make_data(200L, "logistic")
        x_t = sum(d$y_bm[d$w_bm == 1L]); n_t = sum(d$w_bm == 1L)
        x_c = sum(d$y_bm[d$w_bm == 0L]); n_c = sum(d$w_bm == 0L)
        list(desc = "mn_ci_cpp(x_t, n_t, x_c, n_c, x_t/n_t, x_c/n_c, 0.05, 1e-7)",
             REPS = 500000L,
             fn   = function() mn_ci_cpp(x_t, n_t, x_c, n_c, x_t/n_t, x_c/n_c, 0.05, 1e-7))
    }),

    poisson_var = local({
        # fast_poisson_regression_with_var_cpp(X_bm, y_bm, j=2L)
        # Class: InferenceCountPoisson (Wald)
        d = make_data(200L, "poisson")
        list(desc = "fast_poisson_regression_with_var_cpp(X_bm, y_bm, j=2L)",
             REPS = 214000L,
             fn   = function() {
                 res = EDI:::fast_poisson_regression_with_var_cpp(d$X_bm, d$y_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    quasi_var = local({
        # fast_quasipoisson_regression_with_var_cpp(X_bm, y_bm, j=2L)
        # Class: InferenceCountQuasiPoisson (Wald)
        d = make_data(200L, "poisson")
        list(desc = "fast_quasipoisson_regression_with_var_cpp(X_bm, y_bm, j=2L)",
             REPS = 188000L,
             fn   = function() {
                 res = EDI:::fast_quasipoisson_regression_with_var_cpp(d$X_bm, d$y_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    negbin_var = local({
        # fast_neg_bin_with_var_cpp(X_bm, as.integer(y_bm))
        # Class: InferenceCountNegBin (Wald)  -- uses negbin data for Wald
        d = make_data(200L, "negbin")
        list(desc = "fast_neg_bin_with_var_cpp(X_bm, as.integer(y_bm))",
             REPS = 54000L,
             fn   = function() {
                 res = fast_neg_bin_with_var_cpp(d$X_bm, as.integer(d$y_bm))
                 vcov = solve(res$hess_fisher_info_matrix); sqrt(vcov[2L, 2L])
             })
    }),

    hurdle_nb_var = local({
        # fast_hurdle_negbin_with_var_cpp(X_bm, as.integer(y_bm), X_bm, j=2L)
        # Class: InferenceCountHurdleNegBin (Wald) -- uses negbin data for Wald
        d = make_data(200L, "negbin")
        list(desc = "fast_hurdle_negbin_with_var_cpp(X_bm, as.integer(y_bm), X_bm, j=2L)",
             REPS = 58000L,
             fn   = function() {
                 res = EDI:::fast_hurdle_negbin_with_var_cpp(d$X_bm, as.integer(d$y_bm), d$X_bm, j = 2L)
                 sqrt(res$ssq_b_j)
             })
    }),

    hurdle_p_var = local({
        # fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=TRUE, estimate_only=FALSE)
        # Class: InferenceCountHurdlePoisson (Wald)
        d = make_data(200L, "poisson")
        list(desc = "fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=TRUE, estimate_only=FALSE)",
             REPS = 68000L,
             fn   = function() {
                 res = fast_zero_augmented_poisson_cpp(d$X_bm, d$y_bm, d$X_bm, is_hurdle = TRUE, estimate_only = FALSE)
                 sqrt(res$vcov[2L, 2L])
             })
    }),

    zip_var = local({
        # fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=FALSE, estimate_only=FALSE)
        # Class: InferenceCountZeroInflatedPoisson (Wald)
        d = make_data(200L, "poisson")
        list(desc = "fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle=FALSE, estimate_only=FALSE)",
             REPS = 6100L,
             fn   = function() {
                 res = fast_zero_augmented_poisson_cpp(d$X_bm, d$y_bm, d$X_bm, is_hurdle = FALSE, estimate_only = FALSE)
                 sqrt(res$vcov[2L, 2L])
             })
    }),

    zinb_var = local({
        # fast_zinb_cpp(X_bm, X_bm, y_bm, estimate_only=FALSE)
        # Class: InferenceCountZeroInflatedNegBin (Wald) -- uses negbin data for Wald
        d = make_data(200L, "negbin")
        list(desc = "fast_zinb_cpp(X_bm, X_bm, y_bm, estimate_only=FALSE)",
             REPS = 10300L,
             fn   = function() {
                 res = fast_zinb_cpp(d$X_bm, d$X_bm, d$y_bm, estimate_only = FALSE)
                 sqrt(res$vcov[2L, 2L])
             })
    }),

    poisson_robust_var = local({
        # fast_poisson_regression_cpp(X_bm, y_bm, estimate_only=FALSE) for sandwich SE
        # Class: InferenceCountRobustPoisson (Wald)
        d = make_data(200L, "poisson")
        list(desc = "fast_poisson_regression_cpp(X_bm, y_bm, estimate_only=FALSE) + sandwich SE",
             REPS = 125000L,
             fn   = function() {
                 res = fast_poisson_regression_cpp(d$X_bm, d$y_bm, estimate_only = FALSE)
                 bread = solve(res$XtWX)
                 resid = d$y_bm - as.numeric(res$mu)
                 meat  = crossprod(d$X_bm, d$X_bm * (resid^2))
                 sqrt((bread %*% meat %*% bread)[2L, 2L])
             })
    }),

    beta_var = local({
        # fast_beta_regression_with_var_cpp(X_bm, y_bm)
        # Class: InferencePropBetaRegr (Wald)
        d = make_data(200L, "beta")
        list(desc = "fast_beta_regression_with_var_cpp(X_bm, y_bm)",
             REPS = 33000L,
             fn   = function() {
                 res = fast_beta_regression_with_var_cpp(d$X_bm, d$y_bm)
                 sqrt(res$vcov[2L, 2L])
             })
    }),

    prop_odds_var = local({
        # fast_ordinal_regression_with_var_cpp(X_ord, y_bm)
        # Class: InferenceOrdinalPropOddsRegr (Wald)
        d = make_data(200L, "ordinal")
        list(desc = "fast_ordinal_regression_with_var_cpp(X_ord, y_bm)",
             REPS = 63000L,
             fn   = function() {
                 res = fast_ordinal_regression_with_var_cpp(d$X_ord, d$y_bm)
                 sqrt(res$ssq_b_j)
             })
    }),

    adj_cat_var = local({
        # fast_adjacent_category_logit_with_var_cpp(X_ord, y_bm)
        # Class: InferenceOrdinalAdjCatLogitRegr (Wald)
        d = make_data(200L, "ordinal")
        list(desc = "fast_adjacent_category_logit_with_var_cpp(X_ord, y_bm)",
             REPS = 100000L,
             fn   = function() {
                 res = fast_adjacent_category_logit_with_var_cpp(d$X_ord, d$y_bm)
                 sqrt(res$ssq_b_j)
             })
    }),

    cont_ratio_var = local({
        # fast_continuation_ratio_regression_with_var_cpp(X_ord, y_bm)
        # Class: InferenceOrdinalContRatioRegr (Wald)
        d = make_data(200L, "ordinal")
        list(desc = "fast_continuation_ratio_regression_with_var_cpp(X_ord, y_bm)",
             REPS = 136000L,
             fn   = function() {
                 res = fast_continuation_ratio_regression_with_var_cpp(d$X_ord, d$y_bm)
                 sqrt(res$ssq_b_j)
             })
    }),

    ridit_var = local({
        # EDI:::fast_ridit_analysis_cpp(w_bm, as.integer(y_bm), reference="control")
        # Class: InferenceOrdinalRidit (Wald)
        d = make_data(200L, "ordinal")
        list(desc = "EDI:::fast_ridit_analysis_cpp(w_bm, as.integer(y_bm), reference='control')",
             REPS = 1500000L,
             fn   = function() {
                 res = EDI:::fast_ridit_analysis_cpp(d$w_bm, as.integer(d$y_bm), reference = "control")
                 res$estimate / res$se
             })
    }),

    jt_var = local({
        # exact_jonckheere_terpstra_pval_cpp(as.integer(y_bm), w_bm)
        # Class: InferenceOrdinalJonckheereTerpstraTest (Wald)
        d = make_data(200L, "ordinal")
        list(desc = "exact_jonckheere_terpstra_pval_cpp(as.integer(y_bm), w_bm)",
             REPS = 500000L,
             fn   = function() exact_jonckheere_terpstra_pval_cpp(as.integer(d$y_bm), d$w_bm)$p_exact)
    }),

    coxph_var = local({
        # fast_coxph_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=FALSE)
        # Class: InferenceSurvivalCoxPHRegr (Wald)
        d = make_data(200L, "cox")
        list(desc = "fast_coxph_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=FALSE)",
             REPS = 125000L,
             fn   = function() {
                 res = fast_coxph_regression_cpp(d$X_ord, d$y_bm, d$dead_bm, estimate_only = FALSE)
                 sqrt(res$vcov[1L, 1L])
             })
    }),

    strat_coxph_var = local({
        # fast_coxph_regression_prebuilt_cpp(cache, estimate_only=FALSE)
        # Class: InferenceSurvivalStratCoxPHRegr (Wald)
        sc = make_strat_cox_cache(200L)
        list(desc = "fast_coxph_regression_prebuilt_cpp(cache, estimate_only=FALSE) [stratified]",
             REPS = 38000L,
             fn   = function() {
                 res = fast_coxph_regression_prebuilt_cpp(sc$cache, estimate_only = FALSE)
                 sqrt(res$vcov[1L, 1L])
             })
    }),

    weibull_var = local({
        # fast_weibull_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=FALSE)
        # Class: InferenceSurvivalWeibullRegr (Wald)
        d = make_data(200L, "cox")
        list(desc = "fast_weibull_regression_cpp(X_ord, y_bm, dead_bm, estimate_only=FALSE)",
             REPS = 94000L,
             fn   = function() {
                 res = fast_weibull_regression_cpp(d$X_ord, d$y_bm, d$dead_bm, estimate_only = FALSE)
                 sqrt(res$vcov[2L, 2L])
             })
    }),

    logrank_var = local({
        # EDI:::fast_logrank_stats_cpp(w_bm, y_bm, dead_bm) at N=200
        # Class: InferenceSurvivalLogRank (Wald)
        d = make_data(200L, "cox")
        list(desc = "EDI:::fast_logrank_stats_cpp(w_bm, y_bm, dead_bm) [N=200]",
             REPS = 1500000L,
             fn   = function() {
                 res = EDI:::fast_logrank_stats_cpp(d$w_bm, d$y_bm, as.integer(d$dead_bm))
                 stats::pchisq(res$score^2 / res$var_score, df = 1L, lower.tail = FALSE)
             })
    }),

    wilcox_var = local({
        # EDI:::wilcox_hl_point_estimate_cpp(w_bm, y_bm) at N=200
        # Class: InferenceAllSimpleWilcox (Wald)
        d = make_data(200L, "continuous")
        list(desc = "EDI:::wilcox_hl_point_estimate_cpp(w_bm, y_bm) [N=200]",
             REPS = 188000L,
             fn   = function() EDI:::wilcox_hl_point_estimate_cpp(d$w_bm, d$y_bm))
    }),

    gcomp_logistic_post_fit_var = local({
        # fast_logistic_regression_cpp + gcomp_logistic_post_fit_cpp (Wald)
        # Classes: InferencePropGCompMeanDiff, InferenceIncidGCompRiskDiff,
        #          InferenceIncidGCompRiskRatio (Wald paths)
        d = make_data(200L, "logistic")
        list(desc = "fast_logistic_regression_cpp + gcomp_logistic_post_fit_cpp [N=200]",
             REPS = 214000L,
             fn   = function() {
                 fit      = fast_logistic_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE)
                 coef_hat = as.numeric(fit$b)
                 mu_hat   = stats::plogis(as.numeric(d$X_bm %*% coef_hat))
                 res      = gcomp_logistic_post_fit_cpp(d$X_bm, d$y_bm, coef_hat, mu_hat, 2L)
                 res$rd / res$se_rd
             })
    }),

    gcomp_ordinal_var = local({
        # fast_ordinal_regression_with_var_cpp + gcomp_ordinal_proportional_odds_post_fit_cpp (Wald)
        # Class: InferenceOrdinalGCompMeanDiff (Wald)
        d = make_data(200L, "ordinal")
        list(desc = "fast_ordinal_regression_with_var_cpp + gcomp_ordinal_proportional_odds_post_fit_cpp [N=200]",
             REPS = 35000L,
             fn   = function() {
                 fit      = fast_ordinal_regression_with_var_cpp(d$X_ord, d$y_bm)
                 coef_hat = as.numeric(fit$b)
                 alpha_hat = as.numeric(fit$alpha)
                 res      = gcomp_ordinal_proportional_odds_post_fit_cpp(d$X_ord, coef_hat, alpha_hat, 1L)
                 res$md
             })
    }),

    # ================================================================
    # WEIGHTED ESTIMATE-ONLY PATHS  (IPW / bootstrap-weight code paths)
    # ================================================================

    logistic_weighted_est = local({
        d = make_data(1000L, "logistic")
        list(desc = "fast_logistic_regression_weighted_cpp(X_bm, y_bm, wts_bm)",
             REPS = 83000L,
             fn   = function() fast_logistic_regression_weighted_cpp(d$X_bm, d$y_bm, d$wts_bm))
    }),

    poisson_weighted_est = local({
        d = make_data(1000L, "poisson")
        list(desc = "fast_poisson_regression_weighted_cpp(X_bm, y_bm, wts_bm)",
             REPS = 50000L,
             fn   = function() fast_poisson_regression_weighted_cpp(d$X_bm, d$y_bm, d$wts_bm))
    }),

    probit_weighted_est = local({
        d = make_data(1000L, "logistic")
        list(desc = "fast_probit_regression_weighted_cpp(X_bm, y_bm, wts_bm)",
             REPS = 20000L,
             fn   = function() EDI:::fast_probit_regression_weighted_cpp(d$X_bm, d$y_bm, d$wts_bm))
    }),

    beta_weighted_est = local({
        d = make_data(1000L, "beta")
        list(desc = "fast_beta_regression_weighted_cpp(X_bm, y_bm, wts_bm, estimate_only=TRUE)",
             REPS = 7000L,
             fn   = function() fast_beta_regression_weighted_cpp(d$X_bm, d$y_bm, d$wts_bm, estimate_only = TRUE))
    }),

    logbin_weighted_est = local({
        d = make_data(1000L, "log-binomial")
        list(desc = "fast_log_binomial_regression_weighted_cpp(X_bm, y_bm, wts_bm, estimate_only=TRUE)",
             REPS = 8000L,
             fn   = function() EDI:::fast_log_binomial_regression_weighted_cpp(d$X_bm, d$y_bm, d$wts_bm, estimate_only = TRUE))
    }),

    identbin_weighted_est = local({
        d = make_data(1000L, "logistic")
        list(desc = "fast_identity_binomial_regression_weighted_cpp(X_bm, y_bm, wts_bm, estimate_only=TRUE)",
             REPS = 83000L,
             fn   = function() EDI:::fast_identity_binomial_regression_weighted_cpp(d$X_bm, d$y_bm, d$wts_bm, estimate_only = TRUE))
    }),

    prop_odds_weighted_est = local({
        d = make_data(1000L, "ordinal")
        list(desc = "fast_ordinal_regression_weighted_cpp(X_ord, y_bm, wts_bm)",
             REPS = 13000L,
             fn   = function() EDI:::fast_ordinal_regression_weighted_cpp(d$X_ord, d$y_bm, d$wts_bm))
    }),

    # ================================================================
    # MISSING FULL-INFERENCE PATHS (with-var not in original 59)
    # ================================================================

    identbin_var_full = local({
        d = make_data(200L, "logistic")
        list(desc = "fast_identity_binomial_regression_with_var_cpp(X_bm, y_bm, j=2L)",
             REPS = 26000L,
             fn   = function() EDI:::fast_identity_binomial_regression_with_var_cpp(d$X_bm, d$y_bm, j = 2L, maxit = 100L, tol = 1e-9))
    }),

    ord_cauchit_var = local({
        d = make_data(200L, "ordinal")
        list(desc = "fast_ordinal_cauchit_regression_with_var_cpp(X_ord, y_bm)",
             REPS = 10000L,
             fn   = function() EDI:::fast_ordinal_cauchit_regression_with_var_cpp(d$X_ord, d$y_bm))
    }),

    ord_cloglog_var = local({
        d = make_data(200L, "ordinal")
        list(desc = "fast_ordinal_cloglog_regression_with_var_cpp(X_ord, y_bm)",
             REPS = 8000L,
             fn   = function() EDI:::fast_ordinal_cloglog_regression_with_var_cpp(d$X_ord, d$y_bm))
    }),

    ord_probit_var = local({
        d = make_data(200L, "ordinal")
        list(desc = "fast_ordinal_probit_regression_with_var_cpp(X_ord, y_bm)",
             REPS = 7000L,
             fn   = function() EDI:::fast_ordinal_probit_regression_with_var_cpp(d$X_ord, d$y_bm))
    }),

    stereotype_est = local({
        d = make_data(1000L, "ordinal")
        list(desc = "fast_stereotype_logit_cpp(X_ord, y_bm, estimate_only=TRUE)",
             REPS = 700L,
             fn   = function() EDI:::fast_stereotype_logit_cpp(d$X_ord, d$y_bm, estimate_only = TRUE))
    }),

    stereotype_var = local({
        d = make_data(200L, "ordinal")
        list(desc = "fast_stereotype_logit_with_var_cpp(X_ord, y_bm, estimate_only=FALSE)",
             REPS = 2000L,
             fn   = function() EDI:::fast_stereotype_logit_with_var_cpp(d$X_ord, d$y_bm, estimate_only = FALSE))
    }),

    # ================================================================
    # NEW MODEL TYPES (not in original kernel list)
    # ================================================================

    trunc_negbin_est = local({
        d = make_data(1000L, "poisson")
        y_pos = pmax(1L, as.integer(d$y_bm))
        list(desc = "fast_truncated_negbin_count_cpp(X_bm, y_pos, estimate_only=TRUE)",
             REPS = 3000L,
             fn   = function() EDI:::fast_truncated_negbin_count_cpp(d$X_bm, y_pos, estimate_only = TRUE))
    }),

    zero_one_inflated_beta_est = local({
        d = make_data(1000L, "beta")
        set.seed(SEED + 1L)
        u = runif(length(d$y_bm))
        y_zoib = ifelse(u < 0.08, 0, ifelse(u < 0.16, 1, d$y_bm))
        X_zi = d$X_bm[, 1:2, drop = FALSE]
        list(desc = "fast_zero_one_inflated_beta_cpp(X_bm, X_zi, y_zoib, estimate_only=TRUE)",
             REPS = 200L,
             fn   = function() EDI:::fast_zero_one_inflated_beta_cpp(d$X_bm, X_zi, y_zoib, estimate_only = TRUE))
    }),

    weibull_frailty_est = local({
        d = make_glmm_data(400L, "cox")
        list(desc = "fast_weibull_frailty_cpp(X_ord, y_bm, dead_bm, group_id_bm, estimate_only=TRUE)",
             REPS = 400L,
             fn   = function() EDI:::fast_weibull_frailty_cpp(d$X_ord, d$y_bm, d$dead_bm, d$group_id_bm, estimate_only = TRUE))
    }),

    weibull_frailty_var = local({
        d = make_glmm_data(200L, "cox")
        list(desc = "fast_weibull_frailty_cpp(X_ord, y_bm, dead_bm, group_id_bm, estimate_only=FALSE)",
             REPS = 320L,
             fn   = function() EDI:::fast_weibull_frailty_cpp(d$X_ord, d$y_bm, d$dead_bm, d$group_id_bm, estimate_only = FALSE))
    }),

    # ================================================================
    # GLMM PATHS  (mixed-model kernels)
    # ================================================================

    logistic_glmm_est = local({
        d = make_glmm_data(400L, "logistic")
        list(desc = "fast_logistic_glmm_cpp(X_ord, y_bm, group_id_bm, j_T=0L, estimate_only=TRUE)",
             REPS = 1000L,
             fn   = function() EDI:::fast_logistic_glmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, j_T = 0L, estimate_only = TRUE))
    }),

    logistic_glmm_var = local({
        d = make_glmm_data(200L, "logistic")
        list(desc = "fast_logistic_glmm_cpp(X_ord, y_bm, group_id_bm, j_T=0L, estimate_only=FALSE)",
             REPS = 800L,
             fn   = function() EDI:::fast_logistic_glmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, j_T = 0L, estimate_only = FALSE))
    }),

    poisson_glmm_est = local({
        d = make_glmm_data(400L, "poisson")
        list(desc = "fast_poisson_glmm_cpp(X_ord, y_bm, group_id_bm, j_T=0L, estimate_only=TRUE)",
             REPS = 200L,
             fn   = function() EDI:::fast_poisson_glmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, j_T = 0L, estimate_only = TRUE))
    }),

    poisson_glmm_var = local({
        d = make_glmm_data(200L, "poisson")
        list(desc = "fast_poisson_glmm_cpp(X_ord, y_bm, group_id_bm, j_T=0L, estimate_only=FALSE)",
             REPS = 800L,
             fn   = function() EDI:::fast_poisson_glmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, j_T = 0L, estimate_only = FALSE))
    }),

    gaussian_lmm_est = local({
        d = make_glmm_data(400L, "continuous")
        list(desc = "fast_gaussian_lmm_cpp(X_ord, y_bm, group_id_bm, estimate_only=TRUE)",
             REPS = 20000L,
             fn   = function() EDI:::fast_gaussian_lmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, estimate_only = TRUE))
    }),

    gaussian_lmm_var = local({
        d = make_glmm_data(200L, "continuous")
        list(desc = "fast_gaussian_lmm_cpp(X_ord, y_bm, group_id_bm, estimate_only=FALSE)",
             REPS = 7000L,
             fn   = function() EDI:::fast_gaussian_lmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, estimate_only = FALSE))
    }),

    ordinal_glmm_est = local({
        d = make_glmm_data(400L, "ordinal")
        list(desc = "fast_ordinal_glmm_cpp(X_ord, as.integer(y_bm), group_id_bm, K=3L, j_T=0L, estimate_only=TRUE)",
             REPS = 70L,
             fn   = function() fast_ordinal_glmm_cpp(d$X_ord, as.integer(d$y_bm), d$group_id_bm, K = 3L, j_T = 0L, estimate_only = TRUE))
    }),

    ordinal_glmm_var = local({
        d = make_glmm_data(200L, "ordinal")
        list(desc = "fast_ordinal_glmm_cpp(X_ord, as.integer(y_bm), group_id_bm, K=3L, j_T=0L, estimate_only=FALSE)",
             REPS = 120L,
             fn   = function() fast_ordinal_glmm_cpp(d$X_ord, as.integer(d$y_bm), d$group_id_bm, K = 3L, j_T = 0L, estimate_only = FALSE))
    }),

    ordinal_clmm_est = local({
        d = make_glmm_data(400L, "ordinal")
        list(desc = "fast_ordinal_clmm_cpp(X_ord, as.integer(y_bm), group_id_bm, K=3L, j_T=0L, link='logit', estimate_only=TRUE)",
             REPS = 160L,
             fn   = function() EDI:::fast_ordinal_clmm_cpp(d$X_ord, as.integer(d$y_bm), d$group_id_bm, K = 3L, j_T = 0L, link = "logit", estimate_only = TRUE))
    }),

    hurdle_p_glmm_est = local({
        d = make_glmm_data(400L, "poisson")
        list(desc = "fast_hurdle_poisson_glmm_cpp(X_ord, y_bm, group_id_bm, j_T=0L, estimate_only=TRUE)",
             REPS = 600L,
             fn   = function() EDI:::fast_hurdle_poisson_glmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, j_T = 0L, estimate_only = TRUE))
    }),

    hurdle_p_glmm_var = local({
        d = make_glmm_data(200L, "poisson")
        list(desc = "fast_hurdle_poisson_glmm_cpp(X_ord, y_bm, group_id_bm, j_T=0L, estimate_only=FALSE)",
             REPS = 100L,
             fn   = function() EDI:::fast_hurdle_poisson_glmm_cpp(d$X_ord, d$y_bm, d$group_id_bm, j_T = 0L, estimate_only = FALSE))
    }),

    # ================================================================
    # GEE PATHS  (pairs/singletons GEE kernels)
    # ================================================================

    gee_pairs_singletons_logistic = local({
        d = make_glmm_data(200L, "logistic")
        list(desc = "gee_pairs_singletons_cpp(X_ord, y_bm, pairs_group_id_bm, 'logistic')",
             REPS = 20000L,
             fn   = function() EDI:::gee_pairs_singletons_cpp(d$X_ord, d$y_bm, d$pairs_group_id_bm, "logistic"))
    }),

    gee_pairs_singletons_weighted_logistic = local({
        d = make_glmm_data(200L, "logistic")
        list(desc = "gee_pairs_singletons_weighted_cpp(X_ord, y_bm, pairs_group_id_bm, 'logistic', wts_bm)",
             REPS = 20000L,
             fn   = function() EDI:::gee_pairs_singletons_weighted_cpp(d$X_ord, d$y_bm, d$pairs_group_id_bm, "logistic", d$wts_bm))
    }),

    # ================================================================
    # KK21 WEIGHT FUNCTIONS
    # ================================================================

    kk21_continuous_wts = local({
        d = make_data(500L, "continuous")
        X_cov = d$X[, -1L]
        list(desc = "kk21_continuous_weights_cpp(X_cov, y_bm)",
             REPS = 200000L,
             fn   = function() EDI:::kk21_continuous_weights_cpp(X_cov, d$y_bm))
    }),

    kk21_logistic_wts = local({
        d = make_data(500L, "logistic")
        X_cov = d$X[, -1L]
        list(desc = "kk21_logistic_weights_cpp(X_cov, y_bm, 25L, 1e-9)",
             REPS = 2000L,
             fn   = function() EDI:::kk21_logistic_weights_cpp(X_cov, d$y_bm, 25L, 1e-9))
    }),

    kk21_beta_wts = local({
        d = make_data(500L, "beta")
        X_cov = d$X[, -1L]
        list(desc = "kk21_beta_weights_cpp(X_cov, y_bm)",
             REPS = 500L,
             fn   = function() EDI:::kk21_beta_weights_cpp(X_cov, d$y_bm))
    }),

    kk21_negbin_wts = local({
        d = make_data(500L, "negbin")
        X_cov = d$X[, -1L]
        list(desc = "kk21_negbin_weights_cpp(X_cov, y_bm)",
             REPS = 2000L,
             fn   = function() EDI:::kk21_negbin_weights_cpp(X_cov, d$y_bm))
    }),

    kk21_ordinal_wts = local({
        d = make_data(500L, "ordinal")
        X_cov = d$X[, -1L]
        list(desc = "kk21_ordinal_weights_cpp(X_cov, y_bm)",
             REPS = 2000L,
             fn   = function() EDI:::kk21_ordinal_weights_cpp(X_cov, d$y_bm))
    }),

    kk21_survival_wts = local({
        d = make_data(500L, "cox")
        X_cov = d$X[, -1L]
        list(desc = "kk21_survival_weights_cpp(X_cov, y_bm, dead_bm)",
             REPS = 2000L,
             fn   = function() EDI:::kk21_survival_weights_cpp(X_cov, d$y_bm, as.numeric(d$dead_bm)))
    }),

    kk21_stepwise_continuous_wts = local({
        d = make_data(500L, "continuous")
        X_cov = d$X[, -1L]
        list(desc = "kk21_stepwise_continuous_weights_cpp(X_cov, y_bm, wts_bm)",
             REPS = 60000L,
             fn   = function() EDI:::kk21_stepwise_continuous_weights_cpp(X_cov, d$y_bm, d$wts_bm))
    }),

    kk21_stepwise_logistic_wts = local({
        d = make_data(500L, "logistic")
        X_cov = d$X[, -1L]
        list(desc = "kk21_stepwise_logistic_weights_cpp(X_cov, y_bm, wts_bm)",
             REPS = 6000L,
             fn   = function() EDI:::kk21_stepwise_logistic_weights_cpp(X_cov, d$y_bm, d$wts_bm))
    }),

    # ================================================================
    # STATS HELPERS (not exercised by existing inference kernels)
    # ================================================================

    newcombe_paired = local({
        list(desc = "newcombe_paired_ci_cpp(n11=40, n10=20, n01=10, n00=30, alpha=0.05)",
             REPS = 500000L,
             fn   = function() EDI:::newcombe_paired_ci_cpp(40, 20, 10, 30, 0.05))
    }),

    mn_ci = local({
        list(desc = "mn_ci_cpp(x_t=60, n_t=100, x_c=40, n_c=100, p_t=0.60, p_c=0.40, alpha=0.05, pval_eps=1e-10)",
             REPS = 500000L,
             fn   = function() EDI:::mn_ci_cpp(60, 100, 40, 100, 0.60, 0.40, 0.05, 1e-10))
    }),

    zhang_binom_pval = local({
        list(desc = "zhang_exact_binom_pval_cpp(d_plus=15, d_minus=5, delta_0=0.0)",
             REPS = 500000L,
             fn   = function() EDI:::zhang_exact_binom_pval_cpp(15L, 5L, 0.0))
    }),

    zhang_fisher_pval = local({
        list(desc = "zhang_exact_fisher_pval_cpp(n11=30, n10=10, n01=5, n00=55, delta_0=0.0)",
             REPS = 500000L,
             fn   = function() EDI:::zhang_exact_fisher_pval_cpp(30L, 10L, 5L, 55L, 0.0))
    }),

    # ================================================================
    # POST-FIT HELPERS (variance-computation helpers)
    # ================================================================

    gcomp_frac_logit_post_fit_var = local({
        d = make_data(200L, "beta")
        fit = fast_logistic_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE)
        coef_hat = as.numeric(fit$b)
        mu_hat   = stats::plogis(as.numeric(d$X_bm %*% coef_hat))
        list(desc = "gcomp_fractional_logit_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, j_treat=2L)",
             REPS = 100000L,
             fn   = function() EDI:::gcomp_fractional_logit_post_fit_cpp(d$X_bm, d$y_bm, coef_hat, mu_hat, 2L))
    }),

    glm_sandwich_post_fit_var = local({
        d = make_data(200L, "logistic")
        fit = fast_logistic_regression_cpp(d$X_bm, d$y_bm, estimate_only = TRUE)
        coef_hat = as.numeric(fit$b)
        mu_hat   = stats::plogis(as.numeric(d$X_bm %*% coef_hat))
        ww       = mu_hat * (1 - mu_hat)
        list(desc = "glm_sandwich_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, working_weights, j_treat=2L)",
             REPS = 100000L,
             fn   = function() EDI:::glm_sandwich_post_fit_cpp(d$X_bm, d$y_bm, coef_hat, mu_hat, ww, 2L))
    }),

    ordinal_gcomp_post_fit_var = local({
        d = make_data(200L, "ordinal")
        fit = fast_ordinal_regression_with_var_cpp(d$X_ord, d$y_bm)
        coef_hat  = as.numeric(fit$b)
        alpha_hat = as.numeric(fit$alpha)
        list(desc = "ordinal_gcomp_post_fit_cpp(X_ord, y_bm, coef_hat, alpha_hat, j_treat=1L)",
             REPS = 30000L,
             fn   = function() EDI:::ordinal_gcomp_post_fit_cpp(d$X_ord, d$y_bm, coef_hat, alpha_hat, 1L))
    }),

    # ================================================================
    # BOOTSTRAP INDEX GENERATORS
    # ================================================================

    bootstrap_indices = local({
        list(desc = "bootstrap_indices_cpp(n=1000L, B=500L)",
             REPS = 500L,
             fn   = function() EDI:::bootstrap_indices_cpp(1000L, 500L))
    }),

    stratified_bootstrap_indices = local({
        set.seed(SEED)
        strata_keys = as.character(rep(paste0("S", 1:5), each = 200L))
        list(desc = "stratified_bootstrap_indices_cpp(strata_keys) [5 strata x 200]",
             REPS = 40000L,
             fn   = function() EDI:::stratified_bootstrap_indices_cpp(strata_keys))
    }),

    bootstrap_m_indices = local({
        set.seed(SEED)
        n_pairs = 250L; n = 500L
        m_vec      = as.integer(rep(seq_len(n_pairs), each = 2L))
        i_reservoir = as.integer(seq_len(n))
        list(desc = "bootstrap_m_indices_cpp(m_vec, i_reservoir, n_reservoir=500L, m=250L, B=500L)",
             REPS = 500L,
             fn   = function() EDI:::bootstrap_m_indices_cpp(m_vec, i_reservoir, n, n_pairs, 500L))
    }),

    # ================================================================
    # RANDOMIZATION PRIMITIVES
    # ================================================================

    complete_randomization_balanced = local({
        list(desc = "complete_randomization_forced_balanced_cpp(n=1000L, r=1000L, seed=42L)",
             REPS = 50L,
             fn   = function() EDI:::complete_randomization_forced_balanced_cpp(1000L, 1000L, 42L))
    }),

    complete_randomization_imbalanced = local({
        list(desc = "complete_randomization_imbalanced_cpp(n=1000L, nT=300L, r=1000L, seed=42L)",
             REPS = 50L,
             fn   = function() EDI:::complete_randomization_imbalanced_cpp(1000L, 300L, 1000L, 42L))
    }),

    shuffle_w = local({
        set.seed(SEED)
        w = as.numeric(c(rep(1, 500L), rep(0, 500L)))
        list(desc = "shuffle_cpp(w) [n=1000]",
             REPS = 80000L,
             fn   = function() EDI:::shuffle_cpp(w))
    }),

    efron_redraw = local({
        list(desc = "efron_redraw_cpp(t=500L, prob_T=0.5, weighted_coin_prob=2/3)",
             REPS = 200000L,
             fn   = function() EDI:::efron_redraw_cpp(500L, 0.5, 2/3))
    }),

    spbr_redraw = local({
        set.seed(SEED)
        strata_keys = as.character(rep(c("A","B","C","D"), each = 4L))
        list(desc = "spbr_redraw_w_cpp(strata_keys_block, block_size=4L, prob_T=0.5)",
             REPS = 400000L,
             fn   = function() EDI:::spbr_redraw_w_cpp(strata_keys, 4L, 0.5))
    }),

    random_block_size_redraw = local({
        set.seed(SEED)
        strata_keys  = as.character(rep(c("A","B","C","D"), each = 4L))
        block_sizes  = as.integer(rep(4L, 4L))
        list(desc = "random_block_size_redraw_w_cpp(strata_keys_block, block_sizes, prob_T=0.5)",
             REPS = 500000L,
             fn   = function() EDI:::random_block_size_redraw_w_cpp(strata_keys, block_sizes, 0.5))
    }),

    redraw_w_kk14 = local({
        set.seed(SEED)
        n_pairs = 250L; n = 500L
        m_vec = as.integer(rep(seq_len(n_pairs), each = 2L))
        w_cur = as.numeric(c(rep(1, n_pairs), rep(0, n - n_pairs)))
        list(desc = "redraw_w_kk14_cpp(m_vec, w) [250 pairs]",
             REPS = 200000L,
             fn   = function() EDI:::redraw_w_kk14_cpp(m_vec, w_cur))
    }),

    atkinson_redraw = local({
        set.seed(SEED)
        n = 100L; p_raw = 4L
        X_raw = matrix(rnorm(n * p_raw), n, p_raw)
        list(desc = "atkinson_redraw_batch_cpp(X_raw, n=100L, p_raw=4L, prob_T=0.5)",
             REPS = 4000L,
             fn   = function() EDI:::atkinson_redraw_batch_cpp(X_raw, n, p_raw, 0.5))
    }),

    pocock_simon_assign = local({
        set.seed(SEED)
        counts = matrix(as.double(c(10L, 8L, 12L, 9L, 11L, 10L)), nrow = 3L, ncol = 2L)
        levels_idx = as.integer(c(0L, 1L, 2L))
        weights    = c(1.0, 1.0, 1.0)
        list(desc = "pocock_simon_assign_cpp(counts, levels_idx, weights, p_best=0.7, prob_T=0.5)",
             REPS = 400000L,
             fn   = function() EDI:::pocock_simon_assign_cpp(counts, levels_idx, weights, 0.7, 0.5))
    }),

    pocock_simon_assign_and_update = local({
        set.seed(SEED)
        counts = matrix(as.double(c(10L, 8L, 12L, 9L, 11L, 10L)), nrow = 3L, ncol = 2L)
        levels_idx = as.integer(c(0L, 1L, 2L))
        weights    = c(1.0, 1.0, 1.0)
        list(desc = "pocock_simon_assign_and_update_cpp(counts, levels_idx, weights, p_best=0.7, prob_T=0.5)",
             REPS = 400000L,
             fn   = function() EDI:::pocock_simon_assign_and_update_cpp(counts, levels_idx, weights, 0.7, 0.5))
    }),

    pocock_simon_redraw_w = local({
        set.seed(SEED)
        n = 200L; n_factors = 3L; n_levels_each = 2L
        x_levels_matrix = matrix(as.integer(sample(0:(n_levels_each - 1L), n * n_factors, replace = TRUE)),
                                  nrow = n, ncol = n_factors)
        num_levels_total = n_factors * n_levels_each
        weights = rep(1.0, n_factors)
        list(desc = "pocock_simon_redraw_w_cpp(x_levels_matrix, num_levels_total=6L, weights, p_best=0.7, prob_T=0.5) [n=200]",
             REPS = 20000L,
             fn   = function() EDI:::pocock_simon_redraw_w_cpp(x_levels_matrix, num_levels_total, weights, 0.7, 0.5))
    }),

    # ================================================================
    # PERMUTATION GENERATORS
    # ================================================================

    generate_permutations_bernoulli = local({
        list(desc = "generate_permutations_bernoulli_cpp(n=1000L, nsim=1000L, prob_T=0.5)",
             REPS = 200L,
             fn   = function() EDI:::generate_permutations_bernoulli_cpp(1000L, 1000L, 0.5))
    }),

    generate_permutations_efron = local({
        list(desc = "generate_permutations_efron_cpp(n=1000L, nsim=1000L, prob_T=0.5, wcoin_prob=2/3)",
             REPS = 200L,
             fn   = function() EDI:::generate_permutations_efron_cpp(1000L, 1000L, 0.5, 2/3))
    }),

    generate_permutations_ibcrd = local({
        list(desc = "generate_permutations_ibcrd_cpp(n=1000L, nsim=1000L, prob_T=0.5)",
             REPS = 200L,
             fn   = function() EDI:::generate_permutations_ibcrd_cpp(1000L, 1000L, 0.5))
    }),

    generate_permutations_blocking = local({
        n_blocks = 10L; block_size = 100L; n = n_blocks * block_size
        strata_indices = lapply(seq_len(n_blocks), function(i) as.integer(((i - 1L) * block_size + 1L):(i * block_size)))
        list(desc = "generate_permutations_blocking_cpp(n=1000L, nsim=500L, prob_T=0.5, 10 blocks x 100)",
             REPS = 300L,
             fn   = function() EDI:::generate_permutations_blocking_cpp(n, 500L, 0.5, strata_indices))
    }),

    generate_permutations_matching = local({
        n_pairs = 500L
        m_vec = as.integer(rep(seq_len(n_pairs), each = 2L))
        list(desc = "generate_permutations_matching_cpp(m_vec, nsim=1000L, prob_T=0.5) [500 pairs]",
             REPS = 400L,
             fn   = function() EDI:::generate_permutations_matching_cpp(m_vec, 1000L, 0.5))
    }),

    generate_permutations_atkinson = local({
        set.seed(SEED)
        n = 100L; p_raw = 4L
        X_raw = matrix(rnorm(n * p_raw), n, p_raw)
        list(desc = "generate_permutations_atkinson_cpp(X_raw, n=100L, p_raw=4L, prob_T=0.5, nsim=100L)",
             REPS = 50L,
             fn   = function() EDI:::generate_permutations_atkinson_cpp(X_raw, n, p_raw, 0.5, 100L))
    }),

    generate_permutations_pocock_simon = local({
        set.seed(SEED)
        n = 200L; n_factors = 3L; n_levels_each = 2L
        x_levels_matrix = matrix(as.integer(sample(0:(n_levels_each - 1L), n * n_factors, replace = TRUE)),
                                  nrow = n, ncol = n_factors)
        num_levels_total = n_factors * n_levels_each
        weights = rep(1.0, n_factors)
        list(desc = "generate_permutations_pocock_simon_cpp(x_levels_matrix, 6L, weights, p_best=0.7, prob_T=0.5, nsim=500L)",
             REPS = 100L,
             fn   = function() EDI:::generate_permutations_pocock_simon_cpp(x_levels_matrix, num_levels_total, weights, 0.7, 0.5, 500L))
    }),

    generate_permutations_cluster = local({
        n_clusters = 20L; cluster_size = 50L; n = n_clusters * cluster_size
        cluster_indices = lapply(seq_len(n_clusters), function(i) as.integer(((i - 1L) * cluster_size + 1L):(i * cluster_size)))
        list(desc = "generate_permutations_cluster_cpp(n=1000L, nsim=500L, prob_T=0.5, 20 clusters x 50)",
             REPS = 700L,
             fn   = function() EDI:::generate_permutations_cluster_cpp(n, 500L, 0.5, cluster_indices))
    }),

    generate_permutations_spbr = local({
        set.seed(SEED)
        strata_keys = as.character(rep(paste0("S", 1:5), each = 200L))
        list(desc = "generate_permutations_spbr_cpp(strata_keys, block_size=4L, prob_T=0.5, nsim=500L) [5 strata x 200]",
             REPS = 100L,
             fn   = function() EDI:::generate_permutations_spbr_cpp(strata_keys, 4L, 0.5, 500L))
    }),

    # ================================================================
    # BOOTSTRAP / RANDOMIZATION LOOPS
    # ================================================================

    draw_binary_match_assignments = local({
        set.seed(SEED)
        n_pairs = 500L; n = 1000L
        indices_pairs = matrix(as.integer(matrix(seq_len(n), nrow = 2L)), ncol = 2L, byrow = FALSE)
        list(desc = "draw_binary_match_assignments_cpp(indices_pairs, n=1000L, r=200L, num_cores=1L)",
             REPS = 700L,
             fn   = function() EDI:::draw_binary_match_assignments_cpp(indices_pairs, n, 200L, 1L))
    }),

    draw_matching_bootstrap_sample = local({
        set.seed(SEED)
        n_pairs = 250L; n = 500L
        i_reservoir = as.integer(seq_len(n))
        pair_rows = matrix(as.integer(matrix(seq_len(n), nrow = 2L)), ncol = 2L, byrow = FALSE)
        list(desc = "draw_matching_bootstrap_sample_cpp(i_reservoir, pair_rows, n_reservoir=500L)",
             REPS = 30000L,
             fn   = function() EDI:::draw_matching_bootstrap_sample_cpp(i_reservoir, pair_rows, n))
    }),

    randomization_loop = local({
        dup_design_fn    = function() list()
        dup_inference_fn = function() list()
        run_rand_iter_fn = function(objects) runif(1L)
        list(desc = "randomization_loop_cpp(r=200L, stub_fns, num_cores=1L)",
             REPS = 100000L,
             fn   = function() EDI:::randomization_loop_cpp(200L, dup_design_fn, dup_inference_fn, run_rand_iter_fn, 1L))
    }),

    base_bootstrap_loop = local({
        d = make_data(100L, "continuous")
        indices = EDI:::bootstrap_indices_cpp(100L, 50L)
        dup_inference_fn    = function() NULL
        compute_estimate_fn = function(data) as.numeric(fast_ols_cpp(data$X, data$y)$b[2L])
        list(desc = "base_bootstrap_loop_cpp(X_bm, y_bm, dead=NA, w_bm, indices, ols_fns, num_cores=1L) [n=100, B=50]",
             REPS = 2000L,
             fn   = function() EDI:::base_bootstrap_loop_cpp(d$X_bm, d$y_bm, as.numeric(rep(NA, 100L)),
                                                              as.numeric(d$w_bm), indices,
                                                              dup_inference_fn, compute_estimate_fn, 1L))
    }),

    # ================================================================
    # MISSING COVERAGE PATHS — added 2026-07-04
    # ================================================================

    # --- ClogitPlusGLMM (fast_clogit_plus_glmm.cpp) ---
    clogit_glmm_est = local({
        set.seed(SEED)
        n_disc = 200L; n_conc = 200L; p = 4L; n_grp = 80L
        X_disc = matrix(rnorm(n_disc * (p + 1L)), n_disc, p + 1L)
        X_disc[, 1L] = 1; X_disc[, 2L] = sample(c(-1L, 1L), n_disc, replace = TRUE)
        y_disc = as.numeric(rbinom(n_disc, 1L, 0.4))
        X_conc = matrix(rnorm(n_conc * (p + 1L)), n_conc, p + 1L)
        X_conc[, 1L] = 1; X_conc[, 2L] = as.numeric(rep(c(0L, 1L), n_conc / 2L))
        y_conc = as.numeric(rbinom(n_conc, 1L, 0.4))
        group_conc = as.integer(rep(seq_len(n_grp), length.out = n_conc))
        list(desc = "fast_clogit_plus_glmm_cpp(..., estimate_only=TRUE) [n_disc=200, n_conc=200, p=4, G=80]",
             REPS = 1000L,
             fn   = function() EDI:::fast_clogit_plus_glmm_cpp(
                 X_disc, y_disc, X_conc, y_conc, group_conc,
                 has_discordant = TRUE, has_concordant = TRUE, estimate_only = TRUE))
    }),

    clogit_glmm_var = local({
        set.seed(SEED)
        n_disc = 200L; n_conc = 200L; p = 4L; n_grp = 80L
        X_disc = matrix(rnorm(n_disc * (p + 1L)), n_disc, p + 1L)
        X_disc[, 1L] = 1; X_disc[, 2L] = sample(c(-1L, 1L), n_disc, replace = TRUE)
        y_disc = as.numeric(rbinom(n_disc, 1L, 0.4))
        X_conc = matrix(rnorm(n_conc * (p + 1L)), n_conc, p + 1L)
        X_conc[, 1L] = 1; X_conc[, 2L] = as.numeric(rep(c(0L, 1L), n_conc / 2L))
        y_conc = as.numeric(rbinom(n_conc, 1L, 0.4))
        group_conc = as.integer(rep(seq_len(n_grp), length.out = n_conc))
        list(desc = "fast_clogit_plus_glmm_cpp(..., estimate_only=FALSE) [n_disc=200, n_conc=200, p=4, G=80]",
             REPS = 300L,
             fn   = function() EDI:::fast_clogit_plus_glmm_cpp(
                 X_disc, y_disc, X_conc, y_conc, group_conc,
                 has_discordant = TRUE, has_concordant = TRUE, estimate_only = FALSE))
    }),

    # --- Dep-cens-transform survival (fast_survival_models_optim.cpp) ---
    dep_cens_transform_est = local({
        d = make_data(500L, "cox")
        list(desc = "fast_dep_cens_transform_optim_cpp(X_ord, y_bm, dead_bm, smart_cold_start=TRUE, estimate_only=TRUE) [n=500]",
             REPS = 200L,
             fn   = function() EDI:::fast_dep_cens_transform_optim_cpp(
                 d$X_ord, d$y_bm, as.numeric(d$dead_bm),
                 smart_cold_start = TRUE, estimate_only = TRUE))
    }),

    dep_cens_transform_var = local({
        d = make_data(500L, "cox")
        list(desc = "fast_dep_cens_transform_optim_cpp(X_ord, y_bm, dead_bm, smart_cold_start=TRUE, estimate_only=FALSE) [n=500]",
             REPS = 100L,
             fn   = function() EDI:::fast_dep_cens_transform_optim_cpp(
                 d$X_ord, d$y_bm, as.numeric(d$dead_bm),
                 smart_cold_start = TRUE, estimate_only = FALSE))
    }),

    # --- D-optimal design search (optimal_design_search.cpp) ---
    d_optimal_search = local({
        set.seed(SEED)
        n = 60L; p = 5L; n_T = 30L; nsim = 500L
        X = cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
        XtX_inv = solve(crossprod(X) + diag(1e-8, p))
        P = X %*% XtX_inv %*% t(X)
        list(desc = "d_optimal_search_cpp(P, nsim=500, n_T=30) [n=60, p=5]",
             REPS = 200L,
             fn   = function() EDI:::d_optimal_search_cpp(P, nsim, n_T))
    }),

    # --- KK compound distribution (kk_compound_distr_parallel.cpp) ---
    kk_compound_distr = local({
        set.seed(SEED)
        n = 400L; nsim = 2000L; n_pairs = n / 2L
        y = rnorm(n)
        w_mat = matrix(as.integer(sample(c(0L, 1L), n * nsim, replace = TRUE)), nrow = n, ncol = nsim)
        m_mat = matrix(as.integer(rep(seq_len(n_pairs), each = 2L)), nrow = n, ncol = nsim)
        list(desc = "compute_matching_compound_distr_parallel_cpp(y, w_mat, m_mat, num_cores=1) [n=400, nsim=2000]",
             REPS = 200L,
             fn   = function() EDI:::compute_matching_compound_distr_parallel_cpp(y, w_mat, m_mat, 1L))
    }),

    # --- BAI parallel statistic (fast_bai_parallel.cpp) ---
    bai_distr = local({
        set.seed(SEED)
        n = 400L; nsim = 2000L; n_pairs = n / 2L; n_halves = n_pairs / 2L
        y = rnorm(n)
        w_mat = matrix(as.integer(sample(c(0L, 1L), n * nsim, replace = TRUE)), nrow = n, ncol = nsim)
        m_mat = matrix(as.integer(rep(seq_len(n_pairs), each = 2L)), nrow = n, ncol = nsim)
        pair_ids = seq_len(n_pairs)
        halves_idx = matrix(as.integer(sample(pair_ids, n_halves * 2L, replace = FALSE)), nrow = n_halves, ncol = 2L)
        list(desc = "compute_bai_distr_parallel_cpp(w_mat, m_mat, y, delta=0, halves_idx, convex_flag=TRUE, num_cores=1) [n=400, nsim=2000]",
             REPS = 200L,
             fn   = function() compute_bai_distr_parallel_cpp(w_mat, m_mat, y, 0, halves_idx, TRUE, 1L))
    }),

    # --- Rerandomization search (rerandomization_helpers.cpp) ---
    rerandomization_search = local({
        set.seed(SEED)
        n = 200L; p = 4L; r = 500L
        X = matrix(rnorm(n * p), n, p)
        # cutoff of ~1.0 accepts ~10% of randomizations for abs_sum_diff
        list(desc = "rerandomization_search_cpp(X, r=500, 'abs_sum_diff', cutoff=1.0, max_draws=50000) [n=200, p=4]",
             REPS = 50L,
             fn   = function() EDI:::rerandomization_search_cpp(X, r, "abs_sum_diff", 1.0, 50000L))
    }),

    rerandomization_obj_vals = local({
        set.seed(SEED)
        n = 200L; p = 4L; r = 5000L
        X = matrix(rnorm(n * p), n, p)
        indicTs = matrix(as.integer(sample(c(0L, 1L), n * r, replace = TRUE)), nrow = r, ncol = n)
        list(desc = "compute_objective_vals_cpp(X, indicTs, 'abs_sum_diff') [n=200, p=4, r=5000]",
             REPS = 100L,
             fn   = function() EDI:::compute_objective_vals_cpp(X, indicTs, "abs_sum_diff"))
    }),

    # --- CMH block SE (cmh_speedups.cpp) ---
    cmh_block_se = local({
        set.seed(SEED)
        n = 2000L; n_blocks = 500L
        y    = as.numeric(rbinom(n, 1L, 0.4))
        m_vec = as.integer(rep(seq_len(n_blocks), each = n / n_blocks))
        list(desc = "compute_cmh_block_se_cpp(y, m_vec, n_total=2000) [n=2000, B=500]",
             REPS = 20000L,
             fn   = function() EDI:::compute_cmh_block_se_cpp(y, m_vec, n))
    }),

    # Unknown kernel
    stop("Unknown kernel: '", name, "'. Run without args to see list.")
    )
}

k = get_kernel(KERNEL)

# Validation
cat("Kernel:     ", KERNEL, "\n")
cat("Seed:       ", SEED, "\n")
cat("Desc:       ", k$desc, "\n")
cat("REPS:       ", k$REPS, "\n")
r = k$fn()
cat("Result[1]:  ", if (is.numeric(r)) round(r[1L], 6) else "non-numeric", "\n")

# Timing warmup (1000 reps or REPS if smaller)
WARM = min(1000L, k$REPS)
t0 = proc.time()["elapsed"]
for (i in seq_len(WARM)) k$fn()
ms_per_call = (proc.time()["elapsed"] - t0) / WARM * 1000
cat("ms/call:    ", round(ms_per_call, 4), "\n")
cat("Est. loop:  ", round(ms_per_call * k$REPS / 1000, 1), "s\n")

# Perf target loop
cat("Starting perf loop (", k$REPS, " reps)...\n", sep = "")
for (i in seq_len(k$REPS)) k$fn()
cat("Done.\n")
