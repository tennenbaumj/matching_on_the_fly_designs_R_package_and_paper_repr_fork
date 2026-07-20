#!/usr/bin/env Rscript
# BRT scenario timing probe: diverse designs + inference classes
suppressPackageStartupMessages({
  library(EDI)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})
set_num_cores(1L); toggle_asserts(FALSE)

B  = 151L
timed = function(lbl, expr) {
  t0 = proc.time()
  r  = tryCatch(expr, error = function(e) { cat("  ERR:", conditionMessage(e), "\n"); NULL })
  cat(sprintf("  [%-58s] %.2fs\n", lbl, (proc.time() - t0)[3]))
  invisible(r)
}

# ── Common data ────────────────────────────────────────────
set.seed(42)
N = 148L
ds = ggplot2::diamonds %>% na.omit() %>%
  slice_sample(n = N, replace = TRUE) %>%
  mutate_if(where(is.factor), as.character)
X1  = data.frame(x1 = scale(as.numeric(ds$carat)))       # 1-covariate (fast designs)
X_f = ds %>% model.matrix(price ~ 0 + ., .) %>%
        apply(2, scale) %>% `/`(ncol(.)) %>%
        data.table() %>% select(where(~ !any(is.na(.)))) %>% as.data.frame()
y_cont = log(ds$price)
y_bin  = as.integer(ds$cut %in% c("Good","Very Good","Premium","Ideal"))
y_cnt  = pmax(0L, as.integer(round(ds$price / 500)))
y_ord  = as.character(as.integer(cut(y_cont, 4, labels = 1:4)))
set.seed(99)
surv_t = pmax(0.01, rexp(N, 0.1) * exp(-0.2 * X1$x1))
surv_d = rbinom(N, 1, 0.7)

# ── Design builders ────────────────────────────────────────
bern_des = function(resp, X, y, extra = list()) {
  des = DesignSeqOneByOneBernoulli$new(n = N, response_type = resp)
  for (i in seq_len(N)) {
    des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
    do.call(des$add_one_subject_response, c(list(i, y[i]), extra))
  }
  des
}
greedy_des = function(resp, X, y, extra = list()) {
  des = DesignFixedGreedy$new(response_type = resp, n = N, design_formula = ~ .)
  des$add_all_subjects_to_experiment(X)
  des$assign_w_to_all_subjects()
  for (i in seq_len(N))
    do.call(des$add_one_subject_response, c(list(i, y[i]), extra))
  des
}

# ── BRT suite (pval only) ──────────────────────────────────
brt_pval = function(inf, lbl) {
  cat(sprintf("\n--- %s ---\n", lbl))
  timed(paste0(lbl, " | percentile"),             inf$compute_rand_bootstrap_two_sided_pval(B=B, show_progress=FALSE))
  timed(paste0(lbl, " | studentized"),            inf$compute_rand_bootstrap_two_sided_pval(B=B, type="studentized",           show_progress=FALSE))
  timed(paste0(lbl, " | symmetric-percentile-t"), inf$compute_rand_bootstrap_two_sided_pval(B=B, type="symmetric-percentile-t", show_progress=FALSE))
  timed(paste0(lbl, " | smoothed"),               inf$compute_rand_bootstrap_two_sided_pval(B=B, type="smoothed",              show_progress=FALSE))
}
brt_ci = function(inf, lbl) {
  timed(paste0(lbl, " | CI percentile"),             inf$compute_rand_bootstrap_confidence_interval(B=B, pval_epsilon=0.05, show_progress=FALSE))
  timed(paste0(lbl, " | CI studentized"),            inf$compute_rand_bootstrap_confidence_interval(B=B, type="studentized",           pval_epsilon=0.05, show_progress=FALSE))
  timed(paste0(lbl, " | CI symmetric-percentile-t"), inf$compute_rand_bootstrap_confidence_interval(B=B, type="symmetric-percentile-t", pval_epsilon=0.05, show_progress=FALSE))
  timed(paste0(lbl, " | CI smoothed"),               inf$compute_rand_bootstrap_confidence_interval(B=B, type="smoothed",              pval_epsilon=0.05, show_progress=FALSE))
}

cat("\n### PVAL timings (B=151) ###\n")

# 1. MeanDiffPooledVar + Bernoulli (fast design, C++ kernel)
set.seed(1)
brt_pval(InferenceAllSimpleMeanDiffPooledVar$new(bern_des("continuous", X1, y_cont)),
         "MeanDiffPooledVar/Bernoulli")

# 2. MeanDiffPooledVar + FixedGreedy (slow design, C++ kernel)
set.seed(2)
brt_pval(InferenceAllSimpleMeanDiffPooledVar$new(greedy_des("continuous", X_f, y_cont)),
         "MeanDiffPooledVar/FixedGreedy")

# 3. OLS (worker=TRUE) + Bernoulli
set.seed(3)
brt_pval(InferenceContinOLS$new(bern_des("continuous", X1, y_cont), model_formula = ~ 1),
         "ContinOLS(~1)/Bernoulli")

# 4. OLS (worker=TRUE) + FixedGreedy
set.seed(4)
brt_pval(InferenceContinOLS$new(greedy_des("continuous", X_f, y_cont), model_formula = ~ 1),
         "ContinOLS(~1)/FixedGreedy")

# 5. Poisson (worker=TRUE) + Bernoulli
set.seed(5)
brt_pval(InferenceCountPoisson$new(bern_des("count", X1, y_cnt), model_formula = ~ 1),
         "CountPoisson(~1)/Bernoulli")

# 6. HurdlePoisson (worker=FALSE) + Bernoulli
set.seed(6)
brt_pval(InferenceCountHurdlePoisson$new(bern_des("count", X1, y_cnt), model_formula = ~ 1),
         "CountHurdlePoisson(~1)/Bernoulli")

# 7. ZIP (worker=FALSE) + Bernoulli
set.seed(7)
brt_pval(InferenceCountZeroInflatedPoisson$new(bern_des("count", X1, y_cnt), model_formula = ~ 1),
         "CountZIP(~1)/Bernoulli")

# 8. RiskDiff (worker=TRUE) + Bernoulli
set.seed(8)
brt_pval(InferenceIncidRiskDiff$new(bern_des("incidence", X1, y_bin), model_formula = ~ 1),
         "IncidRiskDiff(~1)/Bernoulli")

# 9. OrdinalAdjCatLogit (worker=TRUE) + Bernoulli
set.seed(9)
brt_pval(InferenceOrdinalAdjCatLogitRegr$new(bern_des("ordinal", X1, y_ord), model_formula = ~ 1),
         "OrdinalAdjCatLogit(~1)/Bernoulli")

# 10. CoxPH (check worker path) + Bernoulli
set.seed(10)
des_surv = DesignSeqOneByOneBernoulli$new(n = N, response_type = "survival")
for (i in seq_len(N)) {
  des_surv$add_one_subject_to_experiment_and_assign(X1[i, , drop=FALSE])
  des_surv$add_one_subject_response(i, surv_t[i], dead = surv_d[i])
}
brt_pval(InferenceSurvivalCoxPHRegr$new(des_surv, model_formula = ~ 1),
         "SurvivalCoxPH(~1)/Bernoulli")

cat("\n### CI timings (B=151, pval_eps=0.05) ###\n")

# MeanDiff + Bernoulli (fast)
set.seed(20)
brt_ci(InferenceAllSimpleMeanDiffPooledVar$new(bern_des("continuous", X1, y_cont)),
       "MeanDiffPooledVar/Bernoulli")

# MeanDiff + FixedGreedy (slow design draw)
set.seed(21)
brt_ci(InferenceAllSimpleMeanDiffPooledVar$new(greedy_des("continuous", X_f, y_cont)),
       "MeanDiffPooledVar/FixedGreedy")

# OLS + Bernoulli
set.seed(22)
brt_ci(InferenceContinOLS$new(bern_des("continuous", X1, y_cont), model_formula = ~ 1),
       "ContinOLS(~1)/Bernoulli")

# Poisson + Bernoulli
set.seed(23)
brt_ci(InferenceCountPoisson$new(bern_des("count", X1, y_cnt), model_formula = ~ 1),
       "CountPoisson(~1)/Bernoulli")

# HurdlePoisson + Bernoulli (worker=FALSE — slow path)
set.seed(24)
brt_ci(InferenceCountHurdlePoisson$new(bern_des("count", X1, y_cnt), model_formula = ~ 1),
       "CountHurdlePoisson(~1)/Bernoulli")

cat("\nDone.\n")
