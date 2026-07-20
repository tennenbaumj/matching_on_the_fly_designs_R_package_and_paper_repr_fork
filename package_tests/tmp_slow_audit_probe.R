#!/usr/bin/env Rscript
# Timing probe: audit all currently-SLOW operations in path_audits_source.R
suppressPackageStartupMessages({ library(EDI) })
set_num_cores(1L); toggle_asserts(FALSE)
B = 151L; r = 151L; pval_eps = 0.05

timed = function(lbl, expr) {
  t0 = proc.time()
  res = tryCatch(suppressWarnings(expr), error = function(e) {
    cat("  ERR:", conditionMessage(e), "\n"); NULL
  })
  elapsed = (proc.time() - t0)[3]
  status = if (elapsed < 30) "OK" else "SLOW"
  cat(sprintf("  [%s] %-65s %.1fs\n", status, lbl, elapsed))
  invisible(res)
}

N = 148L
set.seed(42)
X1 = data.frame(x1 = rnorm(N))

# ── Responses ────────────────────────────────────────────────
y_cnt  = pmax(0L, as.integer(rpois(N, 2)))
y_inc  = as.integer(rbinom(N, 1, 0.5))
y_ord  = as.character(sample(1:4, N, replace = TRUE))
y_prop = pmin(0.99, pmax(0.01, rbeta(N, 2, 3)))

# ── Design builders ──────────────────────────────────────────
bern_des = function(resp, y) {
  des = DesignSeqOneByOneBernoulli$new(n = N, response_type = resp)
  for (i in seq_len(N)) {
    des$add_one_subject_to_experiment_and_assign(X1[i, , drop = FALSE])
    des$add_one_subject_response(i, y[i])
  }
  des
}

# ─────────────────────────────────────────────────────────────
cat("\n=== 1. boot_bayes/ci (currently skip_bbt_ci_slow) ===\n")
# CountHurdlePoisson — our worker=TRUE fix should help
{
  des = bern_des("count", y_cnt)
  inf = InferenceCountHurdlePoisson$new(des, model_formula = ~ 1)
  timed("CountHurdlePoisson boot_bayes/ci",
        inf$compute_bayesian_bootstrap_confidence_interval(B = B, na.rm = TRUE, show_progress = FALSE))
}
# CountZIP — our worker=TRUE fix should help
{
  des = bern_des("count", y_cnt)
  inf = InferenceCountZeroInflatedPoisson$new(des, model_formula = ~ 1)
  timed("CountZeroInflatedPoisson boot_bayes/ci",
        inf$compute_bayesian_bootstrap_confidence_interval(B = B, na.rm = TRUE, show_progress = FALSE))
}
# IncidKKCondLogitPlusGLMMOneLik — no change from our code
{
  des = bern_des("incidence", y_inc)
  inf = tryCatch(InferenceIncidKKCondLogitPlusGLMMOneLik$new(des), error = function(e) { cat("  SKIP (ctor err):", conditionMessage(e), "\n"); NULL })
  if (!is.null(inf))
    timed("IncidKKCondLogitPlusGLMM boot_bayes/ci",
          inf$compute_bayesian_bootstrap_confidence_interval(B = B, na.rm = TRUE, show_progress = FALSE))
}

# ─────────────────────────────────────────────────────────────
cat("\n=== 2. boot_stud/pval + boot_stud/ci (currently skip_boot_stud_slow) ===\n")
# IncidRiskDiff
{
  des = bern_des("incidence", y_inc)
  inf = InferenceIncidRiskDiff$new(des, model_formula = ~ 1)
  timed("IncidRiskDiff boot_stud/pval",
        inf$compute_bootstrap_two_sided_pval(B = B, type = "studentized", na.rm = TRUE, show_progress = FALSE))
  timed("IncidRiskDiff boot_stud/ci",
        inf$compute_bootstrap_confidence_interval(B = B, type = "studentized", na.rm = TRUE, show_progress = FALSE))
}
# CountNegBin
{
  des = bern_des("count", y_cnt)
  inf = InferenceCountNegBin$new(des, model_formula = ~ 1)
  timed("CountNegBin boot_stud/pval",
        inf$compute_bootstrap_two_sided_pval(B = B, type = "studentized", na.rm = TRUE, show_progress = FALSE))
  timed("CountNegBin boot_stud/ci",
        inf$compute_bootstrap_confidence_interval(B = B, type = "studentized", na.rm = TRUE, show_progress = FALSE))
}
# OrdinalRidit
{
  des = bern_des("ordinal", y_ord)
  inf = InferenceOrdinalRidit$new(des)
  timed("OrdinalRidit boot_stud/pval",
        inf$compute_bootstrap_two_sided_pval(B = B, type = "studentized", na.rm = TRUE, show_progress = FALSE))
  timed("OrdinalRidit boot_stud/ci",
        inf$compute_bootstrap_confidence_interval(B = B, type = "studentized", na.rm = TRUE, show_progress = FALSE))
}

# ─────────────────────────────────────────────────────────────
cat("\n=== 3. boot_*/ci — pctile/basic/bca CI (currently skip_boot_ci_slow) ===\n")
# IncidMiettinenNurminenRiskDiff
{
  des = bern_des("incidence", y_inc)
  inf = tryCatch(InferenceIncidMiettinenNurminenRiskDiff$new(des), error = function(e) { cat("  SKIP (ctor err):", conditionMessage(e), "\n"); NULL })
  if (!is.null(inf)) {
    timed("MiettinenNurminenRiskDiff boot_pctile/ci",
          inf$compute_bootstrap_confidence_interval(B = B, na.rm = TRUE, show_progress = FALSE))
    timed("MiettinenNurminenRiskDiff boot_basic/ci",
          inf$compute_bootstrap_confidence_interval(B = B, type = "basic", na.rm = TRUE, show_progress = FALSE))
    timed("MiettinenNurminenRiskDiff boot_bca/ci",
          inf$compute_bootstrap_confidence_interval(B = B, type = "bca", na.rm = TRUE, show_progress = FALSE))
  }
}

# ─────────────────────────────────────────────────────────────
cat("\n=== 4. brt_*/ci + rand/ci (currently skip_rand_ci_slow) ===\n")
# PropKKGEE
{
  des = bern_des("proportion", y_prop)
  inf = tryCatch(InferencePropKKGEE$new(des, model_formula = ~ 1), error = function(e) { cat("  SKIP PropKKGEE (ctor err):", conditionMessage(e), "\n"); NULL })
  if (!is.null(inf)) {
    timed("PropKKGEE brt_pctile/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropKKGEE brt_stud/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, type = "studentized", pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropKKGEE brt_sym_t/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, type = "symmetric-percentile-t", pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropKKGEE brt_smooth/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, type = "smoothed", pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropKKGEE rand/ci",
          inf$compute_rand_confidence_interval(r = r, pval_epsilon = pval_eps, show_progress = FALSE))
  }
}
# PropBetaRegr
{
  des = bern_des("proportion", y_prop)
  inf = tryCatch(InferencePropBetaRegr$new(des, model_formula = ~ 1), error = function(e) { cat("  SKIP PropBetaRegr (ctor err):", conditionMessage(e), "\n"); NULL })
  if (!is.null(inf)) {
    timed("PropBetaRegr brt_pctile/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropBetaRegr brt_stud/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, type = "studentized", pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropBetaRegr brt_sym_t/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, type = "symmetric-percentile-t", pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropBetaRegr brt_smooth/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, type = "smoothed", pval_epsilon = pval_eps, show_progress = FALSE))
    timed("PropBetaRegr rand/ci",
          inf$compute_rand_confidence_interval(r = r, pval_epsilon = pval_eps, show_progress = FALSE))
  }
}

# ─────────────────────────────────────────────────────────────
# Reference: ContinRobustRegr brt_*/ci (structural SLOW — skip_brt="ci" in audit)
cat("\n=== 5. ContinRobustRegr brt_*/ci (reference — expect SLOW) ===\n")
{
  y_cont = rnorm(N)
  des = bern_des("continuous", y_cont)
  inf = tryCatch(InferenceContinRobustRegr$new(des, model_formula = ~ 1), error = function(e) { cat("  SKIP (ctor err):", conditionMessage(e), "\n"); NULL })
  if (!is.null(inf)) {
    timed("ContinRobustRegr brt_pctile/ci",
          inf$compute_rand_bootstrap_confidence_interval(B = B, pval_epsilon = pval_eps, show_progress = FALSE))
  }
}

cat("\nDone.\n")
