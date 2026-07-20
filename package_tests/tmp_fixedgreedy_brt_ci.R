#!/usr/bin/env Rscript
# Check BRT CI for PropKKGEE/PropBetaRegr with FixedGreedy design (design-draw-dependent)
suppressPackageStartupMessages({ library(EDI); library(dplyr); library(data.table) })
set_num_cores(1L); toggle_asserts(FALSE)
B = 151L; pval_eps = 0.05

timed = function(lbl, expr) {
  t0 = proc.time()
  res = tryCatch(suppressWarnings(expr), error = function(e) { cat("  ERR:", conditionMessage(e), "\n"); NULL })
  elapsed = (proc.time() - t0)[3]
  cat(sprintf("  [%s] %-60s %.1fs\n", if (elapsed < 30) "OK  " else "SLOW", lbl, elapsed))
  invisible(res)
}

set.seed(42); N = 148L
ds = ggplot2::diamonds |> slice_sample(n = N, replace = TRUE) |> mutate_if(where(is.factor), as.character)
X_f = ds |> model.matrix(price ~ 0 + ., data = _) |> apply(2, scale) |>
      as.data.frame() |> select(where(~ !any(is.na(.)))) |> as.data.frame()
X_f = X_f / ncol(X_f)
y_prop = pmin(0.99, pmax(0.01, (as.numeric(ds$price) - min(as.numeric(ds$price))) /
                               diff(range(as.numeric(ds$price)))))

greedy_des_prop = function(y) {
  des = DesignFixedGreedy$new(response_type = "proportion", n = N, design_formula = ~ .)
  des$add_all_subjects_to_experiment(X_f)
  des$assign_w_to_all_subjects()
  for (i in seq_len(N)) des$add_one_subject_response(i, y[i])
  des
}

cat("\n=== PropKKGEE + FixedGreedy ===\n")
des = greedy_des_prop(y_prop)
inf = tryCatch(InferencePropKKGEE$new(des, model_formula = ~ 1), error = function(e) { cat("ERR:", e$message, "\n"); NULL })
if (!is.null(inf)) {
  timed("PropKKGEE/FixedGreedy brt_pctile/ci",
        inf$compute_rand_bootstrap_confidence_interval(B = B, pval_epsilon = pval_eps, show_progress = FALSE))
  timed("PropKKGEE/FixedGreedy rand/ci",
        inf$compute_rand_confidence_interval(r = B, pval_epsilon = pval_eps, show_progress = FALSE))
}

cat("\n=== PropBetaRegr + FixedGreedy ===\n")
des = greedy_des_prop(y_prop)
inf = tryCatch(InferencePropBetaRegr$new(des, model_formula = ~ 1), error = function(e) { cat("ERR:", e$message, "\n"); NULL })
if (!is.null(inf)) {
  timed("PropBetaRegr/FixedGreedy brt_pctile/ci",
        inf$compute_rand_bootstrap_confidence_interval(B = B, pval_epsilon = pval_eps, show_progress = FALSE))
  timed("PropBetaRegr/FixedGreedy rand/ci",
        inf$compute_rand_confidence_interval(r = B, pval_epsilon = pval_eps, show_progress = FALSE))
}
cat("Done.\n")
