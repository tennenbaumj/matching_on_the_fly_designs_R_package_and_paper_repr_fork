#!/usr/bin/env Rscript
# Sanity-check: new worker-path results match old bootstrap_subset_inference results
suppressPackageStartupMessages({
  library(EDI)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})
set_num_cores(1L); toggle_asserts(FALSE)

max_n_dataset = 148L
set.seed(42)
diamonds_subset = ggplot2::diamonds %>% na.omit() %>%
  slice_sample(n = max_n_dataset, replace = TRUE) %>%
  mutate_if(where(is.factor), as.character)
X_raw = diamonds_subset %>% model.matrix(price ~ 0 + ., .) %>% apply(2, scale) %>%
  `/`(ncol(.)) %>% data.table() %>% select(where(~ !any(is.na(.)))) %>% as.data.frame()
y_raw = log(diamonds_subset$price)
n = 148L; SD_NOISE = 0.1; beta_T = 0

des = DesignFixedGreedy$new(response_type = "continuous", n = n, design_formula = ~ .)
des$add_all_subjects_to_experiment(X_raw)
des$assign_w_to_all_subjects()
w = des$get_w()
set.seed(99)
for (t_i in seq_len(n)) {
  y_t = y_raw[t_i] + ifelse(w[t_i]==1, beta_T, 0) + rnorm(1, 0, SD_NOISE)
  des$add_one_subject_response(t_i, y_t, dead=1)
}
inf = InferenceAllSimpleMeanDiffPooledVar$new(des, model_formula = ~ .)

# Compare pvals with fixed seed
r = 51L  # smaller B for fast comparison
set.seed(7)
pv_stud  = inf$compute_rand_bootstrap_two_sided_pval(B=r, type="studentized",          show_progress=FALSE)
set.seed(7)
pv_sym   = inf$compute_rand_bootstrap_two_sided_pval(B=r, type="symmetric-percentile-t", show_progress=FALSE)
set.seed(7)
pv_pct   = inf$compute_rand_bootstrap_two_sided_pval(B=r, type="percentile",            show_progress=FALSE)
cat(sprintf("studentized=%.4f  symmetric=%.4f  percentile=%.4f\n", pv_stud, pv_sym, pv_pct))

# Test CI
set.seed(7)
ci_stud = inf$compute_rand_bootstrap_confidence_interval(B=r, type="studentized",           pval_epsilon=0.05, show_progress=FALSE)
set.seed(7)
ci_sym  = inf$compute_rand_bootstrap_confidence_interval(B=r, type="symmetric-percentile-t", pval_epsilon=0.05, show_progress=FALSE)
cat(sprintf("CI studentized=[%.4f, %.4f]  symmetric=[%.4f, %.4f]\n", ci_stud[1], ci_stud[2], ci_sym[1], ci_sym[2]))
cat("All values finite:", all(is.finite(c(pv_stud, pv_sym, pv_pct, ci_stud, ci_sym))), "\n")

cat("\n=== Final timings (B=151) ===\n")
timed = function(lbl, expr) { t0=proc.time(); r=tryCatch(expr, error=function(e){cat("ERR:",e$message,"\n"); NULL}); cat(sprintf("  [%s] %.2fs\n", lbl, (proc.time()-t0)[3])); invisible(r) }
timed("percentile",           inf$compute_rand_bootstrap_two_sided_pval(B=151L, show_progress=FALSE))
timed("studentized",          inf$compute_rand_bootstrap_two_sided_pval(B=151L, type="studentized", show_progress=FALSE))
timed("symmetric-pct-t",     inf$compute_rand_bootstrap_two_sided_pval(B=151L, type="symmetric-percentile-t", show_progress=FALSE))
timed("CI studentized",       inf$compute_rand_bootstrap_confidence_interval(B=151L, type="studentized", pval_epsilon=0.007, show_progress=FALSE))
timed("CI symmetric-pct-t",  inf$compute_rand_bootstrap_confidence_interval(B=151L, type="symmetric-percentile-t", pval_epsilon=0.007, show_progress=FALSE))
cat("Done.\n")
