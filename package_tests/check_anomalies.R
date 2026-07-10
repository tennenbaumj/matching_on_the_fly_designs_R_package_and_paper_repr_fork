#!/usr/bin/env Rscript
# Anomaly checker for comprehensive_tests results CSV.
# Usage: Rscript check_anomalies.R [csv_file]

library(data.table)

args = commandArgs(trailingOnly = TRUE)
csv_file = if (length(args) >= 1) args[1] else {
	here = normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
	file.path(here, "comprehensive_tests_results_nc_4_continuous.csv")
}

cat("Reading:", csv_file, "\n")
dt = fread(csv_file)
cat(sprintf("Rows: %d  Reps: %s  Status ok: %d  Status error: %d\n",
	nrow(dt),
	paste(sort(unique(dt$rep)), collapse = ","),
	sum(dt$status == "ok", na.rm = TRUE),
	sum(dt$status == "error", na.rm = TRUE)))

issues = list()

# 1. p-values outside [0, 1]
pval_rows = dt[grepl("pval", function_run) & status == "ok"]
pval_rows[, v := suppressWarnings(as.numeric(result_1))]
bad_pval = pval_rows[is.finite(v) & (v < 0 | v > 1)]
if (nrow(bad_pval) > 0) {
	issues[["pval_out_of_range"]] = bad_pval[, .(rep, design, inference_class, function_run, result_1)]
	cat(sprintf("\nISSUE: %d p-values outside [0,1]:\n", nrow(bad_pval)))
	print(bad_pval[, .(rep, design, inference_class, function_run, result_1)])
}

# 2. CI lower > CI upper (reversed)
ci_rows = dt[grepl("confidence_interval", function_run) & status == "ok"]
ci_rows[, lo := suppressWarnings(as.numeric(result_1))]
ci_rows[, hi := suppressWarnings(as.numeric(result_2))]
bad_ci_order = ci_rows[is.finite(lo) & is.finite(hi) & lo > hi]
if (nrow(bad_ci_order) > 0) {
	issues[["reversed_ci"]] = bad_ci_order[, .(rep, design, inference_class, function_run, result_1, result_2)]
	cat(sprintf("\nISSUE: %d confidence intervals with lower > upper:\n", nrow(bad_ci_order)))
	print(bad_ci_order[, .(rep, design, inference_class, function_run, result_1, result_2)])
}

# 3. CI that does not contain true beta_T but is extremely wide (>100 units)
ci_rows[, width := hi - lo]
extreme_wide_ci = ci_rows[is.finite(width) & width > 100]
if (nrow(extreme_wide_ci) > 0) {
	issues[["extreme_wide_ci"]] = extreme_wide_ci[, .(rep, design, inference_class, function_run, result_1, result_2, width)]
	cat(sprintf("\nWARN: %d CIs with width > 100:\n", nrow(extreme_wide_ci)))
	print(extreme_wide_ci[, .(rep, design, inference_class, function_run, result_1, result_2, width)])
}

# 4. Estimates that are extreme (|estimate| > 100, for datasets where responses are O(1))
est_rows = dt[function_run == "compute_estimate" & status == "ok" & response_type %in% c("continuous","proportion")]
est_rows[, v := suppressWarnings(as.numeric(result_1))]
bad_est = est_rows[is.finite(v) & abs(v) > 100]
if (nrow(bad_est) > 0) {
	issues[["extreme_estimate"]] = bad_est[, .(rep, design, inference_class, function_run, result_1)]
	cat(sprintf("\nWARN: %d extreme estimates (|est| > 100):\n", nrow(bad_est)))
	print(bad_est[, .(rep, design, inference_class, function_run, result_1)])
}

# 5. Coverage: fraction of CIs that contain beta_T should be ~95% per inference class
if ("beta_T_in_confidence_interval" %in% names(dt)) {
	ci_cov = dt[!is.na(beta_T_in_confidence_interval) & grepl("confidence_interval", function_run) & status == "ok",
		.(coverage = mean(beta_T_in_confidence_interval, na.rm = TRUE), n = .N),
		by = .(inference_class, function_run)]
	bad_cov = ci_cov[n >= 10 & (coverage < 0.70 | coverage > 1.0)]
	if (nrow(bad_cov) > 0) {
		issues[["low_coverage"]] = bad_cov
		cat(sprintf("\nWARN: %d CI types with coverage < 70%% (>= 10 reps):\n", nrow(bad_cov)))
		print(bad_cov[order(coverage)])
	}
}

# 6. Error rates: >20% errors for a given inference_class + function_run combo
err_rate = dt[, .(n_err = sum(status == "error", na.rm = TRUE), n_total = .N, err_rate = mean(status == "error", na.rm = TRUE)),
	by = .(inference_class, function_run)]
high_err = err_rate[n_total >= 5 & err_rate > 0.20]
if (nrow(high_err) > 0) {
	issues[["high_error_rate"]] = high_err
	cat(sprintf("\nWARN: %d combos with >20%% error rate:\n", nrow(high_err)))
	print(high_err[order(-err_rate)][1:min(20, nrow(high_err))])
}

if (length(issues) == 0) {
	cat("\nAll checks passed — no anomalies detected.\n")
} else {
	cat(sprintf("\nSummary: %d anomaly categories found.\n", length(issues)))
}
