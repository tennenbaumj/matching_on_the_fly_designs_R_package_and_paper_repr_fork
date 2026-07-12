rm(list = ls())
set.seed(1)

script_file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(script_file_arg)) {
	normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = TRUE)
} else {
	normalizePath("package_tests/comprehensive_tests.R", mustWork = TRUE)
}
script_dir = dirname(script_path)
repo_root = normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
repo_path = function(...) file.path(repo_root, ...)

pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis, devtools)
suppressPackageStartupMessages(library(EDI))

welch_t_stat_cpp = '
	double welch_t_stat(NumericVector y, IntegerVector w) {
		int n = y.size();
		double sum_t = 0, sum_c = 0;
		int n_t = 0, n_c = 0;
		for (int i = 0; i < n; i++) {
			if (w[i] == 1) { sum_t += y[i]; n_t++; }
			else           { sum_c += y[i]; n_c++; }
		}
		double mean_t = sum_t / n_t;
		double mean_c = sum_c / n_c;
		double var_t = 0, var_c = 0;
		for (int i = 0; i < n; i++) {
			if (w[i] == 1) { double d = y[i] - mean_t; var_t += d * d; }
			else           { double d = y[i] - mean_c; var_c += d * d; }
		}
		var_t /= (n_t - 1);
		var_c /= (n_c - 1);
		return (mean_t - mean_c) / sqrt(var_t / n_t + var_c / n_c);
	}
'

source(repo_path("EDI", "tests", "testthat", "helper-likelihood-method-smoke.R"))

max_n_dataset = 148 #needs to be divisible by 4 for some blocking designs
source(repo_path("package_tests", "_dataset_load.R"))
# options(error = recover)
# options(warn=2)

args = commandArgs(trailingOnly = TRUE)
Nrep      = if (length(args) >= 1) as.integer(args[1]) else 40L
NUM_CORES = if (length(args) >= 2) as.integer(args[2]) else 2L
ALL_RESPONSE_TYPES = c("continuous", "incidence", "proportion", "count", "survival", "ordinal")
RESPONSE_TYPE_FILTER = if (length(args) >= 3 && args[3] != "NA") as.character(args[3]) else NA_character_
if (!is.na(RESPONSE_TYPE_FILTER) && !(RESPONSE_TYPE_FILTER %in% ALL_RESPONSE_TYPES)) {
	stop(
		"Unsupported response_type filter: ",
		RESPONSE_TYPE_FILTER,
		". Supported values are: ",
		paste(ALL_RESPONSE_TYPES, collapse = ", ")
	)
}
ALL_DESIGN_TYPES = c("Bernoulli", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise", "SPBR", "PocockSimon", "Urn", "RandomBlockSize", "FixedBernoulli", "FixediBCRD", "FixedBlocking", "FixedCluster", "FixedBlockedCluster", "FixedBinaryMatch", "FixedGreedy", "FixedRerandomization", "FixedMatchingGreedy", "FixedDOptimal", "FixedAOptimal")
DESIGN_TYPE_FILTER = if (length(args) >= 4) as.character(args[4]) else NA_character_
if (!is.na(DESIGN_TYPE_FILTER) && !(DESIGN_TYPE_FILTER %in% ALL_DESIGN_TYPES)) {
	stop(
		"Unsupported design_type filter: ",
		DESIGN_TYPE_FILTER,
		". Supported values are: ",
		paste(ALL_DESIGN_TYPES, collapse = ", ")
	)
}
INFERENCE_CLASS_FILTER = if (length(args) >= 5 && args[5] != "NA") as.character(args[5]) else NA_character_
DATASET_FILTER = if (length(args) >= 6 && args[6] != "NA") as.character(args[6]) else NA_character_
BETA_T_FILTER = if (length(args) >= 7 && args[7] != "NA") as.numeric(args[7]) else NA_real_
REP_FILTER = if (length(args) >= 8 && args[8] != "NA") as.integer(args[8]) else NA_integer_
canonicalize_test_family_filter = function(value){
	value = as.character(value)
	switch(
		value,
		lr = "lik_ratio",
		likelihood_ratio = "lik_ratio",
		bayes_bootstrap = "bayesian_bootstrap",
		randomization = "rand",
		value
	)
}
ALL_TEST_FAMILY_FILTERS = c("estimate", "exact", "asymp", "wald", "score", "lik_ratio", "gradient", "bootstrap", "bayesian_bootstrap", "parametric_bootstrap", "jackknife", "rand", "rand_custom")
TEST_FAMILY_FILTER = if (length(args) >= 9 && args[9] != "NA") canonicalize_test_family_filter(args[9]) else NA_character_
if (!is.na(TEST_FAMILY_FILTER) && !(TEST_FAMILY_FILTER %in% ALL_TEST_FAMILY_FILTERS)) {
	stop(
		"Unsupported test family filter: ",
		TEST_FAMILY_FILTER,
		". Supported values are: ",
		paste(ALL_TEST_FAMILY_FILTERS, collapse = ", ")
	)
}
set_num_cores(NUM_CORES)
toggle_asserts(FALSE)
if (is.na(INFERENCE_CLASS_FILTER)) {
	run_likelihood_method_smoke_suite(RESPONSE_TYPE_FILTER)
} else {
	message("Skipping likelihood method smoke suite because INFERENCE_CLASS_FILTER is set.")
}

prob_censoring = 0.15
r = 151
pval_epsilon = 0.007
test_compute_confidence_interval_rand = TRUE
run_debug_resampling = Sys.getenv("COMPREHENSIVE_DEBUG_RESAMPLING", "0") %in% c("1", "true", "TRUE", "yes", "YES")
run_parametric_bootstrap_ci = Sys.getenv("COMPREHENSIVE_PARAM_BOOT_CI", "0") %in% c("1", "true", "TRUE", "yes", "YES")
param_boot_ci_max_root_iterations = 0L
beta_T_values = c(0, 1)
SD_NOISE = 0.1
pending_rep_header = NULL
pending_beta_header = NULL
pending_dataset_header = NULL
pending_response_header = NULL
pending_design_header = NULL
pending_banner = NULL

sanitize_results_filter = function(value){
	value = as.character(value)
	value = gsub("[^A-Za-z0-9._=-]+", "_", value)
	value = gsub("_+", "_", value)
	gsub("^_|_$", "", value)
}

extra_filter_parts = character(0)
if (!is.na(DESIGN_TYPE_FILTER)) extra_filter_parts = c(extra_filter_parts, paste0("design-", sanitize_results_filter(DESIGN_TYPE_FILTER)))
if (!is.na(INFERENCE_CLASS_FILTER)) extra_filter_parts = c(extra_filter_parts, paste0("class-", sanitize_results_filter(INFERENCE_CLASS_FILTER)))
if (!is.na(DATASET_FILTER)) extra_filter_parts = c(extra_filter_parts, paste0("dataset-", sanitize_results_filter(DATASET_FILTER)))
if (!is.na(BETA_T_FILTER)) extra_filter_parts = c(extra_filter_parts, paste0("beta-", sanitize_results_filter(BETA_T_FILTER)))
if (!is.na(REP_FILTER)) extra_filter_parts = c(extra_filter_parts, paste0("rep-", sanitize_results_filter(REP_FILTER)))
if (!is.na(TEST_FAMILY_FILTER)) extra_filter_parts = c(extra_filter_parts, paste0("family-", sanitize_results_filter(TEST_FAMILY_FILTER)))
filtered_results_suffix = if (length(extra_filter_parts)) {
	paste0("_filtered_", paste(extra_filter_parts, collapse = "_"))
} else {
	""
}

results_file = if (is.na(RESPONSE_TYPE_FILTER)) {
	repo_path("package_tests", paste0("comprehensive_tests_results_nc_", NUM_CORES, filtered_results_suffix, ".csv"))
} else {
	repo_path("package_tests", paste0("comprehensive_tests_results_nc_", NUM_CORES, "_", RESPONSE_TYPE_FILTER, filtered_results_suffix, ".csv"))
}
existing_results_dt = if (file.exists(results_file)) data.table::fread(results_file) else data.table::data.table()
run_row_id = if ("run_row_id" %in% colnames(existing_results_dt) && nrow(existing_results_dt) > 0L) {
	as.integer(max(existing_results_dt$run_row_id, na.rm = TRUE))
} else {
	0L
}

serialize_beta_T = function(value){
	if (is.na(value)) return("NA")
	if (is.numeric(value)) return(format(as.numeric(value), scientific = TRUE, digits = 17))
	as.character(value)
}

round_result_field = function(value){
	if (is.null(value) || is.na(value)) return(NA_character_)
	if (value == "") return("")
	num = suppressWarnings(as.numeric(value))
	if (!is.finite(num)) return(value)
	sprintf("%.3f", num)
}

round_duration_field = function(value){
	if (is.null(value) || is.na(value)) return(NA_real_)
	num = suppressWarnings(as.numeric(value))
	if (!is.finite(num)) return(num)
	round(num, 3)
}

add_assignment_only_cluster_id = function(X_design, strata_cols = character(0), cluster_size = 2L){
	X_out = as.data.frame(X_design)
	cluster_col = ".assignment_only_cluster_id"
	cluster_ids = integer(nrow(X_out))
	next_cluster_id = 1L

	if (!length(strata_cols)) {
		cluster_ids = ((seq_len(nrow(X_out)) - 1L) %/% cluster_size) + 1L
	} else {
		strata_df = X_out[, strata_cols, drop = FALSE]
		strata_key = do.call(paste, c(lapply(strata_df, as.character), sep = "\r"))
		for (key in unique(strata_key)) {
			idx = which(strata_key == key)
			cluster_ids[idx] = ((seq_along(idx) - 1L) %/% cluster_size) + next_cluster_id
			next_cluster_id = max(cluster_ids[idx]) + 1L
		}
	}

	X_out[[cluster_col]] = factor(cluster_ids)
	list(X = X_out, cluster_col = cluster_col)
}

build_result_key = function(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val){
	paste(
		as.integer(rep_val),
		serialize_beta_T(beta_val),
		as.character(dataset_val),
		as.character(response_val),
		as.character(design_val),
		as.character(inference_val),
		as.character(function_run_val),
		sep = "||"
	)
}

completed_rows_cache = new.env(parent = emptyenv())

mark_row_completed = function(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val){
	key = build_result_key(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val)
	completed_rows_cache[[key]] <- TRUE
}

last_results_mtime = NULL

reload_completed_rows = function() {
	if (file.exists(results_file) && file.info(results_file)$size > 0) {
		current_mtime = file.info(results_file)$mtime
		if (!is.null(last_results_mtime) && current_mtime == last_results_mtime) {
			return(invisible(NULL))
		}
		
		# Reset the cache
		completed_rows_cache <<- new.env(parent = emptyenv())
		
		# Define columns needed for building the key and filtering status
		needed_cols = c("rep", "beta_T", "dataset", "response_type", "design", "inference_class", "function_run", "status")
		
		# Try to read the file, handling potential locks or partial writes with retries
		max_retries = 10
		retry_count = 0
		dt = NULL
		
		while (retry_count < max_retries) {
			dt = tryCatch({
				# Check which columns actually exist first
				header = names(data.table::fread(results_file, nrows = 0))
				cols_to_read = intersect(needed_cols, header)
				data.table::fread(results_file, select = cols_to_read)
			}, error = function(e) {
				NULL
			})
			
			if (!is.null(dt)) break
			
			retry_count = retry_count + 1
			if (retry_count < max_retries) {
				Sys.sleep(0.5)
			}
		}
		
		if (is.null(dt)) {
			# If still failing after retries, return early (cache remains empty)
			return(invisible(NULL))
		}
		
		last_results_mtime <<- current_mtime
		
		if (nrow(dt) > 0L) {
			rows_to_cache = if ("status" %in% colnames(dt)) {
				dt[status == "ok"]
			} else {
				dt
			}
			for (row_idx in seq_len(nrow(rows_to_cache))) {
				row = rows_to_cache[row_idx]
				mark_row_completed(
					row$rep,
					row$beta_T,
					row$dataset,
					row$response_type,
					row$design,
					row$inference_class,
					row$function_run
				)
			}
		}
	}
}

is_row_completed = function(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val){
	# Reload from disk each time we check
	reload_completed_rows()
	
	key = build_result_key(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val)
	!is.null(completed_rows_cache[[key]])
}

if (nrow(existing_results_dt) > 0L) {
	rows_to_cache = if ("status" %in% colnames(existing_results_dt)) {
		existing_results_dt[status == "ok"]
	} else {
		existing_results_dt
	}
	for (row_idx in seq_len(nrow(rows_to_cache))) {
		row = rows_to_cache[row_idx]
		mark_row_completed(
			row$rep,
			row$beta_T,
			row$dataset,
			row$response_type,
			row$design,
			row$inference_class,
			row$function_run
		)
	}
}
results_dt = data.table(
	rep = integer(),
	beta_T = numeric(),
	dataset = character(),
	response_type = character(),
	design = character(),
	inference_class = character(),
	function_run = character(),
	timestamp = character(),
	duration_time_sec = numeric(),
	result_1 = character(),
	result_2 = character(),
	beta_T_in_confidence_interval = logical(),
	error_message = character(),
	run_row_id = integer(),
	r = integer(),
	pval_epsilon = numeric(),
	prob_censoring = numeric(),
	sd_noise = numeric(),
	num_cores = integer(),
	dataset_n_rows = integer(),
	dataset_n_cols = integer(),
	result = character(),
	status = character()
)
write_results_if_needed = function(force = FALSE){
	if ((force || nrow(results_dt) > 0L) && nrow(results_dt) > 0L){
		append_mode = file.exists(results_file) && file.info(results_file)$size > 0
		data.table::fwrite(
			results_dt,
			results_file,
			append = append_mode,
			col.names = !append_mode,
			na = "NA"
		)
		results_dt <<- results_dt[0]
	}
}

log_progress = function(msg){
	message(msg)
	flush.console()
}

is_skipped_inference_label = function(inference_label){
	grepl("IVWC", inference_label, fixed = TRUE)
}

should_run_inference_label = function(inference_label){
	!is_skipped_inference_label(inference_label) &&
		(is.na(INFERENCE_CLASS_FILTER) || grepl(INFERENCE_CLASS_FILTER, inference_label))
}

is_test_family_filter_active = function(){
	!is.na(TEST_FAMILY_FILTER)
}

should_run_test_family = function(test_family){
	!is_test_family_filter_active() || identical(TEST_FAMILY_FILTER, test_family)
}

inference_banner = function(inf_name, mf = NULL){
	if (is.null(mf)) {
		mf = if (exists("model_formula", envir = .GlobalEnv)) get("model_formula", envir = .GlobalEnv) else NULL
	}
	pending_banner <<- sprintf("\n\n  == Inference: %s design_type = %s dataset = %s response_type = %s beta_T = [%s] num_cores = [%d] rep = [%d/%d]%s\n", 
	  inf_name, 
	  if (exists("design_type", envir = .GlobalEnv)) get("design_type", envir = .GlobalEnv) else "unknown",
	  if (exists("dataset_name", envir = .GlobalEnv)) get("dataset_name", envir = .GlobalEnv) else "unknown",
	  if (exists("response_type", envir = .GlobalEnv)) get("response_type", envir = .GlobalEnv) else "unknown",
	  format(beta_T), NUM_CORES, rep_curr, Nrep,
	  if (!is.null(mf)) paste0(" formula = [", deparse(mf), "]") else "")
}

record_result = function(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_class, function_run, result, status, duration_time_sec, error_message = NA_character_){
	result_vec = if (is.null(result)) {
		NA_character_
	} else if (length(result) == 0) {
		character(0)
	} else {
		as.character(result)
	}

	result_str = if (is.null(result)) {
		NA_character_
	} else if (length(result) == 0) {
		""
	} else if (length(result) == 1) {
		as.character(result)
	} else {
		paste(as.character(result), collapse = " ")
	}
	result_1 = if (length(result_vec) >= 1) round_result_field(result_vec[1]) else NA_character_
	is_debug_distribution = grepl("_debug$", function_run)
	result_2 = if ((grepl("confidence_interval", function_run, fixed = TRUE) || is_debug_distribution) && length(result_vec) >= 2) round_result_field(result_vec[2]) else NA_character_
	beta_T_in_confidence_interval = NA
	if (grepl("confidence_interval", function_run, fixed = TRUE) && length(result) >= 2 && all(is.finite(result[1:2]))){
		ci_lo = min(result[1:2])
		ci_hi = max(result[1:2])
		beta_T_in_confidence_interval = (beta_T >= ci_lo && beta_T <= ci_hi)
	}
	run_row_id <<- run_row_id + 1L
	results_dt <<- data.table::rbindlist(list(
		results_dt,
		data.table(
			rep = as.integer(rep_curr),
			beta_T = beta_T,
			dataset = dataset_name,
			response_type = response_type,
			design = design_type,
			inference_class = inference_class,
			function_run = function_run,
			timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
			duration_time_sec = round_duration_field(duration_time_sec),
			result_1 = result_1,
			result_2 = result_2,
			beta_T_in_confidence_interval = beta_T_in_confidence_interval,
			error_message = ifelse(is.null(error_message), NA_character_, as.character(error_message)),
			run_row_id = run_row_id,
			r = as.integer(r),
			pval_epsilon = pval_epsilon,
			prob_censoring = prob_censoring,
			sd_noise = SD_NOISE,
			num_cores = as.integer(NUM_CORES),
			dataset_n_rows = as.integer(dataset_n_rows),
			dataset_n_cols = as.integer(dataset_n_cols),
			result = result_str,
			status = status
		)
	), use.names = TRUE)
	if (identical(status, "ok")) {
		mark_row_completed(rep_curr, beta_T, dataset_name, response_type, design_type, inference_class, function_run)
	}
	write_results_if_needed(force = TRUE)
}

record_existing_error_keys_as_skipped = function(inference_class, dataset_n_rows, dataset_n_cols, error_message){
	if (!file.exists(results_file) || file.info(results_file)$size <= 0) return(invisible(NULL))
	dt = tryCatch(data.table::fread(results_file), error = function(e) NULL)
	if (is.null(dt) || nrow(dt) == 0L) return(invisible(NULL))
	key_cols = c("rep", "beta_T", "dataset", "response_type", "design", "inference_class", "function_run")
	if (!all(c(key_cols, "status") %in% names(dt))) return(invisible(NULL))
	latest = dt[, .SD[.N], by = key_cols]
	to_mark = latest[
		status == "error" &
			rep == as.integer(rep_curr) &
			as.numeric(beta_T) == as.numeric(get("beta_T", envir = .GlobalEnv)) &
			dataset == dataset_name &
			response_type == get("response_type", envir = .GlobalEnv) &
			design == design_type &
			inference_class == inference_class,
		unique(function_run)
	]
	if (!length(to_mark)) return(invisible(NULL))
	for (fn in to_mark) {
		record_result(
			dataset_name, dataset_n_rows, dataset_n_cols,
			get("response_type", envir = .GlobalEnv), design_type,
			inference_class, fn, NA_character_, status = "ok",
			duration_time_sec = 0,
			error_message = error_message
		)
	}
	invisible(NULL)
}

run_inference_checks_both_paths = function(class_gen, class_name, des_obj, response_type, design_type, dataset_name, n_rows, n_cols, model_formula = NULL, ...){
	# Univariate path
	univ_label = paste0(class_name, " (Univ)")
	if (should_run_inference_label(univ_label)) {
		inference_banner(univ_label, mf = ~ 1)
		run_inference_checks(
			class_gen$new(des_obj, model_formula = ~ 1, ...),
			response_type, design_type, dataset_name, n_rows, n_cols,
			inference_label = univ_label
		)
	}
	
	# Multivariate path
	multi_label = paste0(class_name, " (Multi)")
	if (should_run_inference_label(multi_label)) {
		inference_banner(multi_label, mf = ~ .)
		run_inference_checks(
			class_gen$new(des_obj, model_formula = ~ ., ...),
			response_type, design_type, dataset_name, n_rows, n_cols,
			inference_label = multi_label
		)
	}
}

run_inference_checks_for_class = function(class_gen, class_name, des_obj, response_type, design_type, dataset_name, n_rows, n_cols, model_formula = NULL, ...){
	if (!should_run_inference_label(class_name)) return(invisible(NULL))
	inference_banner(class_name, mf = model_formula)
	run_inference_checks(
		class_gen$new(des_obj, model_formula = model_formula, ...),
		response_type, design_type, dataset_name, n_rows, n_cols
	)
}

run_inference_checks = function(seq_des_inf, response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols, exhaustive_sweep = FALSE, inference_label = NULL){
	inference_base_label = if (!is.null(inference_label) && nzchar(inference_label)) {
		inference_label
	} else {
		class(seq_des_inf)[1]
	}
	design_formula_suffix = if (exists(".comprehensive_current_model_formula", envir = .GlobalEnv, inherits = FALSE)) {
		current_model_formula = get(".comprehensive_current_model_formula", envir = .GlobalEnv)
		paste0(" [design_formula=", paste(deparse(current_model_formula), collapse = " "), "]")
	} else {
		""
	}
	inference_result_label = paste0(inference_base_label, design_formula_suffix)
	if (!should_run_inference_label(inference_result_label)) {
		return(invisible(NULL))
	}
	skip_slow = exhaustive_sweep
	B_debug = as.integer(r)
	r_debug = as.integer(r)
	is_any_inference_class = function(classes){
		any(vapply(classes, function(cls) is(seq_des_inf, cls), logical(1)))
	}
	skip_bootstrap = is_any_inference_class(c(
		"InferenceAbstractKKGEE",
		"InferenceCountPoissonKKGEE",
		"InferenceAbstractKKGLMM",
		"InferenceCountKKGLMM",
		"InferenceContinMultGLS",
		"InferenceAbstractKKClaytonCopulaIVWC",
		"InferenceSurvivalKKClaytonCopulaOneLik",
		"InferenceAbstractKKWeibullFrailtyIVWC",
		"InferenceAbstractKKWeibullFrailtyOneLik",
		"InferenceAllKKWilcoxIVWC",
		"InferenceAbstractKKWilcoxRegrIVWC",
		"InferenceSurvivalKKRankRegrIVWC",
		"InferenceIncidExactZhangAbstract",
		"InferenceAllSimpleWilcox",
		"InferenceOrdinalPairedSignTest",
		"InferenceOrdinalKKCondAdjCatLogitRegr",
		"InferenceOrdinalGCompMeanDiff",
		"InferenceOrdinalGCompMeanDiff",
		"InferenceOrdinalCloglogRegr",
		"InferenceOrdinalOrderedProbitRegr",
		"InferenceOrdinalOrderedProbitRegr",
		"InferenceOrdinalCauchitRegr",
		"InferenceOrdinalCauchitRegr",
		"InferenceOrdinalKKCondAdjCatLogitRegr",
		"InferenceOrdinalPartialProportionalOddsRegr",
		"InferencePropZeroOneInflatedBetaRegr",
		"InferencePropFractionalLogit",
		"InferenceCountHurdleNegBin",
		"InferenceContinRobustRegr"
	))
	skip_bayesian_bootstrap = skip_bootstrap || is_any_inference_class(c(
		"InferenceBaiAdjustedT",
		"InferenceBaiAdjustedTKK14",
		"InferenceBaiAdjustedTKK21",
		"InferenceIncidKKCondLogitIVWC",
		"InferenceAbstractKKHurdlePoissonIVWC",
		"InferenceAbstractKKStratCoxIVWC",
		"InferenceAbstractKKWeibullFrailtyIVWC",
		"InferenceAbstractKKClaytonCopulaIVWC",
		"InferenceAbstractKKLWACoxIVWC",
		"InferenceAbstractKKSurvivalRankRegrIVWC",
		"InferenceAbstractKKRobustRegrIVWC",
		"InferenceAbstractKKQuantileRegrIVWC",
		"InferenceAbstractKKWilcoxBaseIVWC",
		"InferenceAllKKMeanDiffIVWC",
		"InferenceAllKKWilcoxIVWC",
		"InferenceIncidKKCondLogitPlusGLMMIVWC",
		"InferenceContinKKOLSIVWC",
		"InferenceContinKKRobustRegrIVWC",
		"InferenceContinKKQuantileRegrIVWC",
		"InferencePropKKQuantileRegrIVWC",
		"InferenceSurvivalKKRankRegrIVWC"
	))
	supports_parametric_bootstrap =
		is(seq_des_inf, "InferenceParamBootstrap") &&
		!is_any_inference_class(c("InferenceCountKKGLMM")) &&
		isTRUE(tryCatch(
			seq_des_inf$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(),
			error = function(e) FALSE
		))
	supports_parametric_bootstrap_ci =
		run_parametric_bootstrap_ci &&
		supports_parametric_bootstrap &&
		!is_any_inference_class(c("InferenceCountKKGLMM")) &&
		isTRUE(tryCatch({
			ci_support_fn = seq_des_inf$.__enclos_env__$private$supports_lik_ratio_param_bootstrap_confidence_interval
			if (is.function(ci_support_fn)) ci_support_fn() else TRUE
		}, error = function(e) TRUE))
	skip_rand      = is(seq_des_inf, "InferenceAbstractKKGEE") || is(seq_des_inf, "InferenceAbstractKKGLMM") || is(seq_des_inf, "InferenceIncidExactZhangAbstract") || is(seq_des_inf, "InferencePropGCompMeanDiff") || is(seq_des_inf, "InferencePropGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalPairedSignTest") || is(seq_des_inf, "InferenceOrdinalKKCondAdjCatLogitRegr") || is(seq_des_inf, "InferenceOrdinalGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalCloglogRegr") || is(seq_des_inf, "InferenceOrdinalOrderedProbitRegr") || is(seq_des_inf, "InferenceOrdinalOrderedProbitRegr") || is(seq_des_inf, "InferenceOrdinalCauchitRegr") || is(seq_des_inf, "InferenceOrdinalCauchitRegr") || is(seq_des_inf, "InferenceOrdinalKKCondAdjCatLogitRegr")
	skip_mle_pval  = is(seq_des_inf, "InferenceSurvivalKKWeibullFrailtyOneLik")
	skip_rand_pval = is(seq_des_inf, "InferenceSurvivalKKWeibullFrailtyOneLik") || is(seq_des_inf, "InferenceContinMultGLS") || is(seq_des_inf, "InferencePropGCompMeanDiff") || is(seq_des_inf, "InferencePropGCompMeanDiff") || is_any_inference_class(c(
		"InferenceSurvivalKKRankRegrIVWC",
		"InferenceSurvivalKKClaytonCopulaIVWC",
		"InferenceSurvivalKKClaytonCopulaOneLik"
	))
	skip_ci_rand   = is_any_inference_class(c(
		"InferenceContinMultKKQuantileRegrIVWC",
		"InferenceContinMultKKQuantileRegrOneLik",
		"InferencePropGCompMeanDiff",
		"InferencePropGCompMeanDiff",
		"InferenceContinRobustRegr",
		"InferenceCountHurdleNegBin",
		"InferenceCountNegBin",
		"InferenceCountZeroInflatedNegBin",
		"InferencePropZeroOneInflatedBetaRegr",
		"InferencePropFractionalLogit",
		"InferenceCountHurdlePoisson",
		"InferenceCountZeroInflatedPoisson",
		"InferenceCountHurdlePoisson",
		"InferenceCountZeroInflatedPoisson",
		"InferenceCountZeroInflatedNegBin",
		"InferenceCountKKHurdlePoissonOneLik"
	)) || response_type == "count" ||
		(response_type != "continuous" && (is(seq_des_inf, "InferenceAllSimpleMeanDiff") || is(seq_des_inf, "InferenceAllKKMeanDiffIVWC")))
	skip_ci_rand_custom = FALSE
	supports_jackknife = is(seq_des_inf, "InferenceJackknife") ||
		(
			"compute_jackknife_wald_two_sided_pval" %in% names(seq_des_inf) &&
			"compute_jackknife_wald_confidence_interval" %in% names(seq_des_inf)
		)
	supports_jackknife = supports_jackknife && !is_any_inference_class(c("InferenceCountKKGLMM"))
	
	skip_ci = beta_T == 1 && (
		is(seq_des_inf, "InferenceIncidLogRegr") ||
		is(seq_des_inf, "InferencePropBetaRegr") ||
		is(seq_des_inf, "InferencePropBetaRegr") ||
		is(seq_des_inf, "InferenceSurvivalCoxPHRegr") ||
		is(seq_des_inf, "InferenceSurvivalCoxPHRegr") ||
		is(seq_des_inf, "InferenceSurvivalKKLWACoxPHIVWC") ||
		is(seq_des_inf, "InferenceSurvivalKKStratCoxPHIVWC") ||
		is(seq_des_inf, "InferenceSurvivalKKClaytonCopulaIVWC") ||
		is(seq_des_inf, "InferenceSurvivalKKLWACoxPHOneLik") ||
		is(seq_des_inf, "InferenceSurvivalKKStratCoxPHOneLik") ||
		is(seq_des_inf, "InferenceSurvivalKKClaytonCopulaOneLik") ||
		is(seq_des_inf, "InferenceSurvivalKKWeibullFrailtyIVWC") ||
		is(seq_des_inf, "InferenceSurvivalKKWeibullFrailtyOneLik") ||
		is(seq_des_inf, "InferenceSurvivalKKRankRegrIVWC")
	)
	snap_small_numeric_to_zero = function(x, tol = sqrt(.Machine$double.eps)){
		if (is.null(x)) return(x)
		if (is.list(x)) return(lapply(x, snap_small_numeric_to_zero, tol = tol))
		if (is.atomic(x) && is.numeric(x)){
			x[is.finite(x) & abs(x) < tol] = 0
			return(x)
		}
		x
	}

	has_invalid_numeric = function(x){
		if (is.null(x)) return(FALSE)
		if (is.list(x)) return(any(vapply(x, has_invalid_numeric, logical(1))))
		if (is.atomic(x) && is.numeric(x)) return(any(!is.finite(x) | is.na(x) | is.nan(x)))
		FALSE
	}

	is_zero_zero_confidence_interval = function(label, result){
		if (!grepl("confidence_interval", label, fixed = TRUE)) return(FALSE)
		if (!(is.atomic(result) && is.numeric(result))) return(FALSE)
		if (length(result) < 2) return(FALSE)
		isTRUE(all(result[1:2] == 0))
	}

	is_allowed_missing_output = function(label, result){
		if (!has_invalid_numeric(result)) return(FALSE)
		identical(response_type, "ordinal") &&
			identical(label, "compute_asymp_two_sided_pval")
	}

is_explicitly_nonestimable = function(obj){
	if (is.null(obj) || !is.function(obj$is_nonestimable)) return(FALSE)
	isTRUE(tryCatch(obj$is_nonestimable(), error = function(e) FALSE))
}

supports_direct_testing_type = function(testing_type){
	if (is.null(seq_des_inf) || !is.function(seq_des_inf$get_supported_testing_types)) return(FALSE)
	supported = tryCatch(seq_des_inf$get_supported_testing_types(), error = function(e) character())
	testing_type %in% supported
}

	should_record_nonestimable_as_missing = function(obj, label, result = NULL){
		if (!is_explicitly_nonestimable(obj)) return(FALSE)
		# If the method returned a valid (finite) result despite non-estimable state, record it normally
		if (!is.null(result) && !has_invalid_numeric(result)) return(FALSE)
		stage = if (!is.null(obj) && is.function(obj$get_nonestimable_stage)) {
			tryCatch(obj$get_nonestimable_stage(), error = function(e) NULL)
		} else {
			NULL
		}
		if (identical(stage, "estimate")) return(TRUE)
		if (identical(stage, "se")) {
			# These methods require an original or per-replicate standard error.
			# If the inference object explicitly marked the SE unavailable, a
			# missing result is expected rather than a numerical-output failure.
			return(
				grepl("studentized", label, fixed = TRUE) ||
				grepl("bca", label, fixed = TRUE) ||
				identical(label, "compute_bootstrap_confidence_interval") ||
				identical(label, "compute_bootstrap_two_sided_pval") ||
				identical(label, "compute_bayesian_bootstrap_confidence_interval") ||
				identical(label, "compute_bayesian_bootstrap_two_sided_pval") ||
				grepl("jackknife", label, fixed = TRUE) ||
				grepl("wald", label, ignore.case = TRUE) ||
				grepl("score", label, ignore.case = TRUE) ||
				grepl("gradient", label, ignore.case = TRUE) ||
				grepl("lik_ratio", label, fixed = TRUE) ||
				grepl("asymp", label, ignore.case = TRUE) ||
				grepl("lik_ratio_bootstrap", label, fixed = TRUE) ||
				grepl("rand", label, ignore.case = TRUE)
			)
		}
		TRUE
	}

	nonestimable_error_message = function(obj, label){
		reason = if (!is.null(obj) && is.function(obj$get_nonestimable_reason)) {
			tryCatch(obj$get_nonestimable_reason(), error = function(e) NULL)
		} else {
			NULL
		}
		stage = if (!is.null(obj) && is.function(obj$get_nonestimable_stage)) {
			tryCatch(obj$get_nonestimable_stage(), error = function(e) NULL)
		} else {
			NULL
		}
		parts = c(
			paste0("Explicitly non-estimable in ", label),
			if (!is.null(stage) && nzchar(stage)) paste0("stage=", stage) else NULL,
			if (!is.null(reason) && nzchar(reason)) paste0("reason=", reason) else NULL
		)
		paste(parts, collapse = "; ")
	}

safe_call = function(label, expr){
	if (is_row_completed(
		rep_curr,
		beta_T,
			dataset_name,
			response_type,
			design_type,
			inference_result_label,
			label
		)) {
			return(invisible(NULL))
	}
	
	if (!is.null(pending_rep_header)) { message(pending_rep_header); pending_rep_header <<- NULL }
	if (!is.null(pending_beta_header)) { message(pending_beta_header); pending_beta_header <<- NULL }
	if (!is.null(pending_dataset_header)) { message(pending_dataset_header); pending_dataset_header <<- NULL }
	if (!is.null(pending_response_header)) { message(pending_response_header); pending_response_header <<- NULL }
	if (!is.null(pending_design_header)) { message(pending_design_header); pending_design_header <<- NULL }
	if (!is.null(pending_banner)){
		message(pending_banner)
		pending_banner <<- NULL
	}

	message("          Calling ", label, "()")
		start_elapsed = unname(proc.time()[["elapsed"]])
		tryCatch({
			result <- expr
			if (should_record_nonestimable_as_missing(seq_des_inf, label, result)) {
				msg = nonestimable_error_message(seq_des_inf, label)
				message("Recording missing output for ", label, " as ok (explicitly non-estimable).")
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_result_label, label, NA_character_, status = "ok", duration_time_sec = duration_time_sec, error_message = msg)
				return(invisible(NULL))
			}
			if (is_allowed_missing_output(label, result)) {
				message("Recording missing output for ", label, " as ok (ordinal asymptotic p-value not estimable).")
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_result_label, label, NA_character_, status = "ok", duration_time_sec = duration_time_sec)
				return(invisible(NULL))
			}
			if (has_invalid_numeric(result)) {
				msg = paste0("Invalid output detected (NA/NaN/Inf) in ", label)
				message("Skipping ", label, " (non-fatal): ", msg)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_result_label, label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = msg)
				return(invisible(NULL))
			}
			if (is_zero_zero_confidence_interval(label, result)) {
				msg = paste0("Degenerate confidence interval [0, 0] detected in ", label)
				message("Skipping ", label, " (non-fatal): ", msg)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_result_label, label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = msg)
				return(invisible(NULL))
			}
			result = snap_small_numeric_to_zero(result)
			cat("            ", paste(format(result, digits = 3), collapse = " "), "\n")
			duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
			cat(sprintf("              (Duration: %.3gs)\n", duration_time_sec))
			flush.console()
			record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_result_label, label, result, status = "ok", duration_time_sec = duration_time_sec)
			result
		}, error = function(e){
			if (should_record_nonestimable_as_missing(seq_des_inf, label)) {
				msg = nonestimable_error_message(seq_des_inf, label)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_result_label, label, NA_character_, status = "ok", duration_time_sec = duration_time_sec, error_message = msg)
				return(invisible(NULL))
			}
			msg = if (length(e$message) == 0L) "" else e$message
			is_non_fatal = grepl("not implemented", msg, fixed = TRUE) ||
			                 grepl("must implement", msg, fixed = TRUE) ||
			                 grepl("no informative strata are available", msg, fixed = TRUE) ||
			                 grepl("Matching structure is unavailable", msg, fixed = TRUE) ||
			                 grepl("Exact inference is only supported for exact inference classes.", msg, fixed = TRUE) ||
			                 grepl("does not support parametric-bootstrap LR calibration", msg, fixed = TRUE) ||
			                 grepl("Override private\\$supports_lik_ratio_param_bootstrap\\(\\) and simulate_under_lik_null\\(\\)", msg) ||
			                 grepl("singular matrix in 'backsolve'", msg, fixed = TRUE) ||
			                 grepl("G-computation RD: could not compute a finite delta-method standard error.", msg, fixed = TRUE) ||
			                 grepl("G-computation RR: could not compute a finite delta-method standard error.", msg, fixed = TRUE) ||
			                 grepl("G-computation RR: could not compute a finite delta-method confidence interval.", msg, fixed = TRUE) ||
			                 grepl("KK g-computation RR: could not compute a finite delta-method confidence interval.", msg, fixed = TRUE) ||
			                 grepl("G-computation mean difference: could not compute a finite delta-method standard error.", msg, fixed = TRUE) ||
			                 grepl("Zero/one-inflated beta requires y in [0, 1]", msg, fixed = TRUE) ||
			                 grepl("Zhang incidence inference is only supported", msg, fixed = TRUE) ||

					 grepl("This type of inference is only available for incidence", msg, fixed = TRUE) ||
						 grepl("not enough discordant pairs", msg, ignore.case = TRUE) ||
						 grepl("Degenerate confidence interval", msg, fixed = TRUE) ||
						 grepl("inconsistent estimator units", msg, ignore.case = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned NA bounds", msg, fixed = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned non-finite bounds", msg, fixed = TRUE) ||
						 					 						 grepl("Weibull regression failed to converge", msg, fixed = TRUE) ||
				 grepl("Negative binomial regression failed to converge", msg, fixed = TRUE) ||
						 					 						 grepl("Invalid output detected", msg, fixed = TRUE) ||
						 					 						 grepl("missing value where TRUE/FALSE needed", msg, fixed = TRUE) ||
						 					 						 ((grepl("NA/NaN/Inf", msg, fixed = TRUE) || grepl("non-finite standard error", msg, fixed = TRUE) || grepl("could not compute a finite standard error", msg, fixed = TRUE)) &&
						 					 
						 													 (is(seq_des_inf, "InferenceIncidKKCondLogitIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKPoissonCondPoissonIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKStratCoxIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKWeibullFrailtyIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKWilcoxRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKSurvivalRankRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKGEE") ||
						 													  is(seq_des_inf, "InferenceAbstractKKGLMM") ||
						 													  is(seq_des_inf, "InferenceAbstractKKRobustRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceContinKKRobustRegrOneLik") ||
						 													  is(seq_des_inf, "InferenceAbstractKKQuantileRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKQuantileRegrOneLik") ||
						 													  													  is(seq_des_inf, "InferencePropFractionalLogit") ||
						 													  													  is(seq_des_inf, "InferenceIncidRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "InferenceIncidGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "InferenceIncidKKGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidKKGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidKKGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "InferenceIncidKKGCompRiskRatio") ||
						 													  is(seq_des_inf, "InferenceIncidModifiedPoisson") ||
						 													  is(seq_des_inf, "InferencePropGCompMeanDiff") ||
						 													  is(seq_des_inf, "InferencePropGCompMeanDiff") ||
						 													  is(seq_des_inf, "InferencePropZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "InferencePropZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "InferenceContinMultGLS")))
						 													  										if (isTRUE(is_non_fatal)){
				message("Skipping ", label, " (non-fatal): ", e$message)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_result_label, label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = e$message)
			} else {
				stop(e$message)
			}
	})
}

safe_call_family = function(test_family, label, expr){
	if (!should_run_test_family(test_family)) return(invisible(NULL))
	safe_call(label, expr)
}

ensure_estimate_setup_for_filtered_family = function(){
	if (!is_test_family_filter_active() || identical(TEST_FAMILY_FILTER, "estimate")) return(TRUE)
	if (!"compute_estimate" %in% names(seq_des_inf)) return(TRUE)
	isTRUE(tryCatch({
		seq_des_inf$compute_estimate()
		TRUE
	}, error = function(e){
		msg = if (length(e$message) == 0L) "" else e$message
		message(
			"          Skipping ",
			TEST_FAMILY_FILTER,
			" family for ",
			inference_result_label,
			" because compute_estimate setup failed: ",
			msg
		)
		FALSE
	}))
}

call_direct_asymp = function(method_name, testing_type, ...){
	if (!should_run_test_family(testing_type)) return(invisible(NULL))
	if (!method_name %in% names(seq_des_inf)) return(invisible(NULL))
	if (!supports_direct_testing_type(testing_type)) {
		message("          Skipping ", method_name, " (not implemented for testing_type = ", testing_type, ")")
		return(invisible(NULL))
	}
	method_fn = seq_des_inf[[method_name]]
	if (!is.function(method_fn)) return(invisible(NULL))
	result = tryCatch(
		do.call(method_fn, list(...)),
		error = function(e){
			msg = if (length(e$message) == 0L) "" else e$message
			if (grepl("Must be implemented by concrete class or shared helper.", msg, fixed = TRUE) ||
			    grepl("does not expose a likelihood-test specification", msg, fixed = TRUE) ||
			    grepl("not implement", msg, ignore.case = TRUE) ||
			    grepl("must implement", msg, ignore.case = TRUE)) {
				message("          Skipping ", method_name, " (not implemented)")
				return(structure(list(skipped = TRUE), class = "edi_skip_direct"))
			}
			stop(e)
		}
	)
	if (inherits(result, "edi_skip_direct")) return(invisible(NULL))
	safe_call(method_name, result)
}

	if (!ensure_estimate_setup_for_filtered_family()) return(invisible(NULL))

	if (is(seq_des_inf, "InferenceOrdinalJonckheereTerpstraTest")){
		safe_call_family("exact", "compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		safe_call_family("estimate", "compute_estimate", seq_des_inf$compute_estimate())
		return(invisible(NULL))
	}

	supports_exact_inference = is(seq_des_inf, "InferenceExact") || is(seq_des_inf, "InferenceIncidenceExactZhang") || is(seq_des_inf, "InferenceIncidExactZhang")
	if (should_run_test_family("exact") && response_type == "incidence" && supports_exact_inference){
		safe_call_family("exact", "compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		safe_call_family("exact", "compute_exact_confidence_interval", seq_des_inf$compute_exact_confidence_interval(args_for_type = list(Zhang = list(combination_method = "Fisher", pval_epsilon = pval_epsilon))))
	}
	skip_asymp = is(seq_des_inf, "InferenceExact") &&
		!is(seq_des_inf, "InferenceAsymp") &&
		!is(seq_des_inf, "InferenceAsympLik")

	safe_call_family("estimate", "compute_estimate", seq_des_inf$compute_estimate())
	if (should_run_test_family("asymp") && !skip_asymp && !skip_mle_pval){
		if ("compute_asymp_log_rank_two_sided_pval_for_treatment_effect" %in% names(seq_des_inf)) {
			safe_call_family("asymp", "compute_asymp_log_rank_two_sided_pval_for_treatment_effect", seq_des_inf$compute_asymp_log_rank_two_sided_pval_for_treatment_effect())
		}
		if ("compute_asymp_two_sided_pval" %in% names(seq_des_inf)) {
			safe_call_family("asymp", "compute_asymp_two_sided_pval", seq_des_inf$compute_asymp_two_sided_pval())
		}
	}
	if (should_run_test_family("asymp") && !skip_asymp && !skip_ci){
		safe_call_family("asymp", "compute_asymp_confidence_interval", seq_des_inf$compute_asymp_confidence_interval(0.05))
	}
	if (!skip_asymp && !skip_mle_pval){
			call_direct_asymp("compute_wald_two_sided_pval", "wald")
			call_direct_asymp("compute_score_two_sided_pval", "score")
			call_direct_asymp("compute_lik_ratio_two_sided_pval", "lik_ratio")
			call_direct_asymp("compute_gradient_two_sided_pval", "gradient")
		}
	if (!skip_asymp && !skip_ci){
		call_direct_asymp("compute_wald_confidence_interval", "wald", 0.05)
		call_direct_asymp("compute_score_confidence_interval", "score", 0.05)
		call_direct_asymp("compute_lik_ratio_confidence_interval", "lik_ratio", 0.05)
		call_direct_asymp("compute_gradient_confidence_interval", "gradient", 0.05)
	}
	safe_call_debug = function(label, expr) {
		if (is_row_completed(rep_curr, beta_T, dataset_name, response_type, design_type, inference_result_label, label)) {
			return(invisible(NULL))
		}
		if (!is.null(pending_rep_header)) { message(pending_rep_header); pending_rep_header <<- NULL }
		if (!is.null(pending_beta_header)) { message(pending_beta_header); pending_beta_header <<- NULL }
		if (!is.null(pending_dataset_header)) { message(pending_dataset_header); pending_dataset_header <<- NULL }
		if (!is.null(pending_response_header)) { message(pending_response_header); pending_response_header <<- NULL }
		if (!is.null(pending_design_header)) { message(pending_design_header); pending_design_header <<- NULL }
		if (!is.null(pending_banner)) { message(pending_banner); pending_banner <<- NULL }
		message("          Calling ", label, "()")
		start_elapsed = unname(proc.time()[["elapsed"]])
		debug_result = tryCatch(expr, error = function(e) {
			dur = unname(proc.time()[["elapsed"]]) - start_elapsed
			record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type,
						  inference_result_label, label, NA_character_, status = "error",
						  duration_time_sec = dur, error_message = e$message)
			NULL
		})
		if (is.null(debug_result)) return(invisible(NULL))
		duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
		if (identical(label, "approximate_bootstrap_distribution_beta_hat_T_debug")) {
			stats_vec = c(
				debug_result$prop_illegal_values,
				debug_result$prop_iterations_with_errors,
				debug_result$prop_iterations_with_warnings
			)
			cat(sprintf("            prop_illegal=%.3f  prop_err=%.3f  prop_warn=%.3f\n",
					stats_vec[1], stats_vec[2], stats_vec[3]))
		} else {
			stats_vec = c(
				debug_result$prop_iterations_with_errors,
				debug_result$prop_iterations_with_warnings,
				debug_result$prop_illegal_values
			)
			cat(sprintf("            prop_err=%.3f  prop_warn=%.3f  prop_illegal=%.3f\n",
					stats_vec[1], stats_vec[2], stats_vec[3]))
		}
		cat(sprintf("              (Duration: %.3gs)\n", duration_time_sec))
		record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type,
					  inference_result_label, label, stats_vec, status = "ok",
					  duration_time_sec = duration_time_sec)
	}

	if (should_run_test_family("bootstrap") && run_debug_resampling && !skip_slow && !skip_bootstrap){
		safe_call_debug("approximate_bootstrap_distribution_beta_hat_T_debug",
						seq_des_inf$approximate_bootstrap_distribution_beta_hat_T(B = B_debug, debug = TRUE))
	}
	if (should_run_test_family("bayesian_bootstrap") && run_debug_resampling && !skip_slow && !skip_bayesian_bootstrap){
		safe_call_debug("approximate_bayesian_bootstrap_distribution_beta_hat_T_debug",
						seq_des_inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(B = B_debug, debug = TRUE, show_progress = FALSE))
	}
	# Nonparametric bootstrap CI — default type first (warms the distribution cache), then extra types reuse it
	if (should_run_test_family("bootstrap") && !skip_slow && !skip_ci && !skip_bootstrap){
		safe_call("compute_bootstrap_confidence_interval", seq_des_inf$compute_bootstrap_confidence_interval(B = r, na.rm = TRUE, show_progress = FALSE))
		for (boot_ci_type in c("basic", "bca", "studentized")) {
			safe_call(paste0("compute_bootstrap_confidence_interval_", boot_ci_type),
					  seq_des_inf$compute_bootstrap_confidence_interval(B = r, type = boot_ci_type, na.rm = TRUE, show_progress = FALSE))
		}
	}
	# Bayesian bootstrap CI — default type first (warms the distribution cache), then extra types reuse it
	if (should_run_test_family("bayesian_bootstrap") && !skip_slow && !skip_ci && !skip_bayesian_bootstrap){
		safe_call("compute_bayesian_bootstrap_confidence_interval", seq_des_inf$compute_bayesian_bootstrap_confidence_interval(B = r, na.rm = TRUE, show_progress = FALSE))
		for (bayes_ci_type in c("basic", "wald", "bca", "studentized")) {
			safe_call(paste0("compute_bayesian_bootstrap_confidence_interval_", bayes_ci_type),
					  seq_des_inf$compute_bayesian_bootstrap_confidence_interval(B = r, type = bayes_ci_type, na.rm = TRUE, show_progress = FALSE))
		}
	}
	# Nonparametric bootstrap p-val — default first, extra types reuse distribution cache
	if (should_run_test_family("bootstrap") && !skip_slow && !skip_bootstrap){
		safe_call("compute_bootstrap_two_sided_pval", seq_des_inf$compute_bootstrap_two_sided_pval(B = r, na.rm = TRUE, show_progress = FALSE))
		for (boot_pval_type in c("symmetric", "bca", "studentized")) {
			safe_call(paste0("compute_bootstrap_two_sided_pval_", boot_pval_type),
					  seq_des_inf$compute_bootstrap_two_sided_pval(B = r, type = boot_pval_type, na.rm = TRUE, show_progress = FALSE))
		}
	}
	# Bayesian bootstrap p-val — default first, extra types reuse distribution cache
	if (should_run_test_family("bayesian_bootstrap") && !skip_slow && !skip_bayesian_bootstrap){
		safe_call("compute_bayesian_bootstrap_two_sided_pval", seq_des_inf$compute_bayesian_bootstrap_two_sided_pval(B = r, na.rm = TRUE, show_progress = FALSE))
		for (bayes_pval_type in c("symmetric", "wald", "bca", "studentized")) {
			safe_call(paste0("compute_bayesian_bootstrap_two_sided_pval_", bayes_pval_type),
					  seq_des_inf$compute_bayesian_bootstrap_two_sided_pval(B = r, type = bayes_pval_type, na.rm = TRUE, show_progress = FALSE))
		}
	}
	if (should_run_test_family("parametric_bootstrap") && !skip_slow && supports_parametric_bootstrap){
		safe_call("compute_lik_ratio_bootstrap_two_sided_pval", seq_des_inf$compute_lik_ratio_bootstrap_two_sided_pval(B = r, show_progress = FALSE))
	}
	if (should_run_test_family("parametric_bootstrap") && !skip_slow && !skip_ci && supports_parametric_bootstrap_ci){
		safe_call("compute_lik_ratio_bootstrap_confidence_interval",
			seq_des_inf$compute_lik_ratio_bootstrap_confidence_interval(
				B = r,
				show_progress = FALSE,
				max_root_iterations = param_boot_ci_max_root_iterations
			)
		)
	}
	if (should_run_test_family("jackknife") && !skip_slow && supports_jackknife){
		if (response_type != "count") {
			safe_call("compute_jackknife_estimate", seq_des_inf$compute_jackknife_estimate())
		}
		safe_call("compute_jackknife_wald_two_sided_pval", seq_des_inf$compute_jackknife_wald_two_sided_pval())
	}
	if (should_run_test_family("jackknife") && !skip_slow && supports_jackknife){
		safe_call("compute_jackknife_wald_confidence_interval", seq_des_inf$compute_jackknife_wald_confidence_interval())
	}
	if (should_run_test_family("rand") && !skip_slow && !skip_rand && !skip_rand_pval && response_type %in% c("continuous", "survival", "proportion")){
		if (run_debug_resampling) {
			safe_call_debug("approximate_randomization_distribution_beta_hat_T_debug",
						seq_des_inf$approximate_randomization_distribution_beta_hat_T(r = r_debug, debug = TRUE))
		}
		safe_call("compute_rand_two_sided_pval", seq_des_inf$compute_rand_two_sided_pval(r = r, show_progress = FALSE))
		transform_for_rand = switch(
			response_type,
			continuous = "none",
			proportion = "logit",
			count = "log",
			survival = "log",
			"none"
		)
		delta_for_rand = 0.5
		safe_call("compute_rand_two_sided_pval(delta=0.5)",
				seq_des_inf$compute_rand_two_sided_pval(r = r, delta = delta_for_rand, transform_responses = transform_for_rand, show_progress = FALSE))
	}

	if (should_run_test_family("rand") && !skip_slow && !skip_rand && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous", "proportion", "count")){
		safe_call("compute_rand_confidence_interval", seq_des_inf$compute_rand_confidence_interval(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE))
	}
	if (should_run_test_family("rand_custom") && response_type != "incidence"){
		seq_des_inf$set_custom_randomization_statistic_cpp(welch_t_stat_cpp)
		if (!skip_slow && !skip_rand_pval){
			safe_call("compute_rand_two_sided_pval(custom)", seq_des_inf$compute_rand_two_sided_pval(r = r, show_progress = FALSE))
		}
		if (!skip_slow && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous")){
			if (!skip_ci_rand_custom){
				safe_call("compute_rand_confidence_interval(custom)", seq_des_inf$compute_rand_confidence_interval(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE))
			} else {
				message("    Skipping compute_rand_confidence_interval(custom) (too slow)")
			}
		}
		seq_des_inf$set_custom_randomization_statistic_cpp(NULL)
	}
}

instantiate_inference_generator = function(class_gen, des_obj, model_formula = NULL){
	init_fn = class_gen$public_methods$initialize
	init_formals = if (is.function(init_fn)) names(formals(init_fn)) else character()
	candidate_arg_sets = list(
		list(des_obj = des_obj, model_formula = model_formula, verbose = FALSE),
		list(des_obj = des_obj, model_formula = model_formula),
		list(des_obj = des_obj, verbose = FALSE),
		list(des_obj = des_obj)
	)
	for (args in candidate_arg_sets){
		args = args[names(args) %in% init_formals]
		obj = tryCatch(do.call(class_gen$new, args), error = function(e) NULL)
		if (!is.null(obj)) return(obj)
	}
	NULL
}

run_exhaustive_remaining_inference_classes = function(des_obj, response_type, design_type, dataset_name, n_rows, n_cols, model_formula = NULL){
	ns = asNamespace("EDI")
	nms = ls(ns, all.names = TRUE)
	gen_names = nms[vapply(nms, function(nm) inherits(get(nm, envir = ns), "R6ClassGenerator"), logical(1))]
	gen_names = gen_names[grepl("^Inference", gen_names)]
	gen_names = gen_names[!is_skipped_inference_label(gen_names)]
	gen_names = gen_names[!grepl("Abstract|Suite|Custom|RandCI$|NoParamBootstrap$", gen_names)]
	gen_names = setdiff(gen_names, c(
		"Inference",
		"InferenceAsymp",
		"InferenceBoot",
		"InferenceRand",
		"InferenceExact",
		"InferenceKKPassThrough",
		"InferenceKKPassThroughCompound",
		"InferenceMLEorKMforGLMs",
		"InferenceMLEorKMSummaryTable",
		"InferenceBaiAdjustedT",
		"InferenceCustomAsymp",
		"InferenceCustomBoot",
		"InferenceCustomRand",
		"InferenceIncidKKCondLogitPlusGLMMOneLik",
		"InferenceRandCI",
		"InferenceQuantileRandCI"
	))
	is_response_family_mismatch = function(class_name, response_type){
		switch(
			response_type,
			continuous = grepl("Count|Incid|Survival|Ordinal|^InferenceProp", class_name),
			incidence = grepl("Count|Contin|Survival|Ordinal|^InferenceProp", class_name),
			proportion = grepl("Count|Contin|Incid|Survival|Ordinal", class_name),
			count = grepl("Contin|Incid|Survival|Ordinal|^InferenceProp", class_name),
			survival = grepl("Count|Contin|Incid|Ordinal|^InferenceProp", class_name),
			ordinal = grepl("Count|Contin|Incid|Survival|^InferenceProp", class_name),
			FALSE
		)
	}
	is_kk_design_type = design_type %in% c("KK21", "KK21stepwise", "KK14", "FixedBinaryMatch")
	for (class_name in sort(unique(gen_names))){
		if (is_response_family_mismatch(class_name, response_type)) next
		if (!is_kk_design_type && grepl("KK", class_name, fixed = TRUE)) next
		class_gen = get(class_name, envir = ns)
		inf_obj = instantiate_inference_generator(class_gen, des_obj, model_formula = model_formula)
		if (is.null(inf_obj)) next
		class_label = class(inf_obj)[1L]
		if (!startsWith(class_label, "Inference")) next
		inference_banner(class_label)
		tryCatch({
			run_inference_checks(inf_obj, response_type, design_type, dataset_name, n_rows, n_cols, exhaustive_sweep = TRUE)
		}, error = function(e){
			msg = if (length(e$message) == 0L) "" else e$message
			message("  Skipping ", class_label, " (exhaustive sweep): ", msg)
			return(invisible(NULL))
		})
	}
	invisible(NULL)
}

run_tests_for_response = function(response_type, design_type, dataset_name, model_formula = NULL){
	.comprehensive_current_model_formula <<- model_formula
	on.exit({
		if (exists(".comprehensive_current_model_formula", envir = .GlobalEnv, inherits = FALSE)) {
			rm(".comprehensive_current_model_formula", envir = .GlobalEnv)
		}
	}, add = TRUE)
	apply_treatment_effect_and_noise = function(y_t, w_t, response_type){
		eps = rnorm(1, 0, SD_NOISE)
		bt = ifelse(w_t == 1, beta_T, 0)
		if (response_type == "continuous") return(y_t + bt + eps)
		if (response_type == "incidence") {
			p_base = if (is.finite(y_t) && y_t >= 0 && y_t <= 1) y_t else stats::plogis(y_t)
			p_base = pmin(0.95, pmax(0.05, p_base))
			p_t = plogis(qlogis(p_base) + bt + eps)
			return(as.numeric(stats::rbinom(1, size = 1, prob = p_t)))
		}
		if (response_type == "proportion"){
			return(pmin(1, pmax(0, y_t + bt + eps)))
		}
		if (response_type == "count"){
			lambda_t = pmax(.Machine$double.eps, y_t * exp(bt + eps))
			return(as.numeric(stats::rpois(1, lambda = lambda_t)))
		}
		if (response_type == "survival") return(pmax(.Machine$double.eps, y_t * exp(bt + eps)))
		if (response_type == "ordinal"){
			# For ordinal, we'll just use a simple shift and re-cut if needed, 
			# but here y_t is already the ordinal level.
			# Let's just do a simple shift and round.
			# We use a larger noise multiplier (5x) to ensure some subjects jump categories,
			# otherwise with SD_NOISE=0.1, matched pairs are almost always concordant.
			return(as.integer(max(1, round(y_t + bt + 5 * eps))))
		}
		stop("Unsupported response_type: ", response_type)
	}

	D = datasets_and_response_models[[dataset_name]]
	n_X = nrow(D$X)
	p_X = ncol(D$X)
	X_design = as.data.frame(D$X)
	if (identical(design_type, "KK21stepwise") && ncol(X_design) > 20L){
		message(
			"    Truncating KK21stepwise test covariates from ",
			ncol(X_design),
			" to 20 to keep stepwise-weight tests bounded in runtime."
		)
		X_design = X_design[, seq_len(20L), drop = FALSE]
	}

	if (nrow(X_design) %% 4L != 0L){
		n_keep = nrow(X_design) - (nrow(X_design) %% 4L)
		message(
			"    Truncating test rows from ",
			nrow(X_design),
			" to ",
			n_keep,
			" so the design matrix row count is divisible by 4."
		)
		X_design = X_design[seq_len(n_keep), , drop = FALSE]
	}

	n = nrow(X_design)
	dataset_n_rows = nrow(X_design)
	dataset_n_cols = ncol(X_design)
	dead = rep(1, n)
	y = D$y_original[[response_type]]
	t_f = quantile(y, .95)

	if (response_type == "survival"){
		for (i in 1 : n){
			if (runif(1) < prob_censoring || y[i] >= t_f){
				y[i] = runif(1, 0, y[i])
				dead[i] = 0
			}
		}
	}

	cluster_design_setup = NULL
	if (identical(design_type, "FixedCluster")) {
		cluster_design_setup = add_assignment_only_cluster_id(X_design)
	} else if (identical(design_type, "FixedBlockedCluster")) {
		cluster_design_setup = add_assignment_only_cluster_id(X_design, strata_cols = names(X_design)[2:min(2, ncol(X_design))])
	}
	X_design_for_design = if (is.null(cluster_design_setup)) X_design else cluster_design_setup$X

	# For sequential designs that use stratification, we MUST discretize continuous strata
	# because the DesignSeqOneByOne base class explicitly disallows numeric strata.
	X_design_sequential_strata = X_design
	strata_cols_to_use = names(X_design)[1:min(2, ncol(X_design))]
	if (design_type %in% c("SPBR", "PocockSimon", "RandomBlockSize")) {
		for (col in strata_cols_to_use) {
			if (is.numeric(X_design_sequential_strata[[col]])) {
				med = stats::median(X_design_sequential_strata[[col]], na.rm = TRUE)
				val_vec = ifelse(X_design_sequential_strata[[col]] <= med, "low", "high")
				if (length(unique(val_vec)) < 2) {
					val_vec = ifelse(X_design_sequential_strata[[col]] < max(X_design_sequential_strata[[col]], na.rm = TRUE), "low", "high")
				}
				X_design_sequential_strata[[col]] = factor(val_vec)
			}
		}
	}

	# If model_formula is provided, ensure it's passed to the design if it's NOT NULL.
	# The base Design class initialize handles model_formula.
	design_formula = if (is.null(model_formula)) ~ . else model_formula

	des_obj = tryCatch(switch(design_type,
		KK21 =         DesignSeqOneByOneKK21$new(        response_type = response_type, n = n, model_formula = design_formula),
		KK21stepwise = DesignSeqOneByOneKK21stepwise$new(response_type = response_type, n = n, model_formula = design_formula),
		KK14 =         DesignSeqOneByOneKK14$new(        response_type = response_type, n = n, model_formula = design_formula),
		Bernoulli =    DesignSeqOneByOneBernoulli$new(   response_type = response_type, n = n, model_formula = design_formula),
		Efron =        DesignSeqOneByOneEfron$new(       response_type = response_type, n = n, model_formula = design_formula),
		Atkinson =     DesignSeqOneByOneAtkinson$new(    response_type = response_type, n = n, model_formula = design_formula),
		iBCRD =        DesignSeqOneByOneiBCRD$new(       response_type = response_type, n = n, model_formula = design_formula),
		Urn =          DesignSeqOneByOneUrn$new(         response_type = response_type, n = n, model_formula = design_formula),
		RandomBlockSize = DesignSeqOneByOneRandomBlockSize$new( strata_cols = strata_cols_to_use, response_type = response_type, n = n, model_formula = design_formula),
		SPBR =         DesignSeqOneByOneSPBR$new(        strata_cols = strata_cols_to_use, block_size = 4, response_type = response_type, n = n, model_formula = design_formula),
		PocockSimon =  DesignSeqOneByOnePocockSimon$new( strata_cols = strata_cols_to_use, response_type = response_type, n = n, model_formula = design_formula),
		FixedBernoulli = DesignFixedBernoulli$new( response_type = response_type, n = n, model_formula = design_formula),
		FixediBCRD =     DesignFixediBCRD$new(     response_type = response_type, n = n, model_formula = design_formula),
		FixedBlocking =  DesignFixedBlocking$new(  strata_cols = strata_cols_to_use, response_type = response_type, n = n, model_formula = design_formula),
		FixedCluster =   DesignFixedCluster$new(   cluster_col = cluster_design_setup$cluster_col, response_type = response_type, n = n, model_formula = design_formula),
		FixedBlockedCluster = DesignFixedBlockedCluster$new( strata_cols = names(X_design)[2:min(2, ncol(X_design))], cluster_col = cluster_design_setup$cluster_col, response_type = response_type, n = n, model_formula = design_formula),
		FixedBinaryMatch = DesignFixedBinaryMatch$new( response_type = response_type, n = n, model_formula = design_formula),
		FixedGreedy =    DesignFixedGreedy$new(    response_type = response_type, n = n, model_formula = design_formula),
		FixedRerandomization = DesignFixedRerandomization$new( response_type = response_type, n = n, model_formula = design_formula),
		FixedMatchingGreedy = DesignFixedMatchingGreedyPairSwitching$new( response_type = response_type, n = n, model_formula = design_formula),
		FixedDOptimal =  DesignFixedDOptimal$new(  response_type = response_type, n = n, model_formula = design_formula),
		FixedAOptimal =  DesignFixedAOptimal$new(  response_type = response_type, n = n, model_formula = design_formula),
		stop("Unsupported design_type: ", design_type)
	), error = function(e){ message("    Skipping design (creation error): ", e$message); NULL })
	if (is.null(des_obj)) return(invisible(NULL))

	if (inherits(des_obj, "DesignSeqOneByOne")){
		seq_ok = tryCatch({
			for (t in 1 : n){
				w_t = des_obj$add_one_subject_to_experiment_and_assign(X_design_sequential_strata[t, , drop = FALSE])
				y_t = apply_treatment_effect_and_noise(y[t], w_t, response_type)
				des_obj$add_one_subject_response(t, y_t, dead[t])
			}
			TRUE
		}, error = function(e){ message("    Skipping design (seq error): ", e$message); FALSE })
		if (!seq_ok) return(invisible(NULL))
	} else {
		# It is a DesignFixed but not a DesignSeqOneByOne
		setup_ok = tryCatch({
			des_obj$add_all_subjects_to_experiment(X_design_for_design)
			des_obj$assign_w_to_all_subjects()
			TRUE
		}, error = function(e){ message("    Skipping design: ", e$message); FALSE })
		if (!setup_ok) return(invisible(NULL))
		w = des_obj$get_w()
		for (t in 1 : n){
			y_t = apply_treatment_effect_and_noise(y[t], w[t], response_type)
			des_obj$add_one_subject_response(t, y_t, dead[t])
		}
	}

	is_kk_design = design_type %in% c("KK21", "KK21stepwise", "KK14", "FixedBinaryMatch")
	if (response_type == "continuous"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceAllSimpleMeanDiffPooledVar")
		run_inference_checks(InferenceAllSimpleMeanDiffPooledVar$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		if (is_kk_design){
			if (design_type == "KK14"){
				inference_banner("InferenceBaiAdjustedTKK14")
				run_inference_checks(InferenceBaiAdjustedTKK14$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			}
			if (design_type %in% c("KK21", "KK21stepwise")){
				inference_banner("InferenceBaiAdjustedTKK21")
				run_inference_checks(InferenceBaiAdjustedTKK21$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			}
			run_inference_checks_for_class(InferenceAllKKMeanDiffIVWC, "InferenceAllKKMeanDiffIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_for_class(InferenceAllKKWilcoxIVWC, "InferenceAllKKWilcoxIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinKKRobustRegrIVWC, "InferenceContinKKRobustRegrIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinKKRobustRegrOneLik, "InferenceContinKKRobustRegrOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinKKOLSOneLik, "InferenceContinKKOLSOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinKKGLMM, "InferenceContinKKGLMM", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinKKQuantileRegrIVWC, "InferenceContinKKQuantileRegrIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinKKQuantileRegrOneLik, "InferenceContinKKQuantileRegrOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		} else {
			run_inference_checks_both_paths(InferenceContinRobustRegr, "InferenceContinRobustRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinQuantileRegr, "InferenceContinQuantileRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinLin, "InferenceContinLin", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceContinOLS, "InferenceContinOLS", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		}
	}

	if (response_type == "incidence"){
		supports_exact_fisher_design =
			design_type %in% c("iBCRD", "FixediBCRD", "FixedBlocking", "RandomBlockSize", "SPBR") ||
			is_kk_design
		supports_exact_binomial_design =
			design_type %in% c("KK14", "FixedBinaryMatch")

		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		if (supports_exact_fisher_design) {
			inference_banner("InferenceIncidExactFisher")
			run_inference_checks(InferenceIncidExactFisher$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		}
		if (supports_exact_binomial_design) {
			inference_banner("InferenceIncidExactBinomial")
			run_inference_checks(InferenceIncidExactBinomial$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		}
		inference_banner("InferenceIncidWald")
		run_inference_checks(InferenceIncidWald$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceIncidCMH")
		err_msg_cmh = tryCatch({
			run_inference_checks(InferenceIncidCMH$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			NULL
		}, error = function(e) if (length(e$message) == 0L) "" else e$message)
		if (!is.null(err_msg_cmh)) {
			if (grepl("equal block sizes", err_msg_cmh, fixed = TRUE) || grepl("even allocation", err_msg_cmh, fixed = TRUE)) {
				message("  Skipping InferenceIncidCMH: ", err_msg_cmh)
			} else {
				stop(err_msg_cmh)
			}
		}
		inference_banner("InferenceIncidExtendedRobins")
		err_msg_er = tryCatch({
			run_inference_checks(InferenceIncidExtendedRobins$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			NULL
		}, error = function(e) if (length(e$message) == 0L) "" else e$message)
		if (!is.null(err_msg_er)) {
			if (grepl("equal block sizes", err_msg_er, fixed = TRUE) || grepl("even allocation", err_msg_er, fixed = TRUE)) {
				message("  Skipping InferenceIncidExtendedRobins: ", err_msg_er)
			} else {
				stop(err_msg_er)
			}
		}
		if (is_kk_design){
			run_inference_checks_for_class(InferenceAllKKMeanDiffIVWC, "InferenceAllKKMeanDiffIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidKKCondLogitOneLik, "InferenceIncidKKCondLogitOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidKKCondLogitPlusGLMMOneLik, "InferenceIncidKKCondLogitPlusGLMMOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidKKCondLogitIVWC, "InferenceIncidKKCondLogitIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidKKGEE, "InferenceIncidKKGEE", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			inference_banner("InferenceIncidKKNewcombeRiskDiff")
			run_inference_checks(InferenceIncidKKNewcombeRiskDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			run_inference_checks_both_paths(InferenceIncidKKGCompRiskDiff, "InferenceIncidKKGCompRiskDiff", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidKKGCompRiskRatio, "InferenceIncidKKGCompRiskRatio", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidKKModifiedPoisson, "InferenceIncidKKModifiedPoisson", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidKKCondLogitPlusGLMMOneLik, "InferenceIncidKKCondLogitPlusGLMMOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		}
		run_inference_checks_both_paths(InferenceIncidLogRegr, "InferenceIncidLogRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		if (!is_kk_design){
			run_inference_checks_both_paths(InferenceIncidProbitRegr, "InferenceIncidProbitRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			inference_banner("InferenceIncidMiettinenNurminenRiskDiff")
			run_inference_checks(InferenceIncidMiettinenNurminenRiskDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			inference_banner("InferenceIncidNewcombeRiskDiff")
			run_inference_checks(InferenceIncidNewcombeRiskDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			run_inference_checks_both_paths(InferenceIncidRiskDiff, "InferenceIncidRiskDiff", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidGCompRiskDiff, "InferenceIncidGCompRiskDiff", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidGCompRiskRatio, "InferenceIncidGCompRiskRatio", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidModifiedPoisson, "InferenceIncidModifiedPoisson", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidLogBinomial, "InferenceIncidLogBinomial", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceIncidBinomialIdentityRiskDiff, "InferenceIncidBinomialIdentityRiskDiff", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		}
	}

	if (response_type == "proportion"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		if (is_kk_design){
			run_inference_checks_for_class(InferenceAllKKMeanDiffIVWC, "InferenceAllKKMeanDiffIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_for_class(InferenceAllKKWilcoxIVWC, "InferenceAllKKWilcoxIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferencePropKKGEE, "InferencePropKKGEE", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferencePropKKGLMM, "InferencePropKKGLMM", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferencePropKKQuantileRegrIVWC, "InferencePropKKQuantileRegrIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferencePropKKQuantileRegrOneLik, "InferencePropKKQuantileRegrOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		}
		run_inference_checks_both_paths(InferencePropBetaRegr, "InferencePropBetaRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		if (!is_kk_design){
			run_inference_checks_both_paths(InferencePropZeroOneInflatedBetaRegr, "InferencePropZeroOneInflatedBetaRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferencePropGCompMeanDiff, "InferencePropGCompMeanDiff", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferencePropFractionalLogit, "InferencePropFractionalLogit", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		}
	}

	if (response_type == "count"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		if (is_kk_design){
			run_inference_checks_for_class(InferenceAllKKMeanDiffIVWC, "InferenceAllKKMeanDiffIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_for_class(InferenceAllKKWilcoxIVWC, "InferenceAllKKWilcoxIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountPoissonKKGEE, "InferenceCountPoissonKKGEE", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountKKGLMM, "InferenceCountKKGLMM", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountKKHurdlePoissonOneLik, "InferenceCountKKHurdlePoissonOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountKKCondPoissonOneLik, "InferenceCountKKCondPoissonOneLik", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		}
		if (!is_kk_design){
			run_inference_checks_both_paths(InferenceCountPoisson, "InferenceCountPoisson", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountRobustPoisson, "InferenceCountRobustPoisson", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountQuasiPoisson, "InferenceCountQuasiPoisson", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountZeroInflatedPoisson, "InferenceCountZeroInflatedPoisson", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountZeroInflatedNegBin, "InferenceCountZeroInflatedNegBin", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountHurdlePoisson, "InferenceCountHurdlePoisson", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceCountHurdleNegBin, "InferenceCountHurdleNegBin", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		}
		run_inference_checks_both_paths(InferenceCountNegBin, "InferenceCountNegBin", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
	}

	if (response_type == "survival"){
		if (is_kk_design){
			is_censoring_skip_error = function(msg){
				grepl("only available for uncensored", msg, fixed = TRUE) ||
				grepl("does not currently support censored", msg, fixed = TRUE) ||
				grepl("does not support censored", msg, fixed = TRUE)
			}
			for (kk_surv_class in Filter(function(class_gen) should_run_inference_label(class_gen$classname), list(InferenceAllKKMeanDiffIVWC, InferenceAllKKWilcoxIVWC))){
				class_name = kk_surv_class$classname
				inference_banner(class_name)
				err_msg = tryCatch({
					run_inference_checks(kk_surv_class$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
					NULL
				}, error = function(e) if (length(e$message) == 0L) "" else e$message)
				if (!is.null(err_msg)){
					if (is_censoring_skip_error(err_msg)) message("  Skipping ", class_name, " (censored data): ", err_msg)
					else stop(err_msg)
				}
			}
			for (kk_surv_class in Filter(function(class_gen) should_run_inference_label(class_gen$classname), list(
				InferenceSurvivalKKClaytonCopulaIVWC,
				InferenceSurvivalKKLWACoxPHIVWC,
				InferenceSurvivalKKStratCoxPHIVWC,
				InferenceSurvivalKKRankRegrIVWC,
				InferenceSurvivalKKWeibullFrailtyIVWC,
				InferenceSurvivalKKClaytonCopulaOneLik,
				InferenceSurvivalKKLWACoxPHOneLik,
				InferenceSurvivalKKStratCoxPHOneLik,
				InferenceSurvivalKKWeibullFrailtyOneLik
			))){
				class_name = kk_surv_class$classname
				inference_banner(class_name)
				err_msg = tryCatch({
					run_inference_checks(kk_surv_class$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
					NULL
				}, error = function(e) if (length(e$message) == 0L) "" else e$message)
				if (!is.null(err_msg)){
					if (grepl("only available for uncensored", err_msg, fixed = TRUE)) message("  Skipping ", class_name, " (censored data): ", err_msg)
					else stop(err_msg)
				}
			}
		}
		inference_banner("InferenceAllSimpleWilcox")
		err_msg_sw = tryCatch({
			run_inference_checks(InferenceAllSimpleWilcox$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			NULL
		}, error = function(e) if (length(e$message) == 0L) "" else e$message)
		if (!is.null(err_msg_sw)){
			if (grepl("does not support censored", err_msg_sw, fixed = TRUE)) {
				message("  Skipping InferenceAllSimpleWilcox (censored data): ", err_msg_sw)
				sw_label = paste0("InferenceAllSimpleWilcox [design_formula=", paste(deparse(model_formula), collapse = " "), "]")
				record_existing_error_keys_as_skipped(
					sw_label, n_X, p_X,
					paste0("Skipped unsupported censored survival inference: ", err_msg_sw)
				)
			}
			else stop(err_msg_sw)
		}
		inference_banner("InferenceSurvivalGehanWilcox")
		run_inference_checks(InferenceSurvivalGehanWilcox$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceSurvivalLogRank")
		run_inference_checks(InferenceSurvivalLogRank$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceSurvivalRestrictedMeanDiff")
		run_inference_checks(InferenceSurvivalRestrictedMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceSurvivalKMDiff")
		run_inference_checks(InferenceSurvivalKMDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		run_inference_checks_both_paths(InferenceSurvivalWeibullRegr, "InferenceSurvivalWeibullRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		run_inference_checks_both_paths(InferenceSurvivalDepCensTransformRegr, "InferenceSurvivalDepCensTransformRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		run_inference_checks_both_paths(InferenceSurvivalCoxPHRegr, "InferenceSurvivalCoxPHRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
		run_inference_checks_both_paths(InferenceSurvivalStratCoxPHRegr, "InferenceSurvivalStratCoxPHRegr", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
	}

	if (response_type == "ordinal"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		if (is_kk_design){
			run_inference_checks_for_class(InferenceAllKKMeanDiffIVWC, "InferenceAllKKMeanDiffIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_for_class(InferenceAllKKWilcoxIVWC, "InferenceAllKKWilcoxIVWC", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			inference_banner("InferenceOrdinalKKGEE")
			run_inference_checks(InferenceOrdinalKKGEE$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			run_inference_checks_both_paths(InferenceOrdinalKKGLMM, "InferenceOrdinalKKGLMM", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceOrdinalKKCLMM, "InferenceOrdinalKKCLMM", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceOrdinalKKCLMMProbit, "InferenceOrdinalKKCLMMProbit", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceOrdinalKKCLMMCauchit, "InferenceOrdinalKKCLMMCauchit", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			run_inference_checks_both_paths(InferenceOrdinalKKCLMMCloglog, "InferenceOrdinalKKCLMMCloglog", des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
			inference_banner("InferenceOrdinalKKCondAdjCatLogitRegr")
			run_inference_checks(InferenceOrdinalKKCondAdjCatLogitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			inference_banner("InferenceOrdinalKKCondAdjCatLogitRegr")
			run_inference_checks(InferenceOrdinalKKCondAdjCatLogitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
			inference_banner("InferenceOrdinalPairedSignTest")
			run_inference_checks(InferenceOrdinalPairedSignTest$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		}
		inference_banner("InferenceOrdinalAdjCatLogitRegr")
		run_inference_checks(InferenceOrdinalAdjCatLogitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalAdjCatLogitRegr")
		run_inference_checks(InferenceOrdinalAdjCatLogitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalOrderedProbitRegr")
		run_inference_checks(InferenceOrdinalOrderedProbitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalOrderedProbitRegr")
		run_inference_checks(InferenceOrdinalOrderedProbitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalStereotypeLogitRegr")
		run_inference_checks(InferenceOrdinalStereotypeLogitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalPropOddsRegr")
		run_inference_checks(InferenceOrdinalPropOddsRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalPartialProportionalOddsRegr")
		run_inference_checks(InferenceOrdinalPartialProportionalOddsRegr$new(des_obj, verbose = FALSE, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalPartialProportionalOddsRegr")
		run_inference_checks(InferenceOrdinalPartialProportionalOddsRegr$new(des_obj, verbose = FALSE, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalCloglogRegr")
		run_inference_checks(InferenceOrdinalCloglogRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalCloglogRegr")
		run_inference_checks(InferenceOrdinalCloglogRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalGCompMeanDiff")
		run_inference_checks(InferenceOrdinalGCompMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalGCompMeanDiff")
		run_inference_checks(InferenceOrdinalGCompMeanDiff$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalJonckheereTerpstraTest")
		run_inference_checks(InferenceOrdinalJonckheereTerpstraTest$new(des_obj), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalCauchitRegr")
		run_inference_checks(InferenceOrdinalCauchitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalCauchitRegr")
		run_inference_checks(InferenceOrdinalCauchitRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalContRatioRegr")
		run_inference_checks(InferenceOrdinalContRatioRegr$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
		inference_banner("InferenceOrdinalRidit")
		run_inference_checks(InferenceOrdinalRidit$new(des_obj, model_formula = model_formula), response_type, design_type, dataset_name, n_X, p_X)
	}

	run_exhaustive_remaining_inference_classes(des_obj, response_type, design_type, dataset_name, n_X, p_X, model_formula = model_formula)
}


for (dataset_name in names(datasets_and_response_models)){
	if (!is.na(DATASET_FILTER) && !identical(dataset_name, DATASET_FILTER)) {
		next
	}
	pending_dataset_header <<- paste0("\n\n====== DATASET: ", dataset_name, " ======\n")
	for (beta_T_iter_curr in seq_along(beta_T_values)){
		beta_T = beta_T_values[beta_T_iter_curr]
		if (!is.na(BETA_T_FILTER) && !identical(as.numeric(beta_T), BETA_T_FILTER)) {
			next
		}
		pending_beta_header <<- paste0("  === beta_T = [", beta_T, "] ===")
		for (rep_curr in 1:Nrep) {
			if (!is.na(REP_FILTER) && !identical(as.integer(rep_curr), REP_FILTER)) {
				next
			}
			pending_rep_header <<- paste0("    === rep ", rep_curr, " of ", Nrep, " ===")
			for (response_type in ALL_RESPONSE_TYPES) {
				if (!is.na(RESPONSE_TYPE_FILTER) && !identical(response_type, RESPONSE_TYPE_FILTER)) {
					next
				}
				if (!(response_type %in% names(datasets_and_response_models[[dataset_name]]$y_original))) {
					next
				}
				pending_response_header <<- paste0("      === response_type: ", response_type, " ===")

				formulas_to_test = list(~ 1, ~ .) #univariate then multivariate, leave as vector so we can edit later

				for (model_formula in formulas_to_test) {
					for (design_type in ALL_DESIGN_TYPES) {
						if (!is.na(DESIGN_TYPE_FILTER) && !identical(design_type, DESIGN_TYPE_FILTER)) {
							next
						}
						pending_design_header <<- paste0("        === design: ", design_type, " ===")
						tryCatch(
							run_tests_for_response(response_type, design_type = design_type, dataset_name = dataset_name, model_formula = model_formula),
							error = function(e){
								message("  FATAL ERROR in run_tests_for_response(", response_type, ", ", design_type, ", ", dataset_name, "): ", e$message)
							}
						)
					}
				}
			}
		}
	}
}
message("\n\n----------------------All tests complete!")
write_results_if_needed(force = TRUE)

unset_num_cores()
rm(list=ls())
