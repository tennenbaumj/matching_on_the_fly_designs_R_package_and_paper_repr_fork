run_likelihood_method_smoke_suite <- function(response_type_filter = NA_character_){
	should_run = function(response_type) is.na(response_type_filter) || response_type == response_type_filter
	old_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
		get(".Random.seed", envir = .GlobalEnv)
	} else {
		NULL
	}
	on.exit({
		if (!is.null(old_seed)) {
			assign(".Random.seed", old_seed, envir = .GlobalEnv)
		} else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
			rm(".Random.seed", envir = .GlobalEnv)
		}
	}, add = TRUE)
	set.seed(20240508)

	call_all_methods <- function(inf, label){
		is_unsupported_method_error <- function(err){
			msg = conditionMessage(err)
			grepl("does not expose a likelihood-test specification", msg, fixed = TRUE) ||
				grepl("does not support score p-values", msg, fixed = TRUE) ||
				grepl("does not support score confidence intervals", msg, fixed = TRUE) ||
				grepl("does not support gradient p-values", msg, fixed = TRUE) ||
				grepl("does not support gradient confidence intervals", msg, fixed = TRUE) ||
				grepl("does not support likelihood-ratio p-values", msg, fixed = TRUE) ||
				grepl("does not support likelihood-ratio confidence intervals", msg, fixed = TRUE) ||
				grepl("does not support parametric-bootstrap LR calibration", msg, fixed = TRUE)
		}

		pval_methods <- c(
			"compute_wald_two_sided_pval",
			"compute_score_two_sided_pval",
			"compute_gradient_two_sided_pval",
			"compute_lik_ratio_two_sided_pval"
		)
		ci_methods <- c(
			"compute_wald_confidence_interval",
			"compute_score_confidence_interval",
			"compute_gradient_confidence_interval",
			"compute_lik_ratio_confidence_interval"
		)

		for (method_name in pval_methods){
			method_fn = inf[[method_name]]
			stopifnot(is.function(method_fn))
			cat(sprintf("  [%s] calling %s ...\n", label, method_name))
			val = tryCatch(
				method_fn(delta = 0),
				error = function(e) {
					if (is_unsupported_method_error(e)) {
						cat(sprintf("  [%s] skipping %s: %s\n", label, method_name, conditionMessage(e)))
						return(NULL)
					}
					stop(label, ": ", method_name, " failed: ", e$message)
				}
			)
			if (is.null(val)) next
			if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val < 0 || val > 1) {
				cat(sprintf("  [%s] skipping %s: non-finite or invalid p-value output\n", label, method_name))
				next
			}
			cat(sprintf("  [%s] %s = %s\n", label, method_name, paste(val, collapse=", ")))
		}

		for (method_name in ci_methods){
			method_fn = inf[[method_name]]
			stopifnot(is.function(method_fn))
			cat(sprintf("  [%s] calling %s ...\n", label, method_name))
			val = tryCatch(
				method_fn(alpha = 0.2),
				error = function(e) {
					if (is_unsupported_method_error(e)) {
						cat(sprintf("  [%s] skipping %s: %s\n", label, method_name, conditionMessage(e)))
						return(NULL)
					}
					stop(label, ": ", method_name, " failed: ", e$message)
				}
			)
			if (is.null(val)) next
			if (!is.numeric(val) || length(val) != 2L || !all(is.finite(val)) || val[1] > val[2]) {
				cat(sprintf("  [%s] skipping %s: non-finite or invalid confidence interval output\n", label, method_name))
				next
			}
			cat(sprintf("  [%s] %s = %s\n", label, method_name, paste(val, collapse=", ")))
		}

		supports_param_boot = is(inf, "InferenceParamBootstrap") &&
			isTRUE(tryCatch(
				inf$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(),
				error = function(e) FALSE
			))
		if (supports_param_boot){
			param_boot_pval = tryCatch(
				inf$compute_lik_ratio_bootstrap_two_sided_pval(delta = 0, B = 5L, show_progress = FALSE),
				error = function(e) {
					if (is_unsupported_method_error(e)) {
						cat(sprintf("  [%s] skipping compute_lik_ratio_bootstrap_two_sided_pval: %s\n", label, conditionMessage(e)))
						return(NULL)
					}
					stop(label, ": compute_lik_ratio_bootstrap_two_sided_pval failed: ", e$message)
				}
			)
			if (!is.null(param_boot_pval)) {
				if (!is.numeric(param_boot_pval) || length(param_boot_pval) != 1L || !is.finite(param_boot_pval) || param_boot_pval < 0 || param_boot_pval > 1) {
					cat(sprintf("  [%s] skipping compute_lik_ratio_bootstrap_two_sided_pval: non-finite or invalid p-value output\n", label))
					return(invisible(TRUE))
				}
				cat(sprintf("  [%s] compute_lik_ratio_bootstrap_two_sided_pval = %s\n", label, paste(param_boot_pval, collapse = ", ")))
			}

			param_boot_ci = tryCatch(
				inf$compute_lik_ratio_bootstrap_confidence_interval(alpha = 0.2, B = 5L, show_progress = FALSE),
				error = function(e) {
					if (is_unsupported_method_error(e)) {
						cat(sprintf("  [%s] skipping compute_lik_ratio_bootstrap_confidence_interval: %s\n", label, conditionMessage(e)))
						return(NULL)
					}
					stop(label, ": compute_lik_ratio_bootstrap_confidence_interval failed: ", e$message)
				}
			)
			if (!is.null(param_boot_ci)) {
				if (!is.numeric(param_boot_ci) || length(param_boot_ci) != 2L || !all(is.finite(param_boot_ci)) || param_boot_ci[1] > param_boot_ci[2]) {
					cat(sprintf("  [%s] skipping compute_lik_ratio_bootstrap_confidence_interval: non-finite or invalid confidence interval output\n", label))
					return(invisible(TRUE))
				}
				cat(sprintf("  [%s] compute_lik_ratio_bootstrap_confidence_interval = %s\n", label, paste(param_boot_ci, collapse = ", ")))
			}
		}

		invisible(TRUE)
	}

	make_fixed_incidence_design <- function(n = 40L){
		x1 = rnorm(n)
		des = DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
		des$add_all_subjects_to_experiment(data.frame(x1 = x1))
		des$assign_w_to_all_subjects()
		w = des$get_w()
		y = rbinom(n, 1, plogis(-0.3 + 0.6 * w + 0.4 * x1))
		des$add_all_subject_responses(y)
		des
	}

	make_fixed_count_design <- function(n = 40L){
		x1 = rnorm(n)
		des = DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
		des$add_all_subjects_to_experiment(data.frame(x1 = x1))
		des$assign_w_to_all_subjects()
		w = des$get_w()
		y = rpois(n, exp(0.2 + 0.3 * w + 0.25 * x1))
		des$add_all_subject_responses(y)
		des
	}

	make_fixed_survival_design <- function(n = 40L){
		x1 = rnorm(n)
		des = DesignFixedBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
		des$add_all_subjects_to_experiment(data.frame(x1 = x1))
		des$assign_w_to_all_subjects()
		w = des$get_w()
		y = rexp(n, rate = exp(-0.2 + 0.15 * w + 0.2 * x1))
		des$add_all_subject_responses(y, rep(1L, n))
		des
	}

	make_fixed_ordinal_design <- function(n = 40L){
		x1 = rnorm(n)
		x2 = rnorm(n)
		des = DesignFixedBernoulli$new(n = n, response_type = "ordinal", verbose = FALSE)
		des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
		des$assign_w_to_all_subjects()
		w = des$get_w()
		eta = 0.45 * w + 0.35 * x1 - 0.20 * x2
		cut_1 = plogis(-1.0 - eta)
		cut_2 = plogis(0.2 - eta)
		cut_3 = plogis(1.1 - eta)
		u = runif(n)
		y = ifelse(u <= cut_1, 1L, ifelse(u <= cut_2, 2L, ifelse(u <= cut_3, 3L, 4L)))
		des$add_all_subject_responses(y)
		des
	}

	make_kk_incidence_design <- function(n = 16L){
		des = DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
		x1 = rnorm(n)
		x2 = rnorm(n)
		for (i in seq_len(n)){
			w_i = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
			y_i = rbinom(1, 1, plogis(-0.15 + 0.5 * w_i + 0.25 * x1[i] - 0.1 * x2[i]))
			des$add_one_subject_response(i, y_i, 1L)
		}
		des
	}

	make_kk_survival_design <- function(n = 16L){
		des = DesignSeqOneByOneKK14$new(n = n, response_type = "survival", verbose = FALSE)
		x1 = rnorm(n)
		x2 = rnorm(n)
		for (i in seq_len(n)){
			w_i = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
			y_i = rexp(1, rate = exp(-0.1 + 0.2 * w_i + 0.15 * x1[i] - 0.05 * x2[i]))
			des$add_one_subject_response(i, y_i, 1L)
		}
		des
	}

	make_kk_ordinal_design <- function(n = 32L){
		des = DesignSeqOneByOneKK14$new(n = n, response_type = "ordinal", verbose = FALSE)
		x1 = rnorm(n)
		x2 = rnorm(n)
		for (i in seq_len(n)){
			w_i = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
			y_i = sample.int(4L, 1L, prob = c(
				plogis(-0.8 + 0.4 * w_i + 0.2 * x1[i] - 0.1 * x2[i]),
				0.2,
				0.2,
				0.4
			))
			des$add_one_subject_response(i, y_i, 1L)
		}
		des
	}

	results = list()
	if (should_run("count")) {
		results$count_poisson = call_all_methods(
			InferenceCountPoisson$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceCountPoisson"
		)
		results$count_robust_poisson = call_all_methods(
			InferenceCountRobustPoisson$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceCountRobustPoisson"
		)
		results$count_quasi_poisson = call_all_methods(
			InferenceCountQuasiPoisson$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceCountQuasiPoisson"
		)
		results$count_zip = call_all_methods(
			InferenceCountZeroInflatedPoisson$new(make_fixed_count_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE, optimization_alg = "lbfgs"),
			"InferenceCountZeroInflatedPoisson"
		)
		results$count_hurdle_poisson = call_all_methods(
			InferenceCountHurdlePoisson$new(make_fixed_count_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE, optimization_alg = "lbfgs"),
			"InferenceCountHurdlePoisson"
		)
		results$count_hurdle_negbin = call_all_methods(
			InferenceCountHurdleNegBin$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceCountHurdleNegBin"
		)
	}
	if (should_run("incidence")) {
		results$incidence_modified_poisson = call_all_methods(
			InferenceIncidModifiedPoisson$new(make_fixed_incidence_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceIncidModifiedPoisson"
		)
		results$incidence_kk_modified_poisson = call_all_methods(
			InferenceIncidKKModifiedPoisson$new(make_kk_incidence_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceIncidKKModifiedPoisson"
		)
	}
	if (should_run("ordinal")) {
		results$ordinal_prop_odds = call_all_methods(
			InferenceOrdinalPropOddsRegr$new(make_fixed_ordinal_design(), model_formula = ~ x1 + x2, verbose = FALSE),
			"InferenceOrdinalPropOddsRegr"
		)
		results$ordinal_ordered_probit = call_all_methods(
			InferenceOrdinalOrderedProbitRegr$new(make_fixed_ordinal_design(), model_formula = ~ x1 + x2, verbose = FALSE),
			"InferenceOrdinalOrderedProbitRegr"
		)
		results$ordinal_cauchit = call_all_methods(
			InferenceOrdinalCauchitRegr$new(make_fixed_ordinal_design(), model_formula = ~ x1 + x2, verbose = FALSE),
			"InferenceOrdinalCauchitRegr"
		)
		results$ordinal_cloglog = call_all_methods(
			InferenceOrdinalCloglogRegr$new(make_fixed_ordinal_design(), model_formula = ~ x1 + x2, verbose = FALSE),
			"InferenceOrdinalCloglogRegr"
		)
		results$ordinal_kk_glmm = call_all_methods(
			InferenceOrdinalKKGLMM$new(make_kk_ordinal_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE),
			"InferenceOrdinalKKGLMM"
		)
	}
	if (should_run("survival")) {
		results$survival_cox = call_all_methods(
			InferenceSurvivalCoxPHRegr$new(make_fixed_survival_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE),
			"InferenceSurvivalCoxPHRegr"
		)
		results$survival_weibull = call_all_methods(
			InferenceSurvivalWeibullRegr$new(make_fixed_survival_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceSurvivalWeibullRegr"
		)
		results$survival_strat_cox = call_all_methods(
			InferenceSurvivalStratCoxPHRegr$new(make_fixed_survival_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE),
			"InferenceSurvivalStratCoxPHRegr"
		)
		results$kk_survival_strat_cox = call_all_methods(
			InferenceSurvivalKKStratCoxPHOneLik$new(make_kk_survival_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceSurvivalKKStratCoxPHOneLik"
		)
		results$kk_survival_lwa_cox = call_all_methods(
			InferenceSurvivalKKLWACoxPHOneLik$new(make_kk_survival_design(), model_formula = ~ x1, verbose = FALSE),
			"InferenceSurvivalKKLWACoxPHOneLik"
		)
		results$kk_survival_clayton = call_all_methods(
			InferenceSurvivalKKClaytonCopulaOneLik$new(make_kk_survival_design(n = 64L), model_formula = ~ x1, verbose = FALSE),
			"InferenceSurvivalKKClaytonCopulaOneLik"
		)
	}
	invisible(results)
}
