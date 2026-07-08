#' Jackknife-based Inference
#'
#' Abstract class for delete-1 jackknife estimate correction and jackknife-Wald
#' inference layered on top of bootstrap-capable inference classes.
#'
#' @keywords internal
InferenceJackknife = R6::R6Class("InferenceJackknife",
	lock_objects = FALSE,
	inherit = InferenceBayesianBootstrap,
	public = list(
		#' @description Returns the leave-one-out jackknife estimate distribution.
		#'
		#' @param unit Deletion unit. Default `\"auto\"`, which chooses a design-aware
		#'   unit automatically.
		#'
		#' @return A numeric vector of jackknife replicate estimates.
		approximate_jackknife_distribution_beta_hat_T = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			private$assert_jackknife_supported(unit = unit)
			as.numeric(private$approximate_jackknife_distribution_beta_hat_T_private(unit = unit))
		},
		#' @description Computes the delete-1 jackknife bias-corrected treatment estimate.
		#'
		#' For blocking designs, this uses leave-one-block-out deletion units. For
		#' matching designs, it uses leave-match-out deletion units. For KK designs,
		#' it uses leave-match-out for matched pairs and leave-one-out for reservoir
		#' subjects.
		#'
		#' @param unit Deletion unit. Default `\"auto\"`, which chooses a design-aware
		#'   unit automatically.
		#'
		#' @return A numeric jackknife bias-corrected treatment estimate.
		compute_jackknife_estimate = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			private$assert_jackknife_supported(unit = unit)
			private$compute_jackknife_summary(unit = unit)$estimate
		},
		#' @description Alias for `compute_jackknife_estimate()`.
		#'
		#' @param unit Deletion unit. Default `\"auto\"`.
		#'
		#' @return A numeric jackknife bias-corrected treatment estimate.
		compute_jackknife_corrected_estimate = function(unit = "auto"){
			self$compute_jackknife_estimate(unit = unit)
		},
		#' @description Computes the jackknife bias estimate.
		#'
		#' @param unit Deletion unit. Default `\"auto\"`.
		#'
		#' @return A numeric jackknife bias estimate.
		compute_jackknife_bias_estimate = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			private$assert_jackknife_supported(unit = unit)
			private$compute_jackknife_summary(unit = unit)$bias
		},
		#' @description Computes the delete-1 jackknife standard error.
		#'
		#' For blocking designs, this uses leave-one-block-out deletion units. For
		#' matching designs, it uses leave-match-out deletion units. For KK designs,
		#' it uses leave-match-out for matched pairs and leave-one-out for reservoir
		#' subjects.
		#'
		#' @param unit Deletion unit. Default `\"auto\"`, which chooses a design-aware
		#'   unit automatically.
		#'
		#' @return A numeric jackknife standard error.
		compute_jackknife_std_error = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			private$assert_jackknife_supported(unit = unit)
			private$compute_jackknife_summary(unit = unit)$std_error
		},
		#' @description Alias for `compute_jackknife_std_error()`.
		#'
		#' @param unit Deletion unit. Default `\"auto\"`.
		#'
		#' @return A numeric jackknife standard error.
		compute_jackknife_standard_error = function(unit = "auto"){
			self$compute_jackknife_std_error(unit = unit)
		},
		#' @description Computes a two-sided Wald p-value using the jackknife estimate
		#'   and jackknife standard error.
		#'
		#' For blocking designs, this uses leave-one-block-out deletion units. For
		#' matching designs, it uses leave-match-out deletion units. For KK designs,
		#' it uses leave-match-out for matched pairs and leave-one-out for reservoir
		#' subjects.
		#'
		#' @param delta Null treatment-effect value. Default 0.
		#' @param unit Deletion unit. Default `\"auto\"`, which chooses a design-aware
		#'   unit automatically.
		#'
		#' @return A two-sided jackknife-Wald p-value.
		compute_jackknife_wald_two_sided_pval = function(delta = 0, unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			private$assert_jackknife_supported(unit = unit)
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
			}
			jack_summary = private$compute_jackknife_summary(unit = unit)
			se_j = jack_summary$std_error
			if (!is.finite(jack_summary$estimate)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_se("jackknife_estimate_unavailable")
				return(NA_real_)
			}
			if (!is.finite(se_j) || se_j <= 0) {
				if (isTRUE(private$harden)) private$cache_nonestimable_se("jackknife_standard_error_unavailable")
				return(NA_real_)
			}
			theta_hat = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
			if (!is.finite(theta_hat)) return(NA_real_)
			2 * stats::pnorm(-abs((theta_hat - delta) / se_j))
		},
		#' @description Computes a normal-approximation confidence interval using the
		#'   jackknife estimate and jackknife standard error.
		#'
		#' For blocking designs, this uses leave-one-block-out deletion units. For
		#' matching designs, it uses leave-match-out deletion units. For KK designs,
		#' it uses leave-match-out for matched pairs and leave-one-out for reservoir
		#' subjects.
		#'
		#' @param alpha Significance level. Default 0.05.
		#' @param unit Deletion unit. Default `\"auto\"`, which chooses a design-aware
		#'   unit automatically.
		#'
		#' @return A jackknife-Wald confidence interval.
		compute_jackknife_wald_confidence_interval = function(alpha = 0.05, unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			private$assert_jackknife_supported(unit = unit)
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			jack_summary = private$compute_jackknife_summary(unit = unit)
			se_j = jack_summary$std_error
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			if (!is.finite(jack_summary$estimate)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_se("jackknife_estimate_unavailable")
				return(ci)
			}
			if (!is.finite(se_j) || se_j <= 0) {
				if (isTRUE(private$harden)) private$cache_nonestimable_se("jackknife_standard_error_unavailable")
				return(ci)
			}
			theta_hat = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
			if (!is.finite(theta_hat)) return(ci)
			z = stats::qnorm(1 - alpha / 2)
			ci[] = c(theta_hat - z * se_j, theta_hat + z * se_j)
			ci
		}
	),
	private = list(
		normalize_jackknife_unit = function(unit = "auto"){
			unit = tolower(as.character(unit)[1L])
			valid = c("auto", "observation", "cluster", "block", "pair", "matched_set")
			if (should_run_asserts()) {
				assertChoice(unit, valid)
			}
			unit
		},
		assert_jackknife_supported = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			design_obj = private$des_obj
			is_blocking_design = is(design_obj, "DesignBlocking") &&
				isTRUE(tryCatch(design_obj$is_blocking_design(), error = function(e) FALSE))
			is_matching_design = is(design_obj, "DesignMatching") &&
				isTRUE(tryCatch(design_obj$is_matching_design(), error = function(e) FALSE))
			is_cluster_design = is(design_obj, "DesignFixedCluster") || is(design_obj, "DesignFixedBlockedCluster")
			if (should_run_asserts()) {
				if (unit == "block" && !is_blocking_design) {
					stop("jackknife unit = 'block' requires a blocking design.", call. = FALSE)
				}
				if (unit %in% c("pair", "matched_set") && !is_matching_design) {
					stop("jackknife unit = '", unit, "' requires a matching design.", call. = FALSE)
				}
				if (unit == "cluster" && !is_cluster_design) {
					stop("jackknife unit = 'cluster' requires a clustered design.", call. = FALSE)
				}
			}
			invisible(NULL)
		},
		jackknife_cache_key = function(unit = "auto"){
			private$normalize_jackknife_unit(unit)
		},
		compute_jackknife_summary = function(unit = "auto"){
			unit = private$normalize_jackknife_unit(unit)
			cache_key = private$jackknife_cache_key(unit)
			if (is.null(private$cached_values$jackknife_summary)) {
				private$cached_values$jackknife_summary = list()
			}
			if (!is.null(private$cached_values$jackknife_summary[[cache_key]])) {
				return(private$cached_values$jackknife_summary[[cache_key]])
			}
			jack = as.numeric(private$approximate_jackknife_distribution_beta_hat_T_private(unit = unit))
			n_units = length(jack)
			if (n_units <= 1L) {
				summary = list(estimate = NA_real_, bias = NA_real_, std_error = NA_real_, distribution = numeric(0))
				private$cached_values$jackknife_summary[[cache_key]] = summary
				return(summary)
			}
			theta_hat = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
			if (!is.finite(theta_hat) || any(!is.finite(jack))) {
				summary = list(estimate = NA_real_, bias = NA_real_, std_error = NA_real_, distribution = jack)
				private$cached_values$jackknife_summary[[cache_key]] = summary
				return(summary)
			}
			jack_bar = mean(jack)
			bias_j = (n_units - 1) * (jack_bar - theta_hat)
			theta_j = theta_hat - bias_j
			var_j = ((n_units - 1) / n_units) * sum((jack - jack_bar)^2)
			se_j = if (is.finite(var_j) && var_j >= 0) sqrt(var_j) else NA_real_
			summary = list(estimate = theta_j, bias = bias_j, std_error = se_j, distribution = jack)
			private$cached_values$jackknife_summary[[cache_key]] = summary
			summary
		}
	)
)
