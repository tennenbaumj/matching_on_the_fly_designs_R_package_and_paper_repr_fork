library(testthat)
library(EDI)

wrap_likelihood_spec_with_counters <- function(inf){
	priv <- inf$.__enclos_env__$private
	orig_get_spec <- priv$get_likelihood_test_spec
	counts <- new.env(parent = emptyenv())
	counts$fit_null <- 0L
	counts$score <- 0L
	counts$neg_loglik <- 0L
	counts$information <- 0L

	if (bindingIsLocked("get_likelihood_test_spec", priv)) {
		get("unlockBinding", envir = asNamespace("base"))("get_likelihood_test_spec", priv)
	}
	priv[["get_likelihood_test_spec"]] <- function(){
		spec <- orig_get_spec()
		if (is.null(spec)) return(NULL)

		orig_fit_null <- spec$fit_null
		orig_score <- spec$score
		orig_neg_loglik <- spec$neg_loglik
		orig_information <- spec$information
		orig_observed_information <- spec$observed_information
		orig_fisher_information <- spec$fisher_information
		fit_null_formals <- tryCatch(names(formals(orig_fit_null)), error = function(e) character())
		accepts_start <- "start" %in% fit_null_formals

		spec$fit_null <- function(delta, start = NULL){
			counts$fit_null <- counts$fit_null + 1L
			if (accepts_start) {
				orig_fit_null(delta, start = start)
			} else {
				orig_fit_null(delta)
			}
		}
		spec$score <- function(fit){
			counts$score <- counts$score + 1L
			orig_score(fit)
		}
		spec$neg_loglik <- function(fit){
			counts$neg_loglik <- counts$neg_loglik + 1L
			orig_neg_loglik(fit)
		}
		if (!is.null(orig_information)) {
			spec$information <- function(fit){
				counts$information <- counts$information + 1L
				orig_information(fit)
			}
		}
		if (!is.null(orig_observed_information)) {
			spec$observed_information <- function(fit){
				counts$information <- counts$information + 1L
				orig_observed_information(fit)
			}
		}
		if (!is.null(orig_fisher_information)) {
			spec$fisher_information <- function(fit){
				counts$information <- counts$information + 1L
				orig_fisher_information(fit)
			}
		}
		spec
	}

	list(private = priv, counts = counts)
}

test_that("likelihood-test p-values memoize null fits and score-side components by testing_type and delta", {
	set.seed(111)
	n <- 80
	x <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	des <- EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(rbinom(n, 1, plogis(-0.25 + 0.5 * ((w+1)/2) + 0.35 * x)))

	inf <- InferenceIncidLogRegr$new(des, verbose = FALSE, smart_cold_start_default = TRUE)
	wrapped <- wrap_likelihood_spec_with_counters(inf)

	p1 <- inf$compute_score_two_sided_pval(0.1)
	counts_after_first <- as.list.environment(wrapped$counts, all.names = TRUE)
	p2 <- inf$compute_score_two_sided_pval(0.1)
	counts_after_second <- as.list.environment(wrapped$counts, all.names = TRUE)

	expect_true(is.finite(p1))
	expect_equal(p1, p2, tolerance = 0)
	expect_equal(counts_after_first$fit_null, 1L)
	expect_equal(counts_after_first$score, 1L)
	expect_equal(counts_after_first$information, 1L)
	expect_equal(counts_after_second, counts_after_first)

	cache_key <- wrapped$private$likelihood_test_delta_key("score", 0.1)
	cache_entry <- wrapped$private$cached_values$likelihood_test_eval_cache[[cache_key]]
	expect_true(!is.null(cache_entry))
	expect_false(isTRUE(cache_entry$invalid))
	expect_true(!is.null(cache_entry$null_fit))
	expect_true(!is.null(cache_entry$score))
	expect_true(!is.null(cache_entry$information))
	expect_true(is.finite(cache_entry$p_value))
})

test_that("likelihood-ratio CI inversion reuses memoized null fits and neg-loglik values", {
	set.seed(222)
	n <- 90
	x <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	des <- EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(rbinom(n, 1, plogis(-0.15 + 0.45 * ((w+1)/2) + 0.25 * x)))

	inf <- InferenceIncidLogRegr$new(des, verbose = FALSE, smart_cold_start_default = TRUE)
	wrapped <- wrap_likelihood_spec_with_counters(inf)

	ci1 <- inf$compute_lik_ratio_confidence_interval(alpha = 0.2)
	counts_after_first <- as.list.environment(wrapped$counts, all.names = TRUE)
	ci2 <- inf$compute_lik_ratio_confidence_interval(alpha = 0.2)
	counts_after_second <- as.list.environment(wrapped$counts, all.names = TRUE)

	expect_true(all(is.finite(ci1)))
	expect_equal(ci1, ci2, tolerance = 0)
	expect_gt(counts_after_first$fit_null, 0L)
	expect_gt(counts_after_first$neg_loglik, 0L)
	expect_gt(counts_after_first$score, 0L)
	expect_equal(counts_after_second$fit_null, counts_after_first$fit_null)
	expect_equal(counts_after_second$neg_loglik, counts_after_first$neg_loglik)
	expect_equal(counts_after_second$score, counts_after_first$score)
})

test_that("KK OLS gradient CI reuses memoized null fits on repeated evaluation", {
	des <- DesignFixedBinaryMatch$new(response_type = "continuous", n = 6, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = c(-2, -1, 0, 1, 2, 3)))
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(c(-1, 0, 1, 1, 4, 9))

	inf <- InferenceContinKKOLSOneLik$new(des, verbose = FALSE)
	wrapped <- wrap_likelihood_spec_with_counters(inf)

	ci1 <- inf$compute_gradient_confidence_interval(alpha = 0.2)
	counts_after_first <- as.list.environment(wrapped$counts, all.names = TRUE)
	ci2 <- inf$compute_gradient_confidence_interval(alpha = 0.2)
	counts_after_second <- as.list.environment(wrapped$counts, all.names = TRUE)

	expect_true(all(is.finite(ci1)))
	expect_equal(ci1, ci2, tolerance = 0)
	expect_gt(counts_after_first$fit_null, 0L)
	expect_gt(counts_after_first$score, 0L)
	expect_equal(counts_after_second, counts_after_first)
})
