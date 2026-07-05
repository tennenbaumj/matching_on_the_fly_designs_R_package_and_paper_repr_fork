library(testthat)
library(EDI)

mn_ci_r_bisection_reference <- function(x_t, n_t, x_c, n_c, p_t, p_c,
		alpha = 0.05, pval_epsilon = 1e-7) {
	estimate <- p_t - p_c
	pvalue <- function(delta) {
		EDI:::mn_pvalue_cpp(x_t, n_t, x_c, n_c, delta, p_t, p_c)
	}
	find_bound <- function(lower, upper) {
		for (iteration in seq_len(50L)) {
			middle <- (lower + upper) / 2
			p <- pvalue(middle)
			if (!is.finite(p)) p <- 0
			if (abs(p - alpha) < pval_epsilon) return(middle)
			if (p > alpha) {
				if (middle < estimate) upper <- middle else lower <- middle
			} else {
				if (middle < estimate) lower <- middle else upper <- middle
			}
		}
		(lower + upper) / 2
	}
	c(find_bound(-0.999999, estimate), find_bound(estimate, 0.999999))
}

test_that("Miettinen-Nurminen CI wrapper dispatches once to C++", {
	method_body <- paste(deparse(body(
		InferenceIncidMiettinenNurminenRiskDiff$public_methods$compute_asymp_confidence_interval
	)), collapse = "\n")
	calls <- regmatches(method_body, gregexpr("mn_ci_cpp", method_body, fixed = TRUE))[[1L]]
	expect_length(calls, 1L)
	expect_false(grepl("mn_pvalue_cpp", method_body, fixed = TRUE))
	expect_false(grepl("for\\s*\\(|while\\s*\\(|repeat\\s*\\{", method_body))
})

test_that("C++ Miettinen-Nurminen CI matches R-level bisection reference", {
	tables <- list(
		c(x_t = 60, n_t = 100, x_c = 40, n_c = 100),
		c(x_t = 12, n_t = 40, x_c = 20, n_c = 50),
		c(x_t = 95, n_t = 100, x_c = 70, n_c = 100)
	)
	for (counts in tables) {
		p_t <- counts[["x_t"]] / counts[["n_t"]]
		p_c <- counts[["x_c"]] / counts[["n_c"]]
		actual <- EDI:::mn_ci_cpp(
			counts[["x_t"]], counts[["n_t"]], counts[["x_c"]], counts[["n_c"]],
			p_t, p_c, 0.05, 1e-10
		)
		expected <- mn_ci_r_bisection_reference(
			counts[["x_t"]], counts[["n_t"]], counts[["x_c"]], counts[["n_c"]],
			p_t, p_c, 0.05, 1e-10
		)
		expect_equal(as.numeric(actual), expected, tolerance = 1e-14)
	}
})
