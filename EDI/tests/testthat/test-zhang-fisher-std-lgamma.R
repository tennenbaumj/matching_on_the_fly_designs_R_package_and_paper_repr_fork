library(testthat)
library(EDI)

test_that("Zhang Fisher p-value matches fisher.test across odds ratios", {
	cases <- list(
		c(n11 = 30, n10 = 10, n01 = 5, n00 = 55, log_or = 0),
		c(n11 = 12, n10 = 18, n01 = 9, n00 = 21, log_or = log(1.7)),
		c(n11 = 1, n10 = 8, n01 = 6, n00 = 15, log_or = log(0.6)),
		c(n11 = 40, n10 = 5, n01 = 3, n00 = 52, log_or = log(2.5))
	)

	for (counts in cases) {
		table_2x2 <- matrix(counts[c("n11", "n10", "n01", "n00")], nrow = 2L)
		expected <- stats::fisher.test(
			table_2x2,
			or = exp(counts[["log_or"]])
		)$p.value
		actual <- EDI:::zhang_exact_fisher_pval_cpp(
			counts[["n11"]], counts[["n10"]],
			counts[["n01"]], counts[["n00"]],
			counts[["log_or"]]
		)
		expect_equal(actual, expected, tolerance = 5e-14)
	}
})

test_that("Zhang Fisher p-value handles degenerate and invalid tables", {
	expect_equal(EDI:::zhang_exact_fisher_pval_cpp(0L, 0L, 0L, 20L, 0), 1)
	expect_true(is.na(EDI:::zhang_exact_fisher_pval_cpp(-1L, 2L, 3L, 4L, 0)))
	expect_true(is.na(EDI:::zhang_exact_fisher_pval_cpp(1L, 2L, 3L, 4L, Inf)))
})
