library(testthat)
library(EDI)

test_that("simple Wilcox defaults away from BCa bootstrap", {
	expect_identical(
		EDI:::edi_bootstrap_dispatch_policy("InferenceAllSimpleWilcox"),
		"percentile"
	)
})
