library(testthat)
library(EDI)

expect_stratified_bootstrap_invariants <- function(strata) {
	indices <- EDI:::stratified_bootstrap_indices_cpp(strata)
	expect_length(indices, length(strata))
	expect_true(all(indices >= 1L & indices <= length(strata)))
	expect_identical(
		unname(table(strata[indices])),
		unname(table(strata))
	)
}

test_that("stratified bootstrap accepts integer IDs without string comparison", {
	strata_ids <- as.integer(rep(c(10L, 30L, 70L, 100L), c(15L, 25L, 10L, 20L)))
	expect_stratified_bootstrap_invariants(strata_ids)
})

test_that("stratified bootstrap retains character compatibility", {
	strata_keys <- rep(c("alpha", "beta", "gamma"), c(12L, 18L, 10L))
	expect_stratified_bootstrap_invariants(strata_keys)
})

test_that("stratified bootstrap is marginally uniform within strata", {
	strata_ids <- as.integer(rep(1:5, each = 20L))
	draws <- replicate(1000L, EDI:::stratified_bootstrap_indices_cpp(strata_ids))
	frequency <- tabulate(as.vector(draws), nbins = length(strata_ids))
	expect_lt(max(abs(frequency - 1000L)) / 1000, 0.15)
})
