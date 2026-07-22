library(testthat)
library(EDI)

mixin_slot_names = function(mixin_names, slot){
	as.character(unlist(lapply(mixin_names, function(mixin_name) {
		mixin = get(mixin_name, envir = asNamespace("EDI"))
		names(mixin[[slot]])
	}), use.names = FALSE))
}

test_that("every mixin has a documented host contract and is collated after the registry", {
	contracts = EDI:::EDI_MIXIN_CONTRACTS
	mixin_names = ls(asNamespace("EDI"), pattern = "^InferenceMixin")
	expect_setequal(names(contracts), mixin_names)

	for (mixin_name in names(contracts)) {
		contract = contracts[[mixin_name]]
		expect_named(contract, c("file", "private_methods", "private_state"))
		expect_true(is.character(contract$file) && length(contract$file) == 1L)
		expect_true(is.character(contract$private_methods))
		expect_true(is.character(contract$private_state))
	}

	collate = strsplit(utils::packageDescription("EDI")$Collate, "[[:space:]]+")[[1L]]
	collate = gsub("'", "", collate, fixed = TRUE)
	registry_position = match("mixin_contracts.R", collate)
	expect_true(is.finite(registry_position))
	for (file in vapply(contracts, `[[`, character(1), "file")) {
		expect_gt(match(file, collate), registry_position)
	}
})

test_that("mixin composition has no undocumented method-name collisions", {
	for (target in names(EDI:::EDI_MIXIN_COMPOSITIONS)) {
		mixins = EDI:::EDI_MIXIN_COMPOSITIONS[[target]]
		allowed = EDI:::EDI_MIXIN_ALLOWED_COLLISIONS[[target]]
		if (is.null(allowed)) allowed = list(public = character(), private = character())
		for (slot in c("public", "private")) {
			methods = mixin_slot_names(mixins, slot)
			collisions = sort(unique(methods[duplicated(methods)]))
			expect_setequal(collisions, allowed[[slot]])
		}
	}
})

test_that("single-host protected bases are not decomposed into new mixins", {
	expect_false("InferenceCountZeroAugmentedPoissonAbstract" %in% names(EDI:::EDI_MIXIN_COMPOSITIONS))
	expect_false("InferenceMixinHurdlePoissonClosedForm" %in% names(EDI:::EDI_MIXIN_CONTRACTS))
	expect_false(exists("InferenceMixinHurdlePoissonClosedForm", envir = asNamespace("EDI"), inherits = FALSE))
})
