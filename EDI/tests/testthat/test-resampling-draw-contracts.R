test_that("resampling draw contracts define comparable hooks", {
	ops = c("rand", "non_param_boot", "rand_bootstrap", "bayesian_boot")
	expect_setequal(names(EDI_RESAMPLING_DRAW_CONTRACTS), ops)

	required_fields = c(
		"operation",
		"draw_type",
		"loader",
		"estimator",
		"cache_name",
		"cache_key_method"
	)
	contracts = lapply(ops, resampling_draw_contract)
	for (i in seq_along(contracts)) {
		expect_named(contracts[[i]], required_fields)
		expect_identical(contracts[[i]]$operation, ops[[i]])
		expect_true(nzchar(contracts[[i]]$draw_type))
		expect_true(nzchar(contracts[[i]]$loader))
		expect_true(nzchar(contracts[[i]]$estimator))
		expect_true(nzchar(contracts[[i]]$cache_name))
	}
	expect_equal(
		unname(vapply(contracts, `[[`, character(1), "cache_name")),
		c("rand_distr_cache", "boot_distr_cache", "rand_boot_distr_cache", "bayes_boot_distr_cache")
	)
})

test_that("unknown resampling operations fail clearly", {
	expect_error(
		resampling_draw_contract("not_an_operation"),
		"Unknown resampling operation: not_an_operation",
		fixed = TRUE
	)
	expect_error(
		resampling_draw_contract(c("rand", "non_param_boot")),
		"operation must be one resampling operation name",
		fixed = TRUE
	)
})

test_that("distribution cache helpers route through the contract cache slot", {
	cached_values = list(
		rand_distr_cache = list(existing = 1),
		unrelated = TRUE
	)

	cached_values = resampling_distribution_cache_set(
		cached_values,
		"rand_bootstrap",
		"k",
		c(1, 2, 3)
	)
	expect_equal(
		resampling_distribution_cache_get(cached_values, "rand_bootstrap", "k"),
		c(1, 2, 3)
	)
	expect_equal(cached_values$rand_distr_cache$existing, 1)
	expect_true(cached_values$unrelated)

	cached_values = resampling_distribution_cache_set(
		cached_values,
		"rand_bootstrap",
		"k",
		NULL
	)
	expect_null(resampling_distribution_cache_get(cached_values, "rand_bootstrap", "k"))
})

test_that("distribution cache ensure initializes only the requested cache", {
	cached_values = resampling_distribution_cache_ensure(list(), "bayesian_boot")
	expect_named(cached_values, "bayes_boot_distr_cache")
	expect_equal(cached_values$bayes_boot_distr_cache, list())
})
