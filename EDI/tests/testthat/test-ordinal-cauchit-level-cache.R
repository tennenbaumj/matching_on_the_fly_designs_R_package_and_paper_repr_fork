test_that("ordinal cauchit cached level count preserves arbitrary labels", {
	set.seed(7201)
	n <- 180L
	X <- matrix(rnorm(n * 3L), ncol = 3L)
	y_index <- rep(c(1L, 2L, 3L), length.out = n)
	y_labels <- c(30, 10, 20)[y_index]
	y_remapped <- match(y_labels, sort(unique(y_labels)))

	fit_labels <- EDI:::fast_ordinal_cauchit_regression_cpp(X, y_labels)
	fit_remapped <- EDI:::fast_ordinal_cauchit_regression_cpp(X, y_remapped)

	expect_equal(fit_labels$n_params, 2L + ncol(X))
	expect_equal(as.numeric(fit_labels$params), as.numeric(fit_remapped$params), tolerance = 0)
	expect_identical(fit_labels$iterations, fit_remapped$iterations)
})
