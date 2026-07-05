library(testthat)
library(EDI)

test_that("continuation-ratio augmentation matches an R reference", {
	X <- cbind(x1 = seq(-0.6, 0.6, length.out = 8L), x2 = rep(c(-0.25, 0.25), 4L))
	y <- c(10, 20, 30, 40, 40, 30, 20, 10)
	levels <- sort(unique(y))
	n_alpha <- length(levels) - 1L

	rows <- lapply(seq_along(y), function(i) {
		y_level <- match(y[i], levels)
		do.call(rbind, lapply(seq_len(min(y_level, n_alpha)), function(j) {
			c(
				as.numeric(seq_len(n_alpha) == j),
				X[i, ],
				z = as.numeric(y_level == j)
			)
		}))
	})
	reference <- do.call(rbind, rows)
	X_aug_reference <- unname(reference[, seq_len(ncol(reference) - 1L), drop = FALSE])
	z_reference <- unname(reference[, ncol(reference)])

	fit <- fast_continuation_ratio_regression_cpp(X, y, maxit = 1L)
	expect_identical(unname(fit$X_aug), X_aug_reference)
	expect_identical(as.numeric(fit$z), z_reference)

	params <- seq(-0.3, 0.3, length.out = n_alpha + ncol(X))
	eta <- drop(X_aug_reference %*% params)
	mu <- stats::plogis(pmax(-20, pmin(20, eta)))
	score_reference <- drop(crossprod(X_aug_reference, z_reference - mu))
	hessian_reference <- -crossprod(X_aug_reference, X_aug_reference * (mu * (1 - mu)))

	expect_equal(
		as.numeric(EDI:::get_continuation_ratio_regression_score_cpp(X, y, params)),
		score_reference,
		tolerance = 1e-12
	)
	expect_equal(
		unname(EDI:::get_continuation_ratio_regression_hessian_cpp(X, y, params)),
		hessian_reference,
		tolerance = 1e-12
	)
})
