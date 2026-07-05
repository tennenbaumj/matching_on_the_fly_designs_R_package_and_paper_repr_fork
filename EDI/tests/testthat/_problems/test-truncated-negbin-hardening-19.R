# Extracted from test-truncated-negbin-hardening.R:19

# test -------------------------------------------------------------------------
set.seed(101)
X_cov <- matrix(rnorm(120), ncol = 2)
X <- cbind(1, X_cov)
lambda <- exp(0.3 + X_cov[, 1] * 0.2 - X_cov[, 2] * 0.1)
y <- pmax(rnbinom(nrow(X), mu = lambda, size = 3), 1)
expect_error(
		fast_truncated_negbin_count_cpp(X, y, warm_start_params = rep(0, ncol(X))),
		"warm_start_params"
	)
expect_error(
		fast_truncated_negbin_count_cpp(X, y, warm_start_fisher_info = diag(ncol(X))),
		"warm_start_fisher_info"
	)
expect_error(
		fast_truncated_negbin_count_cpp(X, c(y[-1], 1.5)),
		"integer-valued counts"
	)
