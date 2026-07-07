library(EDI)
library(testthat)

# Helpers
make_ordinal_data <- function(n, p, K, seed = 42) {
    set.seed(seed)
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    beta <- rnorm(p)
    eta <- X %*% beta
    cuts <- quantile(eta, probs = seq(0, 1, length.out = K + 1))[2:K]
    y <- as.integer(cut(eta + rnorm(n, sd = 0.5), c(-Inf, cuts, Inf)))
    list(X = X, y = y)
}

# Test 1: basic fit matches score near zero
test_that("gradient near zero at MLE (K=4)", {
    d <- make_ordinal_data(200, 2, 4)
    fit <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    expect_true(fit$converged)
    score <- EDI:::get_stereotype_logit_score_cpp(d$X, d$y, fit$params)
    expect_lt(max(abs(score)), 1e-3)
})

# Test 2: hessian is negative definite at MLE (PSD fisher info, up to numerical noise)
test_that("fisher information is positive semidefinite (K=4)", {
    d <- make_ordinal_data(200, 2, 4)
    fit <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    expect_true(fit$converged)
    ev <- eigen(fit$fisher_information, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > -1e-3))
})

# Test 3: repeated calls give identical results (no side-effect from member field reuse)
test_that("repeated calls give identical results", {
    d <- make_ordinal_data(150, 2, 3)
    fit1 <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    fit2 <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    expect_equal(fit1$b, fit2$b)
    expect_equal(fit1$params, fit2$params)
    expect_equal(fit1$fisher_information, fit2$fisher_information)
})

# Test 4: K=3 (n_gamma=1) works correctly
test_that("K=3 fit converges and score near zero", {
    d <- make_ordinal_data(200, 2, 3)
    fit <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    expect_true(fit$converged)
    score <- EDI:::get_stereotype_logit_score_cpp(d$X, d$y, fit$params)
    expect_lt(max(abs(score)), 5e-3)
})

# Test 5: K=2 (binary, n_gamma=0, no score_v/cum_v used) works
test_that("K=2 fit converges", {
    d <- make_ordinal_data(200, 2, 2)
    fit <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    expect_true(fit$converged)
})

# Test 6: K=5 (n_gamma=3) works and hessian PSD
test_that("K=5 fit converges and fisher info PSD", {
    d <- make_ordinal_data(300, 2, 5)
    fit <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    expect_true(fit$converged)
    ev <- eigen(fit$fisher_information, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > -1e-3))
})

# Test 7: intercept-only model (p=1 → just thresholds + gammas)
test_that("intercept-only model converges (K=4)", {
    set.seed(7)
    X <- matrix(1, nrow = 300, ncol = 1)
    y <- sample(1:4, 300, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
    fit <- EDI:::fast_stereotype_logit_with_var_cpp(X, y)
    expect_true(fit$converged)
    expect_length(fit$b, 1)
})

# Test 8: estimate_only flag works correctly
test_that("estimate_only returns consistent params", {
    d <- make_ordinal_data(200, 2, 4)
    full <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    est  <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y, estimate_only = TRUE)
    expect_equal(full$params, est$params, tolerance = 1e-8)
})

# Test 9: ssq_b_1 is finite and positive
test_that("ssq_b_1 is finite positive (K=4)", {
    d <- make_ordinal_data(300, 2, 4)
    fit <- EDI:::fast_stereotype_logit_with_var_cpp(d$X, d$y)
    expect_true(is.finite(fit$ssq_b_1))
    expect_gt(fit$ssq_b_1, 0)
})

# Test 10: interleaved calls don't corrupt member state
test_that("interleaved calls with different data don't corrupt each other", {
    d1 <- make_ordinal_data(200, 2, 4, seed = 1)
    d2 <- make_ordinal_data(200, 2, 4, seed = 2)
    ref1 <- EDI:::fast_stereotype_logit_with_var_cpp(d1$X, d1$y)
    ref2 <- EDI:::fast_stereotype_logit_with_var_cpp(d2$X, d2$y)
    # After calling d2, d1 should still give same result
    again1 <- EDI:::fast_stereotype_logit_with_var_cpp(d1$X, d1$y)
    expect_equal(ref1$params, again1$params)
})
