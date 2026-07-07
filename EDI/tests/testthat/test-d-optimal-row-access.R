library(EDI)
library(testthat)

make_hat_matrix <- function(n, p, seed = 42) {
    set.seed(seed)
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    list(P = P, X = X)
}

d_optimal_obj <- function(w, P) {
    # Objective: minimize w'Pw (D-optimal equivalent in projection space)
    sum(w * (P %*% w))
}

# Test 1: output dimensions correct
test_that("d_optimal_search_cpp returns correct dimensions", {
    d <- make_hat_matrix(20, 3)
    res <- EDI:::d_optimal_search_cpp(d$P, 10L, 10L)
    expect_equal(dim(res), c(20L, 10L))
})

# Test 2: each column sums to n_T
test_that("d_optimal columns sum to n_T", {
    d <- make_hat_matrix(20, 3)
    res <- EDI:::d_optimal_search_cpp(d$P, 20L, 10L)
    expect_true(all(colSums(res) == 10L))
})

# Test 3: entries are 0/1
test_that("d_optimal entries are 0 or 1", {
    d <- make_hat_matrix(20, 3)
    res <- EDI:::d_optimal_search_cpp(d$P, 10L, 10L)
    expect_true(all(res %in% c(0L, 1L)))
})

# Test 4: objective value is < random assignment average (optimiser improves)
test_that("d_optimal objective better than random average", {
    set.seed(1)
    n <- 30; n_T <- 15
    d <- make_hat_matrix(n, 4, seed = 1)
    res <- EDI:::d_optimal_search_cpp(d$P, 100L, n_T)
    opt_objs <- apply(res, 2, function(w) d_optimal_obj(w, d$P))
    # Correct random expected wPw: diagonal + off-diagonal terms
    p <- 4
    rand_expected <- p*n_T/n + (sum(d$P) - p)*n_T*(n_T-1)/(n*(n-1))
    expect_lt(mean(opt_objs), rand_expected)
})

# Test 5: deterministic with fixed seed (same P, same nsim) → same result each call
test_that("d_optimal results are deterministic in distribution", {
    d <- make_hat_matrix(20, 3)
    r1 <- EDI:::d_optimal_search_cpp(d$P, 50L, 10L)
    r2 <- EDI:::d_optimal_search_cpp(d$P, 50L, 10L)
    # Both should produce valid 0/1 assignments; objectives should be close
    obj1 <- mean(apply(r1, 2, function(w) d_optimal_obj(w, d$P)))
    obj2 <- mean(apply(r2, 2, function(w) d_optimal_obj(w, d$P)))
    expect_lt(abs(obj1 - obj2) / obj1, 0.1)  # within 10% of each other
})

# Test 6: a_optimal_search_cpp dimensions
test_that("a_optimal_search_cpp returns correct dimensions", {
    d <- make_hat_matrix(20, 3)
    P <- d$P
    H <- P  # use same symmetric matrix for H
    res <- EDI:::a_optimal_search_cpp(P, H, 10L, 10L)
    expect_equal(dim(res), c(20L, 10L))
})

# Test 7: a_optimal column sums to n_T
test_that("a_optimal columns sum to n_T", {
    d <- make_hat_matrix(20, 3)
    P <- d$P; H <- P
    res <- EDI:::a_optimal_search_cpp(P, H, 20L, 10L)
    expect_true(all(colSums(res) == 10L))
})

# Test 8: a_optimal entries are 0/1
test_that("a_optimal entries are 0 or 1", {
    d <- make_hat_matrix(20, 3)
    res <- EDI:::a_optimal_search_cpp(d$P, d$P, 10L, 10L)
    expect_true(all(res %in% c(0L, 1L)))
})

# Test 9: optimised assignments beat the correct random expected value
test_that("d_optimal beats random expected wPw", {
    set.seed(42)
    n <- 10; n_T <- 5; p <- 2
    d <- make_hat_matrix(n, p, seed = 42)
    res <- EDI:::d_optimal_search_cpp(d$P, 200L, n_T)
    objs <- apply(res, 2, function(w) d_optimal_obj(w, d$P))
    rand_base <- p*n_T/n + (sum(d$P) - p)*n_T*(n_T-1)/(n*(n-1))
    expect_lt(mean(objs), rand_base)
})

# Test 10: row-access fix correctness — result matches brute-force for small case
test_that("d_optimal objective matches brute-force for tiny n", {
    set.seed(7)
    n <- 6; n_T <- 3
    X <- cbind(1, rnorm(n))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::d_optimal_search_cpp(P, 500L, n_T)
    # Brute force: enumerate all C(6,3) = 20 assignments
    combs <- combn(n, n_T)
    bf_objs <- apply(combs, 2, function(idx) {
        w <- integer(n); w[idx] <- 1L
        d_optimal_obj(w, P)
    })
    best_bf <- min(bf_objs)
    # At least some of the 500 simulations should find the brute-force optimum
    opt_objs <- apply(res, 2, function(w) d_optimal_obj(w, P))
    expect_lt(min(opt_objs), best_bf + 1e-10)
})
