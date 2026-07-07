library(EDI)
library(testthat)

test_that("d_optimal_search_cpp returns matrix with correct dimensions", {
    set.seed(42)
    n <- 30; p <- 4; n_T <- 15; nsim <- 50L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::d_optimal_search_cpp(P, nsim, n_T)
    expect_equal(dim(res), c(n, nsim))
})

test_that("d_optimal_search_cpp each column has exactly n_T ones", {
    set.seed(99)
    n <- 40; p <- 5; n_T <- 20; nsim <- 100L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::d_optimal_search_cpp(P, nsim, n_T)
    col_sums <- colSums(res)
    expect_true(all(col_sums == n_T),
        info = paste("col sums range:", range(col_sums)))
})

test_that("d_optimal_search_cpp values are 0 or 1 only", {
    set.seed(7)
    n <- 20; p <- 3; n_T <- 10; nsim <- 30L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::d_optimal_search_cpp(P, nsim, n_T)
    expect_true(all(res %in% c(0L, 1L)))
})

test_that("d_optimal_search_cpp local optimality: no single swap improves d-criterion", {
    set.seed(123)
    n <- 24; p <- 4; n_T <- 12; nsim <- 20L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::d_optimal_search_cpp(P, nsim, n_T)
    # Check the best column found is a local optimum
    w <- res[, which.min(sapply(seq_len(nsim), function(s) {
        w <- res[, s]
        Xw <- X[w == 1, , drop = FALSE]
        if (nrow(Xw) < p) return(Inf)
        det(t(Xw) %*% Xw)
    }))]
    t_idxs <- which(w == 1)
    c_idxs <- which(w == 0)
    Pw <- P %*% w
    is_local_opt <- TRUE
    for (i in t_idxs) {
        for (j in c_idxs) {
            delta <- -2 * Pw[i] + P[i, i] + 2 * Pw[j] + P[j, j] - 2 * P[i, j]
            if (delta < -1e-9) {
                is_local_opt <- FALSE
                break
            }
        }
        if (!is_local_opt) break
    }
    expect_true(is_local_opt)
})

test_that("d_optimal_search_cpp larger problem", {
    set.seed(55)
    n <- 80; p <- 8; n_T <- 40; nsim <- 30L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::d_optimal_search_cpp(P, nsim, n_T)
    expect_equal(dim(res), c(n, nsim))
    expect_true(all(colSums(res) == n_T))
    expect_true(all(res %in% c(0L, 1L)))
})

# --- a_optimal_search_cpp tests ---

test_that("a_optimal_search_cpp returns matrix with correct dimensions", {
    set.seed(42)
    n <- 30; p <- 4; n_T <- 15; nsim <- 50L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    H <- X %*% solve(t(X) %*% X) %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::a_optimal_search_cpp(P, H, nsim, n_T)
    expect_equal(dim(res), c(n, nsim))
})

test_that("a_optimal_search_cpp each column has exactly n_T ones", {
    set.seed(77)
    n <- 40; p <- 5; n_T <- 20; nsim <- 80L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    H <- X %*% solve(t(X) %*% X) %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::a_optimal_search_cpp(P, H, nsim, n_T)
    expect_true(all(colSums(res) == n_T))
})

test_that("a_optimal_search_cpp values are 0 or 1 only", {
    set.seed(11)
    n <- 20; p <- 3; n_T <- 10; nsim <- 25L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    H <- X %*% solve(t(X) %*% X) %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::a_optimal_search_cpp(P, H, nsim, n_T)
    expect_true(all(res %in% c(0L, 1L)))
})

test_that("a_optimal_search_cpp larger problem", {
    set.seed(88)
    n <- 60; p <- 6; n_T <- 30; nsim <- 40L
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    H <- X %*% solve(t(X) %*% X) %*% solve(t(X) %*% X) %*% t(X)
    res <- EDI:::a_optimal_search_cpp(P, H, nsim, n_T)
    expect_equal(dim(res), c(n, nsim))
    expect_true(all(colSums(res) == n_T))
    expect_true(all(res %in% c(0L, 1L)))
})
