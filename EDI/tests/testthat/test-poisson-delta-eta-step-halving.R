library(testthat)
library(EDI)

# Helper: force step-halving by providing a warm start far from the MLE.
# A warm start that inflates mu exponentially causes the Newton step to overshoot,
# triggering the backtracking loop.
make_halving_scenario <- function(seed = 42, n = 600, p = 4) {
    set.seed(seed)
    X <- cbind(1, matrix(rnorm(n * (p - 1L)), ncol = p - 1L))
    beta_true <- c(0.5, -0.4, 0.3, -0.2)
    y <- stats::rpois(n, exp(drop(X %*% beta_true)))
    list(X = X, y = y, beta_true = beta_true)
}

test_that("Poisson IRLS delta_eta: standard fit matches glm() coefficients", {
    sc <- make_halving_scenario(seed = 1)
    fit <- fast_poisson_regression_cpp(sc$X, sc$y, optimization_alg = "irls")
    ref <- stats::glm.fit(sc$X, sc$y, family = stats::poisson())
    expect_equal(as.numeric(fit$b), as.numeric(ref$coefficients), tolerance = 1e-7)
})

test_that("Poisson IRLS delta_eta: bad warm start (triggers step-halving) matches glm()", {
    sc <- make_halving_scenario(seed = 2)
    # warm start far from MLE forces step-halving on early iterations
    bad_start <- sc$beta_true * 8
    fit <- fast_poisson_regression_cpp(
        sc$X, sc$y,
        warm_start_beta = bad_start,
        optimization_alg = "irls"
    )
    ref <- stats::glm.fit(sc$X, sc$y, family = stats::poisson())
    expect_equal(as.numeric(fit$b), as.numeric(ref$coefficients), tolerance = 1e-6)
    expect_true(fit$converged)
})

test_that("Poisson IRLS delta_eta: weighted fit with bad warm start matches glm()", {
    sc <- make_halving_scenario(seed = 3, n = 400)
    weights <- runif(400, 0.5, 2.0)
    bad_start <- sc$beta_true * 6
    fit <- fast_poisson_regression_weighted_cpp(
        sc$X, sc$y, weights,
        warm_start_beta = bad_start,
        optimization_alg = "irls"
    )
    ref <- stats::glm.fit(sc$X, sc$y, weights = weights, family = stats::poisson())
    expect_equal(as.numeric(fit$b), as.numeric(ref$coefficients), tolerance = 1e-6)
    expect_true(fit$converged)
})

test_that("Poisson IRLS delta_eta: score near zero and information PD at convergence", {
    sc <- make_halving_scenario(seed = 4)
    fit <- fast_poisson_regression_cpp(sc$X, sc$y, optimization_alg = "irls")
    expect_lt(max(abs(fit$score)), 1e-4)
    eigs <- eigen(fit$XtWX, only.values = TRUE)$values
    expect_true(all(eigs > 0))
})

test_that("Poisson IRLS delta_eta: repeated calls are deterministic", {
    sc <- make_halving_scenario(seed = 5)
    bad_start <- sc$beta_true * 5
    fit1 <- fast_poisson_regression_cpp(sc$X, sc$y, warm_start_beta = bad_start,
                                        optimization_alg = "irls")
    fit2 <- fast_poisson_regression_cpp(sc$X, sc$y, warm_start_beta = bad_start,
                                        optimization_alg = "irls")
    expect_identical(fit1$b, fit2$b)
    expect_identical(fit1$converged, fit2$converged)
})

test_that("Poisson IRLS delta_eta: estimate_only path unaffected", {
    sc <- make_halving_scenario(seed = 6)
    bad_start <- sc$beta_true * 7
    fit_full <- fast_poisson_regression_cpp(sc$X, sc$y, warm_start_beta = bad_start,
                                            optimization_alg = "irls")
    fit_est  <- fast_poisson_regression_cpp(sc$X, sc$y, warm_start_beta = bad_start,
                                            optimization_alg = "irls", estimate_only = TRUE)
    expect_equal(as.numeric(fit_est$b), as.numeric(fit_full$b), tolerance = 1e-12)
})
