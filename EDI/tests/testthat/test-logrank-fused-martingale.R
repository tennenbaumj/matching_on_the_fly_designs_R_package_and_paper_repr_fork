library(testthat)
library(EDI)

# Verify fast_logrank_stats_cpp after fusing martingale mean/variance accumulation
# into the main sweep (TODO-47). Reference: survival::survdiff (score/var_score)
# and coxph martingale residuals (beta_hat / se_beta_hat).

skip_if_not_installed("survival")

make_surv_data <- function(n, seed) {
    set.seed(seed)
    w    <- as.integer(sample(0:1, n, replace = TRUE))
    y    <- rexp(n, rate = exp(0.3 * w))
    dead <- as.integer(rbinom(n, 1, 0.8))
    list(w = w, y = y, dead = dead)
}

test_that("logrank score and var_score match survival::survdiff", {
    d <- make_surv_data(120L, seed = 1L)
    res <- EDI:::fast_logrank_stats_cpp(d$w, d$y, d$dead)

    df  <- data.frame(time = d$y, status = d$dead, w = factor(d$w))
    sd  <- survival::survdiff(survival::Surv(time, status) ~ w, data = df)
    expect_equal(res$score,     sd$obs[2] - sd$exp[2], tolerance = 1e-10)
    expect_equal(res$var_score, sd$var[2, 2],           tolerance = 1e-10)
})

test_that("logrank beta_hat and se_beta_hat match coxph martingale residuals", {
    d <- make_surv_data(150L, seed = 2L)
    res <- EDI:::fast_logrank_stats_cpp(d$w, d$y, d$dead)

    df  <- data.frame(time = d$y, status = d$dead)
    fit <- survival::coxph(survival::Surv(time, status) ~ 1, data = df)
    m   <- residuals(fit, type = "martingale")

    beta_ref <- mean(m[d$w == 1]) - mean(m[d$w == 0])
    n1 <- sum(d$w == 1); n0 <- sum(d$w == 0)
    se_ref   <- sqrt(var(m[d$w == 1]) / n1 + var(m[d$w == 0]) / n0)

    expect_equal(res$beta_hat,    beta_ref, tolerance = 1e-10)
    expect_equal(res$se_beta_hat, se_ref,   tolerance = 1e-10)
})

test_that("logrank handles all-dead and no-ties data correctly", {
    d <- make_surv_data(80L, seed = 3L)
    d$dead <- rep(1L, 80L)  # all events
    res <- EDI:::fast_logrank_stats_cpp(d$w, d$y, d$dead)
    expect_true(is.finite(res$score))
    expect_true(is.finite(res$beta_hat))
    expect_true(is.finite(res$se_beta_hat))
})

test_that("logrank output matches across multiple seeds", {
    for (seed in c(10L, 20L, 30L, 40L)) {
        d   <- make_surv_data(200L, seed = seed)
        res <- EDI:::fast_logrank_stats_cpp(d$w, d$y, d$dead)

        df  <- data.frame(time = d$y, status = d$dead, w = factor(d$w))
        sd  <- survival::survdiff(survival::Surv(time, status) ~ w, data = df)

        expect_equal(res$score,     sd$obs[2] - sd$exp[2], tolerance = 1e-10,
                     label = paste("score seed", seed))
        expect_equal(res$var_score, sd$var[2, 2],           tolerance = 1e-10,
                     label = paste("var_score seed", seed))
    }
})
