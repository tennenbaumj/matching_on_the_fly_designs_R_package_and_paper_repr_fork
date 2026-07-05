# Extracted from test-weibull-frailty.R:85

# test -------------------------------------------------------------------------
set.seed(11)
n_pairs <- 20
n <- 2 * n_pairs
group_id <- rep(seq_len(n_pairs), each = 2)
w <- rep(c(1, 0), n_pairs)
x1 <- rnorm(n)
y <- (rexp(n) * exp(0.5 * w + 0.3 * x1))^1.7
dead <- rbinom(n, 1, 0.85)
fit_r <- survival::survreg(survival::Surv(y, dead) ~ w + x1, dist = "weibull")
X <- cbind(`(Intercept)` = 1, w = w, x1 = x1)
params <- c(as.numeric(stats::coef(fit_r)), log(fit_r$scale), -8)
neg_ll_frailty <- get_weibull_frailty_neg_loglik_cpp(X, y, dead, group_id, params)
