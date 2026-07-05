# Extracted from test-weibull-frailty.R:60

# test -------------------------------------------------------------------------
set.seed(10)
n_pairs <- 15
n <- 2 * n_pairs
group_id <- rep(seq_len(n_pairs), each = 2)
w <- rep(c(1, 0), n_pairs)
x1 <- rnorm(n)
X <- cbind(w = w, x1 = x1)
u_g <- rnorm(n_pairs, sd = 0.5)[group_id]
y <- rexp(n) * exp(0.6 * w - 0.2 * x1 + u_g)
dead <- rbinom(n, 1, 0.8)
params <- c(0.5, -0.3, 0.1, -0.5)
score <- as.numeric(get_weibull_frailty_score_cpp(X, y, dead, group_id, params))
