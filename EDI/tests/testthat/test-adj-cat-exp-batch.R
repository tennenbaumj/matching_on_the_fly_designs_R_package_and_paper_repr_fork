library(testthat)
library(EDI)

# Verify fast_adjacent_category_logit_cpp after batching exp(alpha[k]) precompute
# and switching to product-recurrence softmax (TODO-50).
# Reference: VGAM::vglm(family=acat(parallel=TRUE)) for coefficients/vcov;
# finite differences for Hessian.

skip_if_not_installed("VGAM")

make_adj_data <- function(n, K, seed) {
	set.seed(seed)
	p <- 3L
	X <- matrix(rnorm(n * p), n, p)
	latent <- drop(X %*% rnorm(p, sd = 0.4)) + rlogis(n)
	breaks <- quantile(latent, seq(0, 1, length.out = K + 1L))
	breaks[1] <- -Inf; breaks[K + 1L] <- Inf
	y <- as.numeric(cut(latent, breaks))
	list(X = X, y = y, K = K)
}

# ── coefficients and vcov match VGAM ────────────────────────────────────────

test_that("adj-cat coefficients match VGAM for K=3", {
	d      <- make_adj_data(400L, 3L, seed = 1L)
	res    <- fast_adjacent_category_logit_with_var_cpp(d$X, d$y)
	ref    <- VGAM::vglm(d$y ~ d$X, family = VGAM::acat(parallel = TRUE))
	n_a    <- d$K - 1L; p <- ncol(d$X)
	bidx   <- (n_a + 1L):(n_a + p)          # beta positions in params vector
	b_r    <- as.numeric(stats::coef(ref))[bidx]
	b_c    <- as.numeric(res$b)
	if (sum((b_c - b_r)^2) > sum((b_c + b_r)^2)) b_c <- -b_c
	expect_equal(b_c, b_r, tolerance = 1e-4)
	expect_equal(unname(diag(res$vcov)[bidx]),
	             as.numeric(diag(stats::vcov(ref))[bidx]), tolerance = 1e-3)
})

test_that("adj-cat coefficients match VGAM for K=4", {
	d      <- make_adj_data(500L, 4L, seed = 2L)
	res    <- fast_adjacent_category_logit_with_var_cpp(d$X, d$y)
	ref    <- VGAM::vglm(d$y ~ d$X, family = VGAM::acat(parallel = TRUE))
	n_a    <- d$K - 1L; p <- ncol(d$X)
	bidx   <- (n_a + 1L):(n_a + p)
	b_r    <- as.numeric(stats::coef(ref))[bidx]
	b_c    <- as.numeric(res$b)
	if (sum((b_c - b_r)^2) > sum((b_c + b_r)^2)) b_c <- -b_c
	expect_equal(b_c, b_r, tolerance = 1e-4)
})

# ── score at MLE ≈ 0 ────────────────────────────────────────────────────────

test_that("score at MLE is near zero for K=3", {
	d   <- make_adj_data(300L, 3L, seed = 3L)
	fit <- fast_adjacent_category_logit_cpp(d$X, d$y)
	expect_true(fit$converged)
	sc  <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, fit$params)
	expect_lt(sqrt(sum(sc^2)), 0.1)  # lbfgs converges on function value, not gradient
})

test_that("score at MLE is near zero for K=4", {
	d   <- make_adj_data(400L, 4L, seed = 4L)
	fit <- fast_adjacent_category_logit_cpp(d$X, d$y)
	expect_true(fit$converged)
	sc  <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, fit$params)
	expect_lt(sqrt(sum(sc^2)), 0.1)
})

# ── score norm is smaller at MLE than at perturbed point ────────────────────

test_that("score norm is smaller at MLE than at perturbed point (K=3)", {
	d      <- make_adj_data(200L, 3L, seed = 5L)
	fit    <- fast_adjacent_category_logit_cpp(d$X, d$y)
	p0     <- fit$params
	sc_mle <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, p0)
	set.seed(99L)
	sc_pert <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, p0 + rnorm(length(p0), sd = 0.2))
	expect_lt(sqrt(sum(sc_mle^2)), sqrt(sum(sc_pert^2)))
})

# ── Hessian matches finite-difference second derivative ──────────────────────

test_that("Hessian matches finite differences for K=3", {
	d   <- make_adj_data(300L, 3L, seed = 6L)
	fit <- fast_adjacent_category_logit_cpp(d$X, d$y)
	p0  <- fit$params
	np  <- length(p0)
	H   <- EDI:::get_adjacent_category_logit_hessian_cpp(d$X, d$y, p0)

	h <- 1e-4
	H_fd <- matrix(0, np, np)
	for (j in seq_len(np)) {
		pp <- pm <- p0; pp[j] <- pp[j] + h; pm[j] <- pm[j] - h
		sc_p <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, pp)
		sc_m <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, pm)
		H_fd[, j] <- (sc_p - sc_m) / (2 * h)
	}
	expect_equal(H, H_fd, tolerance = 1e-4)
})

test_that("Hessian matches finite differences for K=4", {
	d   <- make_adj_data(400L, 4L, seed = 7L)
	fit <- fast_adjacent_category_logit_cpp(d$X, d$y)
	p0  <- fit$params
	np  <- length(p0)
	H   <- EDI:::get_adjacent_category_logit_hessian_cpp(d$X, d$y, p0)

	h <- 1e-4
	H_fd <- matrix(0, np, np)
	for (j in seq_len(np)) {
		pp <- pm <- p0; pp[j] <- pp[j] + h; pm[j] <- pm[j] - h
		sc_p <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, pp)
		sc_m <- EDI:::get_adjacent_category_logit_score_cpp(d$X, d$y, pm)
		H_fd[, j] <- (sc_p - sc_m) / (2 * h)
	}
	expect_equal(H, H_fd, tolerance = 1e-4)
})
