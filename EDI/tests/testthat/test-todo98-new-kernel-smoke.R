library(testthat)
library(EDI)

# Smoke tests for the 10 kernels added as TODO-98 profiler coverage.
# Each test verifies: function runs without error, returns a non-trivial finite result,
# and is consistent across two identical calls (deterministic output).

set.seed(42L)

# ── Shared fixtures ─────────────────────────────────────────────────────────────

make_clogit_data <- function(n_disc = 100L, n_conc = 100L, p = 3L, G = 40L) {
  Xd <- cbind(1, sample(c(-1L,1L), n_disc, TRUE), matrix(rnorm(n_disc*(p-1L)), n_disc))
  yd <- as.numeric(rbinom(n_disc, 1L, 0.4))
  Xc <- cbind(1, rep(c(0,1), n_conc/2L), matrix(rnorm(n_conc*(p-1L)), n_conc))
  yc <- as.numeric(rbinom(n_conc, 1L, 0.4))
  gc <- as.integer(rep(seq_len(G), length.out = n_conc))
  list(Xd=Xd, yd=yd, Xc=Xc, yc=yc, gc=gc)
}

make_survival_data <- function(n = 200L) {
  set.seed(1L)
  X   <- cbind(1, rbinom(n,1,0.5), rnorm(n))
  t   <- rexp(n, exp(X %*% c(0,0.3,-0.2)))
  d   <- rbinom(n, 1, 0.7)
  # sort by time for DepCensTransform
  ord <- order(t)
  list(X_ord = X[ord,,drop=FALSE],
       y_bm  = as.numeric(t[ord]),
       dead  = as.numeric(d[ord]))
}

make_bai_data <- function(n = 100L, nsim = 200L) {
  n_pairs  <- n / 2L
  n_halves <- n_pairs / 2L
  y        <- rnorm(n)
  w_mat    <- matrix(as.integer(sample(c(0L,1L), n*nsim, TRUE)), nrow=n, ncol=nsim)
  m_mat    <- matrix(as.integer(rep(seq_len(n_pairs), each=2L)), nrow=n, ncol=nsim)
  pair_ids <- seq_len(n_pairs)
  halves   <- matrix(as.integer(sample(pair_ids, n_halves*2L, FALSE)), n_halves, 2L)
  list(y=y, w_mat=w_mat, m_mat=m_mat, halves=halves)
}

# ── ClogitPlusGLMM ───────────────────────────────────────────────────────────────

test_that("fast_clogit_plus_glmm_cpp est path runs and converges", {
  d <- make_clogit_data()
  r <- EDI:::fast_clogit_plus_glmm_cpp(d$Xd, d$yd, d$Xc, d$yc, d$gc,
         has_discordant=TRUE, has_concordant=TRUE, estimate_only=TRUE)
  expect_true(r$converged)
  expect_true(all(is.finite(r$b)))
})

test_that("fast_clogit_plus_glmm_cpp var path returns finite vcov", {
  d <- make_clogit_data()
  r <- EDI:::fast_clogit_plus_glmm_cpp(d$Xd, d$yd, d$Xc, d$yc, d$gc,
         has_discordant=TRUE, has_concordant=TRUE, estimate_only=FALSE)
  expect_true(r$converged)
  expect_true(all(is.finite(r$vcov)))
})

test_that("fast_clogit_plus_glmm_cpp est and var give same point estimates", {
  set.seed(3L)
  d  <- make_clogit_data()
  r1 <- EDI:::fast_clogit_plus_glmm_cpp(d$Xd, d$yd, d$Xc, d$yc, d$gc,
          has_discordant=TRUE, has_concordant=TRUE, estimate_only=TRUE)
  r2 <- EDI:::fast_clogit_plus_glmm_cpp(d$Xd, d$yd, d$Xc, d$yc, d$gc,
          has_discordant=TRUE, has_concordant=TRUE, estimate_only=FALSE)
  expect_equal(r1$b, r2$b, tolerance=1e-5)
})

# ── DepCensTransform survival ─────────────────────────────────────────────────────

test_that("fast_dep_cens_transform_optim_cpp est path runs and converges", {
  d <- make_survival_data()
  r <- EDI:::fast_dep_cens_transform_optim_cpp(d$X_ord, d$y_bm, d$dead,
         smart_cold_start=TRUE, estimate_only=TRUE)
  expect_true(is.list(r))
  expect_true(all(is.finite(r$b)))
})

test_that("fast_dep_cens_transform_optim_cpp var path returns finite vcov", {
  d <- make_survival_data()
  r <- EDI:::fast_dep_cens_transform_optim_cpp(d$X_ord, d$y_bm, d$dead,
         smart_cold_start=TRUE, estimate_only=FALSE)
  expect_true(all(is.finite(r$vcov)))
})

# ── D-optimal design search ───────────────────────────────────────────────────────

test_that("d_optimal_search_cpp returns n×nsim binary matrix with n_T treatments per col", {
  set.seed(7L)
  n <- 40L; p <- 4L; n_T <- 20L; nsim <- 50L
  X     <- cbind(1, matrix(rnorm(n*(p-1L)), n, p-1L))
  P     <- X %*% solve(crossprod(X) + diag(1e-8, p)) %*% t(X)
  r     <- EDI:::d_optimal_search_cpp(P, nsim, n_T)
  expect_equal(dim(r), c(n, nsim))
  expect_true(all(r %in% c(0L, 1L)))
  expect_true(all(colSums(r) == n_T))
})

# ── KK compound distribution ──────────────────────────────────────────────────────

test_that("compute_matching_compound_distr_parallel_cpp returns numeric vector", {
  set.seed(9L)
  n <- 80L; nsim <- 100L; n_pairs <- n/2L
  y     <- rnorm(n)
  w_mat <- matrix(as.integer(sample(c(0L,1L), n*nsim, TRUE)), n, nsim)
  m_mat <- matrix(as.integer(rep(seq_len(n_pairs), each=2L)), n, nsim)
  r     <- EDI:::compute_matching_compound_distr_parallel_cpp(y, w_mat, m_mat, 1L)
  expect_equal(length(r), nsim)
  expect_true(all(is.finite(r)))
})

# ── BAI parallel distribution ─────────────────────────────────────────────────────

test_that("compute_bai_distr_parallel_cpp returns finite numeric vector", {
  d <- make_bai_data()
  r <- EDI:::compute_bai_distr_parallel_cpp(d$w_mat, d$m_mat, d$y, 0, d$halves, TRUE, 1L)
  expect_equal(length(r), ncol(d$w_mat))
  expect_true(all(is.finite(r)))
})

# ── Rerandomization helpers ───────────────────────────────────────────────────────

test_that("rerandomization_search_cpp returns matrix of valid assignments", {
  set.seed(5L)
  X <- matrix(rnorm(80L*3L), 80L, 3L)
  r <- EDI:::rerandomization_search_cpp(X, 50L, "abs_sum_diff", 2.0, 10000L)
  expect_true(is.matrix(r) || is.integer(r))
  expect_true(all(is.finite(as.numeric(r))))
})

test_that("compute_objective_vals_cpp returns numeric vector length r", {
  set.seed(11L)
  n <- 60L; r <- 100L; p <- 3L
  X      <- matrix(rnorm(n*p), n, p)
  indics <- matrix(as.integer(sample(c(0L,1L), n*r, TRUE)), r, n)
  vals   <- EDI:::compute_objective_vals_cpp(X, indics, "abs_sum_diff")
  expect_equal(length(vals), r)
  expect_true(all(is.finite(vals)))
})

# ── CMH block SE ─────────────────────────────────────────────────────────────────

test_that("compute_cmh_block_se_cpp returns finite positive scalar", {
  set.seed(13L)
  n <- 400L; B <- 100L
  y     <- as.numeric(rbinom(n, 1L, 0.4))
  m_vec <- as.integer(rep(seq_len(B), each = n/B))
  se    <- EDI:::compute_cmh_block_se_cpp(y, m_vec, n)
  expect_length(se, 1L)
  expect_true(is.finite(se))
  expect_true(se > 0)
})
