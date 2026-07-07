# Reference R implementation of the compound KK estimator (one simulation column).
kk_compound_ref <- function(y, w, m_vec) {
  max_m <- max(m_vec)

  d_bar    <- NA_real_
  ssqD_bar <- NA_real_
  if (max_m > 0) {
    treated_idx <- integer(max_m)
    control_idx <- integer(max_m)
    for (i in seq_along(y)) {
      mid <- m_vec[i]
      if (mid > 0) {
        if (w[i] == 1L) treated_idx[mid] <- i
        else            control_idx[mid] <- i
      }
    }
    diffs <- y[treated_idx] - y[control_idx]
    d_bar <- mean(diffs)
    if (max_m > 1) ssqD_bar <- var(diffs) / max_m
  }

  unm_T <- y[m_vec == 0 & w == 1L]
  unm_C <- y[m_vec == 0 & w == 0L]
  nRT <- length(unm_T); nRC <- length(unm_C)

  r_bar <- NA_real_
  ssqR  <- NA_real_
  if (nRT > 0 && nRC > 0) {
    r_bar <- mean(unm_T) - mean(unm_C)
    if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
      nR    <- nRT + nRC
      var_T <- var(unm_T)
      var_C <- var(unm_C)
      ssqR  <- (var_T * (nRT - 1) + var_C * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)
    }
  }

  if (nRT <= 1 || nRC <= 1) return(d_bar)
  if (max_m == 0)            return(r_bar)
  if (!is.finite(ssqD_bar) || ssqD_bar <= 0) return(r_bar)
  if (!is.finite(ssqR)     || ssqR     <= 0) return(d_bar)
  w_star <- ssqR / (ssqR + ssqD_bar)
  w_star * d_bar + (1 - w_star) * r_bar
}

test_that("merged-pass matches reference on mixed matched/unmatched data", {
  set.seed(7L)
  n    <- 200L
  nsim <- 500L
  n_pairs <- 60L
  y <- rnorm(n)

  make_col <- function() {
    w <- integer(n); m <- integer(n)
    for (k in seq_len(n_pairs)) {
      w[2L*k-1L] <- 1L; w[2L*k] <- 0L
      m[2L*k-1L] <- k;  m[2L*k] <- k
    }
    idx_u <- (2L*n_pairs+1L):n
    w[idx_u] <- sample(c(0L,1L), length(idx_u), replace=TRUE)
    list(w=w, m=m)
  }
  cols  <- lapply(seq_len(nsim), function(b) make_col())
  w_mat <- do.call(cbind, lapply(cols, `[[`, "w"))
  m_mat <- do.call(cbind, lapply(cols, `[[`, "m"))
  storage.mode(w_mat) <- "integer"
  storage.mode(m_mat) <- "integer"

  cpp_res <- EDI:::compute_matching_compound_distr_parallel_cpp(y, w_mat, m_mat, 1L)
  ref_res <- vapply(seq_len(nsim),
    function(b) kk_compound_ref(y, w_mat[, b], m_mat[, b]),
    numeric(1))

  expect_equal(as.numeric(cpp_res), ref_res, tolerance = 1e-10)
})

test_that("merged-pass handles all-matched (no unmatched obs)", {
  set.seed(13L)
  n <- 100L; n_pairs <- 50L
  y <- rnorm(n)
  w <- rep(c(1L, 0L), n_pairs)
  m <- rep(seq_len(n_pairs), each = 2L)
  w_mat <- matrix(w, ncol = 1L); m_mat <- matrix(m, ncol = 1L)
  storage.mode(w_mat) <- "integer"; storage.mode(m_mat) <- "integer"

  cpp_res <- EDI:::compute_matching_compound_distr_parallel_cpp(y, w_mat, m_mat, 1L)
  ref     <- kk_compound_ref(y, w, m)
  expect_equal(as.numeric(cpp_res), ref, tolerance = 1e-10)
})

test_that("merged-pass handles all-unmatched (no matched pairs)", {
  set.seed(21L)
  n <- 80L
  y <- rnorm(n)
  w <- sample(c(0L, 1L), n, replace = TRUE)
  m <- integer(n)
  w_mat <- matrix(w, ncol = 1L); m_mat <- matrix(m, ncol = 1L)
  storage.mode(w_mat) <- "integer"; storage.mode(m_mat) <- "integer"

  cpp_res <- EDI:::compute_matching_compound_distr_parallel_cpp(y, w_mat, m_mat, 1L)
  ref     <- kk_compound_ref(y, w, m)
  expect_equal(as.numeric(cpp_res), ref, tolerance = 1e-10)
})
