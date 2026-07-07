bai_reference <- function(w_mat, m_mat, y, delta, halves_idx, convex_flag) {
  nsim <- ncol(w_mat)
  out <- numeric(nsim)
  for (b in seq_len(nsim)) {
    w_col <- w_mat[, b]
    m_col <- m_mat[, b]
    max_match <- max(m_col)
    d_i <- numeric()
    y_r <- numeric()
    w_r <- integer()
    match_T <- match_C <- numeric(max_match)
    has_T <- has_C <- logical(max_match)

    if (max_match > 0L) {
      for (i in seq_along(y)) {
        y_val <- y[i] + if (w_col[i] == 1L) delta else 0
        m <- m_col[i]
        if (m > 0L) {
          if (w_col[i] == 1L) {
            match_T[m] <- y_val
            has_T[m] <- TRUE
          } else {
            match_C[m] <- y_val
            has_C[m] <- TRUE
          }
        } else {
          y_r <- c(y_r, y_val)
          w_r <- c(w_r, w_col[i])
        }
      }
      for (m in seq_len(max_match)) {
        if (has_T[m] && has_C[m]) {
          d_i <- c(d_i, match_T[m] - match_C[m])
        }
      }
    } else {
      y_r <- y + ifelse(w_col == 1L, delta, 0)
      w_r <- w_col
    }

    m_size <- length(d_i)
    if (m_size == 0L && length(y_r) == 0L) {
      out[b] <- NA_real_
      next
    }

    d_bar <- if (m_size > 0L) mean(d_i) else 0
    r_bar <- 0
    ssqR <- 0
    nRT <- nRC <- 0L
    if (length(y_r) > 0L) {
      yT <- y_r[w_r == 1L]
      yC <- y_r[w_r != 1L]
      nRT <- length(yT)
      nRC <- length(yC)
      if (nRT > 0L && nRC > 0L) {
        r_bar <- mean(yT) - mean(yC)
        if (nRT > 1L && nRC > 1L) {
          ssqR <- stats::var(yT) / nRT + stats::var(yC) / nRC
        }
      }
    }

    bai_var_d_bar <- 0
    if (m_size > 0L) {
      tau_sq <- mean(d_i^2)
      lambda_squ <- 0
      if (nrow(halves_idx) > 0L) {
        for (i in seq_len(nrow(halves_idx))) {
          id1 <- halves_idx[i, 1]
          id2 <- halves_idx[i, 2]
          if (
            id1 > 0L && id1 <= max_match &&
              id2 > 0L && id2 <= max_match &&
              has_T[id1] && has_C[id1] && has_T[id2] && has_C[id2]
          ) {
            lambda_squ <- lambda_squ +
              (match_T[id1] - match_C[id1]) * (match_T[id2] - match_C[id2])
          }
        }
        lambda_squ <- lambda_squ / nrow(halves_idx)
      }
      bai_var_d_bar <- max(1e-8, tau_sq - (lambda_squ + d_bar^2) / 2) / m_size
    }

    if (convex_flag && nRT > 1L && nRC > 1L && m_size > 0L && ssqR > 0) {
      w_star <- ssqR / (ssqR + bai_var_d_bar)
      out[b] <- w_star * d_bar + (1 - w_star) * r_bar
    } else if (m_size > 0L) {
      out[b] <- d_bar
    } else {
      out[b] <- r_bar
    }
  }
  out
}

make_bai_workspace_data <- function(nsim = 80L) {
  set.seed(10201)
  y <- c(1.2, -0.4, 0.7, 1.8, -1.1, 0.3, 2.1, -0.8, 0.5, -1.7, 1.4, -0.2)
  n <- length(y)
  base_m <- c(1L, 1L, 2L, 2L, 3L, 3L, 0L, 0L, 0L, 0L, 0L, 0L)
  w_mat <- matrix(0L, n, nsim)
  m_mat <- matrix(base_m, n, nsim)
  for (b in seq_len(nsim)) {
    w <- integer(n)
    for (pair in 1:3) {
      idx <- which(base_m == pair)
      treated_first <- ((b + pair) %% 2L) == 0L
      w[idx] <- if (treated_first) c(1L, 0L) else c(0L, 1L)
    }
    reservoir <- which(base_m == 0L)
    w[reservoir] <- as.integer(((seq_along(reservoir) + b) %% 3L) == 0L)
    if (b %% 5L == 0L) {
      m_mat[, b] <- 0L
    }
    w_mat[, b] <- w
  }
  halves <- matrix(as.integer(c(1L, 2L, 2L, 3L)), ncol = 2L, byrow = TRUE)
  list(y = y, w_mat = w_mat, m_mat = m_mat, halves = halves)
}

test_that("BAI parallel distribution reuses thread-local workspaces without changing results", {
  d <- make_bai_workspace_data()
  expected_convex <- bai_reference(d$w_mat, d$m_mat, d$y, 0.25, d$halves, TRUE)
  expected_plain <- bai_reference(d$w_mat, d$m_mat, d$y, -0.1, d$halves, FALSE)

  actual_convex_1 <- EDI:::compute_bai_distr_parallel_cpp(
    d$w_mat, d$m_mat, d$y, 0.25, d$halves, TRUE, 1L
  )
  actual_convex_2 <- EDI:::compute_bai_distr_parallel_cpp(
    d$w_mat, d$m_mat, d$y, 0.25, d$halves, TRUE, 2L
  )
  actual_plain <- EDI:::compute_bai_distr_parallel_cpp(
    d$w_mat, d$m_mat, d$y, -0.1, d$halves, FALSE, 2L
  )

  expect_equal(actual_convex_1, expected_convex, tolerance = 1e-12)
  expect_equal(actual_convex_2, expected_convex, tolerance = 1e-12)
  expect_equal(actual_plain, expected_plain, tolerance = 1e-12)
  expect_equal(actual_convex_2, actual_convex_1, tolerance = 0)
})
