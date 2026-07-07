stereotype_scores_reference <- function(gamma, K) {
  scores <- numeric(K)
  if (K >= 2L) {
    scores[K] <- 1
  }
  if (length(gamma) > 0L) {
    v <- exp(gamma)
    denom <- 1 + sum(v)
    scores[2:(K - 1L)] <- cumsum(v) / denom
  }
  scores
}

test_that("stereotype logit GEMV gradient matches independent likelihood derivatives", {
  skip_if_not_installed("numDeriv")
  set.seed(9001)
  n <- 160L
  K <- 4L
  X <- matrix(rnorm(n * 3L, sd = 0.5), ncol = 3L)
  y <- rep(seq_len(K), length.out = n)
  params <- c(-0.9, 0.15, 1.1, 0.35, -0.2, 0.1, -0.3, 0.25)

  loglik_r <- function(par) {
    n_alpha <- K - 1L
    p <- ncol(X)
    alpha <- par[seq_len(n_alpha)]
    beta <- par[n_alpha + seq_len(p)]
    gamma <- par[(n_alpha + p + 1L):length(par)]
    scores <- stereotype_scores_reference(gamma, K)
    eta <- drop(X %*% beta)
    logits <- cbind(
      0,
      alpha[1L] + scores[2L] * eta,
      alpha[2L] + scores[3L] * eta,
      alpha[3L] + scores[4L] * eta
    )
    max_logits <- apply(logits, 1L, max)
    shifted <- logits - max_logits
    log_denom <- max_logits + log(rowSums(exp(shifted)))
    sum(logits[cbind(seq_len(n), y)] - log_denom)
  }

  score <- EDI:::get_stereotype_logit_score_cpp(X, y, params)
  hessian <- EDI:::get_stereotype_logit_hessian_cpp(X, y, params)

  expect_equal(
    as.numeric(score),
    as.numeric(numDeriv::grad(loglik_r, params)),
    tolerance = 1e-7
  )
  expect_equal(
    unname(hessian),
    unname(numDeriv::hessian(loglik_r, params)),
    tolerance = 1e-4
  )
})

test_that("stereotype logit mutable gradient and Hessian buffers are repeatable", {
  set.seed(9002)
  n <- 120L
  X <- matrix(rnorm(n * 2L), ncol = 2L)
  y <- rep(1:3, length.out = n)
  params <- c(-0.75, 0.6, 0.25, -0.35, 0.1)

  score1 <- EDI:::get_stereotype_logit_score_cpp(X, y, params)
  score2 <- EDI:::get_stereotype_logit_score_cpp(X, y, params)
  hessian1 <- EDI:::get_stereotype_logit_hessian_cpp(X, y, params)
  hessian2 <- EDI:::get_stereotype_logit_hessian_cpp(X, y, params)

  expect_equal(score1, score2, tolerance = 0)
  expect_equal(hessian1, hessian2, tolerance = 0)
})
