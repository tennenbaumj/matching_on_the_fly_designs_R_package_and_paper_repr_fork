test_that("stereotype logit Hessian workspace is symmetric and repeatable", {
  set.seed(5901)
  n <- 180L
  p <- 4L
  K <- 5L
  X <- matrix(rnorm(n * p, sd = 0.6), n, p)
  y <- sample(rep(seq_len(K), length.out = n))
  d <- (K - 1L) + p + (K - 2L)
  params <- rnorm(d, sd = 0.12)
  params[seq_len(K - 1L)] <- seq(-0.6, 0.6, length.out = K - 1L)

  H1 <- EDI:::get_stereotype_logit_hessian_cpp(X, y, params)
  H2 <- EDI:::get_stereotype_logit_hessian_cpp(X, y, params)

  expect_equal(H1, H2, tolerance = 0)
  expect_equal(unname(H1), unname(t(H1)), tolerance = 1e-12)
  expect_true(all(is.finite(H1)))
})

test_that("stereotype logit Hessian workspace is reset across parameter calls", {
  set.seed(5902)
  n <- 140L
  p <- 3L
  K <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- sample(rep(seq_len(K), length.out = n))
  d <- (K - 1L) + p + (K - 2L)
  params1 <- rnorm(d, sd = 0.1)
  params2 <- params1 + seq(-0.05, 0.05, length.out = d)

  H1_before <- EDI:::get_stereotype_logit_hessian_cpp(X, y, params1)
  invisible(EDI:::get_stereotype_logit_hessian_cpp(X, y, params2))
  H1_after <- EDI:::get_stereotype_logit_hessian_cpp(X, y, params1)

  expect_equal(H1_after, H1_before, tolerance = 0)
  expect_false(isTRUE(all.equal(H1_before, EDI:::get_stereotype_logit_hessian_cpp(X, y, params2))))
})
