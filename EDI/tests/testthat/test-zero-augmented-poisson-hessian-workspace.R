zap_workspace_log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

make_zap_workspace_data <- function(is_hurdle) {
  set.seed(6000 + as.integer(is_hurdle))
  n <- 90L
  p <- 3L
  X <- cbind(1, matrix(rnorm(n * (p - 1L), sd = 0.7), n, p - 1L))
  Xzi <- cbind(1, matrix(rnorm(n * (p - 1L), sd = 0.6), n, p - 1L))
  beta <- c(0.2, 0.25, -0.15)
  gamma <- c(-0.35, -0.2, 0.1)
  lambda <- exp(drop(X %*% beta))
  pi <- plogis(drop(Xzi %*% gamma))
  y <- integer(n)
  for (i in seq_len(n)) {
    if (is_hurdle) {
      if (runif(1) < pi[i]) {
        y[i] <- 0L
      } else {
        repeat {
          yi <- rpois(1, lambda[i])
          if (yi > 0L) {
            y[i] <- yi
            break
          }
        }
      }
    } else {
      y[i] <- if (runif(1) < pi[i]) 0L else rpois(1, lambda[i])
    }
  }
  list(X = X, Xzi = Xzi, y = as.numeric(y), params = c(beta, gamma))
}

zap_workspace_neg_loglik <- function(params, X, y, Xzi, is_hurdle) {
  p_cond <- ncol(X)
  beta <- params[seq_len(p_cond)]
  gamma <- params[-seq_len(p_cond)]
  eta_c <- drop(X %*% beta)
  eta_z <- drop(Xzi %*% gamma)
  lambda <- exp(eta_c)
  pi <- plogis(eta_z)
  lse <- zap_workspace_log1pexp(eta_z)
  y0 <- y == 0
  ll <- numeric(length(y))

  if (is_hurdle) {
    ll[y0] <- log(pmax(pi[y0], 1e-15))
    if (any(!y0)) {
      exp_ml <- exp(-lambda[!y0])
      log1m_exp_ml <- ifelse(
        lambda[!y0] > log(2),
        log1p(-exp_ml),
        log(-expm1(-lambda[!y0]))
      )
      ll[!y0] <- -lse[!y0] + y[!y0] * eta_c[!y0] - lambda[!y0] - log1m_exp_ml
    }
  } else {
    exp_ml <- exp(-lambda)
    ll[y0] <- log(pmax(pi[y0] + (1 - pi[y0]) * exp_ml[y0], 1e-15))
    ll[!y0] <- -lse[!y0] + y[!y0] * eta_c[!y0] - lambda[!y0]
  }

  -sum(ll)
}

test_that("zero-augmented Poisson Hessian workspace matches independent finite differences", {
  for (is_hurdle in c(FALSE, TRUE)) {
    d <- make_zap_workspace_data(is_hurdle)
    hessian_cpp <- EDI:::get_zero_augmented_poisson_hessian_cpp(
      d$X, d$y, d$Xzi, d$params, is_hurdle
    )
    hessian_ref <- -numDeriv::hessian(
      zap_workspace_neg_loglik,
      d$params,
      X = d$X,
      y = d$y,
      Xzi = d$Xzi,
      is_hurdle = is_hurdle
    )

    expect_equal(unname(hessian_cpp), unname(hessian_ref), tolerance = 1e-4)
    expect_equal(unname(hessian_cpp), unname(t(hessian_cpp)), tolerance = 1e-12)
  }
})

test_that("zero-augmented Poisson variance path recomputes Hessian workspace after optimizer use", {
  for (is_hurdle in c(FALSE, TRUE)) {
    d <- make_zap_workspace_data(is_hurdle)
    fit <- EDI:::fast_zero_augmented_poisson_cpp(
      d$X,
      d$y,
      d$Xzi,
      is_hurdle = is_hurdle,
      warm_start_params = d$params,
      smart_cold_start = FALSE,
      estimate_only = FALSE,
      maxit = 200L,
      tol = 1e-8
    )
    observed_information <- -EDI:::get_zero_augmented_poisson_hessian_cpp(
      d$X, d$y, d$Xzi, fit$params, is_hurdle
    )

    expect_true(fit$converged)
    expect_equal(unname(fit$observed_information), unname(observed_information), tolerance = 0)
  }
})
