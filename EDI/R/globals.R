# Internal environment for the EDI package to store global state
edi_env = new.env(parent = emptyenv())

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Closure to encapsulate the internal assertion override flag (used by SimulationFramework)
.assert_manager = (function() {
  internal_run_asserts = TRUE
  list(
    toggle = function(on = TRUE) {
      internal_run_asserts <<- isTRUE(on)
      invisible(internal_run_asserts)
    },
    should_run = function() {
      internal_run_asserts && isTRUE(getOption("edi.run_asserts", TRUE))
    }
  )
})()
#' Toggle the execution of assertions throughout the package
#' 
#' @description This function enables or disables the internal input validation checks (assertions)
#' by setting the \code{options(edi.run_asserts = ...)} value.
#' Disabling assertions can provide a significant performance boost in heavy
#' simulations (often 10x-20x speedup), but it removes the safety rails that
#' prevent invalid data from reaching the internal algorithms.
#' 
#' \strong{Warning:} If assertions are disabled, passing malformed or invalid
#' data to package functions may result in cryptic R errors, incorrect
#' statistical results, or even hard system crashes (SEGFAULTs) at the C++ layer.
#' Only disable assertions if you are certain your data is pre-validated and
#' follows the package requirements exactly.
#' 
#' @param on Logical scalar. If TRUE (default), assertions are executed. If FALSE, they are skipped.
#' @keywords internal
#' @export
toggle_asserts = function(on = TRUE) {
  options(edi.run_asserts = isTRUE(on))
  invisible(isTRUE(on))
}
# private method
should_run_asserts = .assert_manager$should_run

stop_bayesian_bootstrap_for_ivwc = function(self_obj = NULL) {
  # Relaxing this to allow Jackknife if implemented.
  # For now, if we reach here, it means it is not yet implemented for this specific IVWC operation.
  cls = if (is.null(self_obj)) "This IVWC class" else class(self_obj)[1]
  stop(
    cls,
    " does not support this Bayesian bootstrap operation. Implementation is missing for this IVWC path."
  )
}

stop_bayesian_bootstrap_for_bai = function(self_obj = NULL) {
  cls = if (is.null(self_obj)) "This Bai-adjusted KK class" else class(self_obj)[1]
  stop(
    cls,
    " does not support this Bayesian bootstrap operation. Implementation is missing for this Bai path."
  )
}

# Bayesian bootstrap helpers were removed to allow Jackknife support across all paths.

weighted_ordinal_bootstrap_surrogate_fit = function(X, y, row_weights, method = c("logistic", "probit", "cauchit", "cloglog"), warm_start_params = NULL) {
  method = match.arg(method)
  X = as.matrix(X)
  y_num = as.integer(y)
  row_weights = as.numeric(row_weights)
  ok = is.finite(row_weights) & row_weights > 0 & is.finite(y_num)
  if (!any(ok)) return(NULL)
  X_fit = X[ok, , drop = FALSE]
  y_fit = y_num[ok]
  w_fit = row_weights[ok]
  if (is.null(colnames(X_fit))) {
    colnames(X_fit) = paste0("x", seq_len(ncol(X_fit)))
  }
  if (!("treatment" %in% colnames(X_fit)) && ncol(X_fit) >= 1L) {
    colnames(X_fit)[1L] = "treatment"
  }
  dat = as.data.frame(X_fit, check.names = FALSE)
  y_levels = sort(unique(y_fit))
  dat$y_ord = ordered(y_fit, levels = y_levels)
  fit_polr = function(start = NULL) {
    tryCatch(
      suppressWarnings(
        MASS::polr(
          y_ord ~ .,
          data = dat,
          weights = w_fit,
          method = method,
          Hess = FALSE,
          start = start
        )
      ),
      error = function(e) NULL
    )
  }
  polr_start = NULL
  warm_start_params = as.numeric(warm_start_params)
  n_coef = ncol(X_fit)
  n_zeta = length(y_levels) - 1L
  if (length(warm_start_params) == n_coef + n_zeta && all(is.finite(warm_start_params))) {
    beta_start = utils::tail(warm_start_params, n_coef)
    zeta_start = utils::head(warm_start_params, n_zeta)
    if (n_zeta == 0L || all(diff(zeta_start) > 0)) {
      polr_start = c(beta_start, zeta_start)
      names(polr_start) = c(colnames(X_fit), paste0("zeta", seq_len(n_zeta)))
    }
  }
  fit = fit_polr(polr_start)
  if (is.null(fit) && !is.null(polr_start)) {
    fit = fit_polr(NULL)
  }
  coef_vec = tryCatch(stats::coef(fit), error = function(e) NULL)
  beta_hat = if (!is.null(coef_vec) && ("treatment" %in% names(coef_vec))) as.numeric(coef_vec[["treatment"]]) else NA_real_
  if (is.finite(beta_hat)) {
    return(list(beta_hat = beta_hat, coefficients = coef_vec, fit_type = paste0("polr_", method)))
  }
  X_lm = cbind(`(Intercept)` = 1, X_fit)
  lm_fit = tryCatch(
    stats::lm.wfit(x = X_lm, y = as.numeric(y_fit), w = w_fit),
    error = function(e) NULL
  )
  coef_lm = if (is.null(lm_fit)) NULL else as.numeric(lm_fit$coefficients)
  names(coef_lm) = if (is.null(coef_lm)) NULL else colnames(X_lm)
  beta_hat_lm = if (!is.null(coef_lm) && ("treatment" %in% names(coef_lm))) as.numeric(coef_lm[["treatment"]]) else NA_real_
  if (!is.finite(beta_hat_lm)) return(NULL)
  list(beta_hat = beta_hat_lm, coefficients = coef_lm, fit_type = "weighted_lm_surrogate")
}

weights_are_effectively_constant = function(weights, tol = sqrt(.Machine$double.eps)) {
  weights = as.numeric(weights)
  ok = is.finite(weights)
  if (!any(ok)) return(FALSE)
  weights = weights[ok]
  (max(weights) - min(weights)) <= tol
}

kk_pair_and_reservoir_bootstrap_weights = function(private_env, row_weights) {
  row_weights = as.numeric(row_weights)
  m_vec = private_env$m
  n = length(row_weights)
  if (is.null(m_vec)) {
    m_vec = rep(NA_integer_, n)
  } else {
    m_vec = as.integer(m_vec)
    if (length(m_vec) != n) {
      m_vec = rep_len(m_vec, n)
    }
  }
  pair_ids = sort(unique(m_vec[is.finite(m_vec) & m_vec > 0L]))
  pair_weights = if (length(pair_ids) > 0L) {
    vapply(pair_ids, function(pid) {
      idx = which(m_vec == pid)
      mean(row_weights[idx], na.rm = TRUE)
    }, numeric(1))
  } else {
    numeric(0)
  }
  reservoir_idx = which(is.na(m_vec) | m_vec <= 0L)
  list(
    pair_ids = pair_ids,
    pair_weights = as.numeric(pair_weights),
    reservoir_idx = reservoir_idx,
    reservoir_weights = as.numeric(row_weights[reservoir_idx])
  )
}

weighted_cox_bootstrap_surrogate_fit = function(time, dead, X, row_weights, strata = NULL, cluster = NULL, warm_start_beta = NULL) {
  X = as.matrix(X)
  row_weights = as.numeric(row_weights)
  ok = is.finite(time) & is.finite(dead) & is.finite(row_weights) & row_weights > 0
  if (!is.null(strata)) ok = ok & !is.na(strata)
  if (!is.null(cluster)) ok = ok & !is.na(cluster)
  if (!any(ok)) return(NULL)
  X_fit = X[ok, , drop = FALSE]
  if (is.null(colnames(X_fit))) colnames(X_fit) = paste0("x", seq_len(ncol(X_fit)))
  if (!("treatment" %in% colnames(X_fit)) && ncol(X_fit) >= 1L) colnames(X_fit)[1L] = "treatment"
  dat = as.data.frame(X_fit, check.names = FALSE)
  dat$.time__ = as.numeric(time[ok])
  dat$.dead__ = as.numeric(dead[ok])
  dat$.wgt__ = as.numeric(row_weights[ok])
  rhs_terms = colnames(X_fit)
  if (!is.null(strata)) {
    dat$.strata__ = factor(as.integer(strata[ok]))
    rhs_terms = c(rhs_terms, "strata(.strata__)")
  }
  if (!is.null(cluster)) {
    dat$.cluster__ = factor(as.integer(cluster[ok]))
    rhs_terms = c(rhs_terms, "cluster(.cluster__)")
  }
  warm_start_beta = as.numeric(warm_start_beta)
  cox_init = if (length(warm_start_beta) == ncol(X_fit) && all(is.finite(warm_start_beta))) {
    warm_start_beta
  } else {
    NULL
  }
  fit_cox = function(init = NULL) {
    tryCatch(
      suppressWarnings(
        survival::coxph(
          stats::as.formula(paste0("survival::Surv(.time__, .dead__) ~ ", paste(rhs_terms, collapse = " + "))),
          data = dat,
          weights = .wgt__,
          init = init
        )
      ),
      error = function(e) NULL
    )
  }
  fit = tryCatch(
    fit_cox(cox_init),
    error = function(e) NULL
  )
  if (is.null(fit) && !is.null(cox_init)) {
    fit = fit_cox(NULL)
  }
  coef_vec = tryCatch(stats::coef(fit), error = function(e) NULL)
  beta_hat = if (!is.null(coef_vec) && ("treatment" %in% names(coef_vec))) as.numeric(coef_vec[["treatment"]]) else NA_real_
  if (!is.finite(beta_hat)) return(NULL)
  list(beta_hat = beta_hat, coefficients = coef_vec, fit = fit)
}

weighted_weibull_bootstrap_surrogate_fit = function(time, dead, X, row_weights, cluster = NULL, warm_start_params = NULL) {
  X = as.matrix(X)
  row_weights = as.numeric(row_weights)
  ok = is.finite(time) & is.finite(dead) & is.finite(row_weights) & row_weights > 0
  if (!is.null(cluster)) ok = ok & !is.na(cluster)
  if (!any(ok)) return(NULL)
  X_fit = X[ok, , drop = FALSE]
  if (is.null(colnames(X_fit))) colnames(X_fit) = paste0("x", seq_len(ncol(X_fit)))
  if (!("treatment" %in% colnames(X_fit)) && ncol(X_fit) >= 1L) colnames(X_fit)[1L] = "treatment"
  dat = as.data.frame(X_fit, check.names = FALSE)
  dat$.time__ = as.numeric(time[ok])
  dat$.dead__ = as.numeric(dead[ok])
  dat$.wgt__ = as.numeric(row_weights[ok])
  rhs_terms = colnames(X_fit)
  if (!is.null(cluster)) {
    dat$.cluster__ = factor(as.integer(cluster[ok]))
  }
  warm_start_params = as.numeric(warm_start_params)
  survreg_init = NULL
  if (length(warm_start_params) >= ncol(X_fit) + 1L && all(is.finite(warm_start_params))) {
    survreg_init = utils::head(warm_start_params, ncol(X_fit) + 1L)
    names(survreg_init) = c("(Intercept)", colnames(X_fit))
  }
  fit_survreg = function(init = NULL) {
    tryCatch(
      suppressWarnings(
        survival::survreg(
          stats::as.formula(paste0("survival::Surv(.time__, .dead__) ~ ", paste(rhs_terms, collapse = " + "))),
          data = dat,
          weights = .wgt__,
          dist = "weibull",
          init = init
        )
      ),
      error = function(e) NULL
    )
  }
  fit = fit_survreg(survreg_init)
  if (is.null(fit) && !is.null(survreg_init)) {
    fit = fit_survreg(NULL)
  }
  coef_vec = tryCatch(stats::coef(fit), error = function(e) NULL)
  beta_hat = if (!is.null(coef_vec) && ("treatment" %in% names(coef_vec))) as.numeric(coef_vec[["treatment"]]) else NA_real_
  if (!is.finite(beta_hat)) return(NULL)
  list(beta_hat = beta_hat, coefficients = coef_vec, fit = fit)
}

# Creates a fork cluster and caps OMP/BLAS threads on each worker to 1.
# Returns the cluster without storing it — callers decide where it lives.
make_configured_fork_cluster = function(n_cores) {
  default_port = tryCatch(parallel:::getClusterOption("port"), error = function(e) NA_integer_)
  candidate_ports = unique(as.integer(c(default_port, sample.int(20000L, 20L) + 10000L)))
  candidate_ports = candidate_ports[is.finite(candidate_ports) & candidate_ports > 0L & candidate_ports <= 65535L]
  last_error = NULL
  cl = NULL
  for (port in candidate_ports) {
    cl = tryCatch(
      parallel::makeForkCluster(n_cores, port = port),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(cl)) break
  }
  if (is.null(cl)) {
    stop(
      "Could not create fork cluster after trying ",
      length(candidate_ports),
      " ports",
      if (!is.null(last_error)) paste0(": ", conditionMessage(last_error)) else ".",
      call. = FALSE
    )
  }
  tryCatch(
    parallel::clusterCall(cl, function() {
      Sys.setenv(
        OMP_NUM_THREADS        = 1L,
        MKL_NUM_THREADS        = 1L,
        OPENBLAS_NUM_THREADS   = 1L,
        GOTO_NUM_THREADS       = 1L,
        VECLIB_MAXIMUM_THREADS = 1L,
        NUMEXPR_NUM_THREADS    = 1L
      )
      options(mc.cores = 1L)
      if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
      if (requireNamespace("fixest", quietly = TRUE)) suppressWarnings(try(fixest::setFixest_nthreads(1L), silent = TRUE))
      invisible(NULL)
    }),
    error = function(e) invisible(NULL)
  )
  cl
}
#' Set the number of cores for parallelization
#'
#' This function initializes a persistent parallel cluster (either a fork cluster
#' on Unix-like systems or a mirai cluster on others) to be used by all Design
#' and Inference objects. This avoids the overhead of creating clusters
#' repeatedly.
#'
#' @details
#' \code{set_num_cores()} sets a global upper bound for parallel work. It does
#' not guarantee that every inference routine will use all requested workers.
#' EDI's inference dispatcher applies a blocklist-first heuristic informed by
#' package benchmarks: workloads that have shown consistent multicore slowdowns
#' are forced to run serially, while the remaining workloads are allowed to use
#' their method-specific warmup heuristics and native thread caps.
#'
#' The default forced-serial blocklist covers incidence randomization confidence
#' intervals, bootstrap for non-regression KK Wilcoxon inference, bootstrap for
#' non-KK survival procedures, and bootstrap for incidence procedures. Do not
#' expect a universal "more cores is faster" rule.
#'
#' If you want to change the default policy, use
#' \code{set_parallel_dispatch_policy()}.
#' 
#' @param num_cores Integer number of worker processes to make available.
#' @param force_mirai If \code{TRUE}, forces the use of the \code{mirai} package
#'   even on systems where forking is available.
#' 
#' @return Invisible \code{NULL}.
#' 
#' @examples
#' set_num_cores(2)
#' unset_num_cores()
#' @export
set_num_cores = function(num_cores, force_mirai = FALSE) {
  checkmate::assertCount(num_cores, positive = TRUE)
  checkmate::assertFlag(force_mirai)
  
  # Clear any existing clusters first
  unset_num_cores()
  # Set package/native thread budgets before workers are created so forked
  # children inherit the intended global thread settings.
  set_package_threads(num_cores)
  if (as.integer(num_cores) <= 1L) {
    return(invisible(NULL))
  }
  
  if (force_mirai || .Platform$OS.type != "unix") {
    if (!check_package_installed("mirai")) {
      stop("The 'mirai' package is required for parallelization on this system or when force_mirai = TRUE. Please install it.")
    }
    edi_env$mirai_has_been_used = TRUE
    # Initialize mirai daemons
    mirai::daemons(num_cores)
    edi_env$global_mirai_num_cores = num_cores
    # Each daemon inherited OMP_NUM_THREADS=N from the parent; cap them at 1 so
    # N daemons don't each spawn N OMP threads.  Main-process Rcpp OMP functions
    # use omp_set_num_threads() explicitly and are unaffected by this reset.
    tryCatch(
      mirai::everywhere({
        Sys.setenv(
          OMP_NUM_THREADS        = 1L,
          MKL_NUM_THREADS        = 1L,
          OPENBLAS_NUM_THREADS   = 1L,
          GOTO_NUM_THREADS       = 1L,
          VECLIB_MAXIMUM_THREADS = 1L,
          NUMEXPR_NUM_THREADS    = 1L
        )
        options(mc.cores = 1L)
        if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
        if (requireNamespace("fixest", quietly = TRUE)) suppressWarnings(try(fixest::setFixest_nthreads(1L), silent = TRUE))
      }),
      error = function(e) invisible(NULL)
    )
  } else {
    if (isTRUE(edi_env$mirai_has_been_used)) {
      stop(
        "Cannot switch from mirai-backed parallelism to fork-based parallelism in the same R session. ",
        "Restart R or keep using force_mirai = TRUE. This avoids the nng is not fork-reentrant safe panic."
      )
    }
    # Unix-like system, use forking
    cl = make_configured_fork_cluster(num_cores)
    edi_env$global_fork_cluster = cl
  }
  
  invisible(NULL)
}
#' Unset the number of cores and stop parallel clusters
#' 
#' This function stops any global fork or mirai clusters stored in the package
#' environment and resets the core count to serial execution.
#' @return Invisible \code{NULL}.
#' 
#' @examples
#' set_num_cores(2)
#' unset_num_cores()
#' @export
unset_num_cores = function() {
  # Handle fork cluster
  if (!is.null(edi_env$global_fork_cluster)) {
    cl = edi_env$global_fork_cluster
    edi_env$global_fork_cluster = NULL
    tryCatch(parallel::stopCluster(cl), error = function(e) invisible(NULL))
  }
  
  # Handle mirai cluster
  if (!is.null(edi_env$global_mirai_num_cores)) {
    if (check_package_installed("mirai")) {
      tryCatch(mirai::daemons(0), error = function(e) invisible(NULL)) # Stop daemons
    }
    edi_env$global_mirai_num_cores = NULL
    edi_env$mirai_has_been_used = FALSE
  }
  
  # Reset package threads to 1
  set_package_threads(1L)
  
  invisible(NULL)
}
#' Get the default parallel dispatch policy
#'
#' Returns EDI's built-in blocklist-first dispatch policy. This is the policy
#' @examples
#' get_bootstrap_dispatch_policy()
#' @return A named list describing the built-in dispatch policy.
#'
#' @export
get_parallel_dispatch_policy = function() {
  list(
    bootstrap = list(
      serial_inference_class_patterns = c(
        "^InferenceIncid",
        "^InferenceSurvival(?!.*KK)",
        "^InferenceAllKKWilcoxIVWC$"
      ),
      serial_response_types = c("incidence")
    ),
    rand_ci = list(
      serial_inference_class_patterns = c("^InferenceIncid"),
      serial_response_types = c("incidence")
    )
  )
}
edi_env$parallel_dispatch_policy_config = get_parallel_dispatch_policy()
#' Get the default bootstrap dispatch policy
#'
#' Returns EDI's built-in policy for choosing a default bootstrap type. Each
#' inference class can override the standard \code{"bca"} default via a regular
#' expression pattern.
#'
#' @details
#' The built-in overrides are empirical. Inference classes with repeated
#' @return A named list describing the default bootstrap type configuration.
#' @examples
#' get_bootstrap_dispatch_policy()
#' @export
get_bootstrap_dispatch_policy = function() {
  list(
    default_type = "bca",
    inference_class_overrides = c(
      "^InferenceContinLin$" = "percentile",
      "^InferenceIncidGCompRisk(Diff|Ratio)$" = "percentile",
      "^InferenceIncidKKGCompRisk(Diff|Ratio)$" = "percentile",
      "^InferenceIncidBinomialIdentityRiskDiff$" = "percentile",
      "^InferencePropGCompMeanDiff$" = "percentile",
      "^InferenceSurvivalDepCensTransformRegr$" = "percentile",
      "^InferenceSurvivalKKRankRegrIVWC$" = "percentile",
      "^InferenceOrdinalAdjCatLogitRegr$" = "percentile",
      "^InferenceSurvivalKKStratCoxPHOneLik$" = "percentile",
	      "^InferenceCountPoisson$" = "percentile",
	      "^InferenceCountRobustPoisson$" = "percentile",
	      "^InferenceCountQuasiPoisson$" = "percentile",
	      "^InferenceCountNegBin$" = "percentile",
	      "^InferenceCountZeroInflatedPoisson$" = "percentile",
	      "^InferenceCountZeroInflatedNegBin$" = "percentile",
	      "^InferenceCountHurdlePoisson$" = "percentile",
	      "^InferenceCountKKHurdlePoissonOneLik$" = "percentile",
	      "^InferenceCountKKCondPoissonOneLik$" = "percentile",
	      "^InferencePropZeroOneInflatedBetaRegr$" = "percentile",
	      "^InferencePropFractionalLogit$" = "percentile",
	      "^InferenceCountHurdleNegBin$" = "percentile",
      "^InferenceContinRobustRegr$" = "percentile"
    ),
    design_class_overrides = list(
      DesignFixedBlockedCluster = c(
        "^InferenceContinRobustRegr$" = "percentile",
        "^InferenceContinLin$" = "percentile",
        "^InferenceContinOLS$" = "percentile",
        "^InferenceContinKKOLSIVWC$" = "percentile",
        "^InferenceContinKKOLSOneLik$" = "percentile"
      )
    )
  )
}
edi_env$bootstrap_dispatch_policy_config = get_bootstrap_dispatch_policy()
#' Get the default optimization dispatch policy
#'
#' Returns EDI's built-in policy for choosing a default optimization algorithm.
#' @details
#' The built-in overrides are empirical. Regression models generally default to
#' @examples
#' get_optimization_dispatch_policy()
#' @return A named list describing the default optimization algorithm configuration.
#' @export
get_optimization_dispatch_policy = function() {
  list(
    default_alg = "newton_raphson",
    inference_class_overrides = c(
      "^InferenceCountKKGLMM$"  = "newton_raphson",
      "KKGLMM$"                 = "lbfgs",
      "KKWeibullFrailtyIVWC$"   = "lbfgs",
      "KKWeibullFrailtyOneLik$" = "lbfgs",
      "KKHurdlePoissonIVWC$"    = "lbfgs",
      "KKHurdlePoissonOneLik$"  = "lbfgs",
      "KKCondPoissonOneLik$"       = "lbfgs",
      "InferenceCountPoisson$"  = "irls",
      "InferenceCountQuasiPoisson$" = "irls",
      "InferenceCountRobustPoisson$" = "irls",
      "InferenceCountNegBin$" = "lbfgs",
      "InferenceIncidModifiedPoisson$" = "irls",
      "InferenceIncidLogRegr$" = "irls",
      "InferenceIncidProbitRegr$" = "irls",
      "InferencePropFractionalLogit$" = "irls",
      "InferencePropGCompAbstract$" = "irls",
      "InferencePropBetaRegr$" = "lbfgs",
      "InferenceOrdinalAdjCatLogitRegr$" = "lbfgs",
      "InferenceOrdinalContRatioRegr$" = "lbfgs",
      "InferenceOrdinalPropOddsRegr$" = "lbfgs",
      "InferenceOrdinalOrderedProbitRegr$" = "lbfgs",
      "InferenceOrdinalCauchitRegr$" = "lbfgs",
      "InferenceOrdinalCloglogRegr$" = "lbfgs",
      "InferenceSurvivalWeibullRegr$" = "lbfgs",
      "InferenceSurvivalStratCoxPHRegr$" = "newton_raphson",
      "InferenceSurvivalKKClaytonCopulaOneLik$" = "lbfgs",
      "InferenceIncidKKCondLogitPlusGLMMOneLik$"   = "lbfgs"
    )
  )
}
edi_env$optimization_dispatch_policy_config = get_optimization_dispatch_policy()
#' Get the default cold-start dispatch policy
#'
#' Returns EDI's built-in policy for the \code{smart_cold_start} default used
#' by each inference class.  A \code{TRUE} entry means the solver initialises
#' via an OLS warm-up before IRLS; \code{FALSE} means a plain zero-vector cold
#' start.  Benchmarks show the OLS warm-up is net-negative for logistic and
#' Poisson IRLS at typical trial sizes (the one extra OLS solve costs more than
#' the IRLS iterations it saves), so those families default to \code{FALSE}.
#' @return A named list with \code{default} (logical) and
#'   \code{inference_class_overrides} (named logical vector of regex ->
#'   logical).
#' @examples
#' get_cold_start_dispatch_policy()
#' @export
get_cold_start_dispatch_policy = function() {
  list(
    default = TRUE,
    inference_class_overrides = c(
      "^InferenceIncidLogRegr$"          = FALSE,
      "^InferencePropFractionalLogit$"   = FALSE,
      "^InferencePropGComp"              = FALSE,
      "^InferenceIncidGComp"             = FALSE,
      "^InferenceIncidKKGComp"           = FALSE,
      "^InferenceCountPoisson$"          = FALSE,
      "^InferenceCountQuasiPoisson$"     = FALSE,
      "^InferenceCountRobustPoisson$"    = FALSE,
      "^InferenceIncidModifiedPoisson$"  = FALSE
    )
  )
}
edi_env$cold_start_dispatch_policy_config = get_cold_start_dispatch_policy()
#' Update the cold-start dispatch policy
#'
#' @param policy Either \code{NULL} or a named list of policy overrides.
#' @param reset If \code{TRUE}, restore the built-in default policy.
#' @return Invisible \code{NULL} or the current policy configuration.
#' @examples
#' set_cold_start_dispatch_policy(reset = TRUE)
#' @export
set_cold_start_dispatch_policy = function(policy = NULL, reset = FALSE) {
  checkmate::assertFlag(reset)
  if (isTRUE(reset)) {
    edi_env$cold_start_dispatch_policy_config = get_cold_start_dispatch_policy()
    return(invisible(edi_env$cold_start_dispatch_policy_config))
  }
  if (is.null(policy)) {
    return(invisible(edi_env$cold_start_dispatch_policy_config))
  }
  checkmate::assertList(policy, names = "named")
  edi_env$cold_start_dispatch_policy_config = utils::modifyList(edi_env$cold_start_dispatch_policy_config, policy)
  invisible(NULL)
}
edi_cold_start_dispatch_policy = function(inference_class) {
  config = edi_env$cold_start_dispatch_policy_config
  inference_class = as.character(inference_class[[1]])
  overrides = config$inference_class_overrides
  default_val = if (isTRUE(config$default)) TRUE else FALSE
  if (!is.null(overrides) && length(overrides) > 0L) {
    for (pattern in names(overrides)) {
      if (is.na(pattern) || pattern == "") next
      if (grepl(pattern, inference_class, perl = TRUE)) {
        return(isTRUE(overrides[[pattern]]))
      }
    }
  }
  default_val
}
#' Get the default warm-start dispatch policy
#'
#' Returns EDI's built-in policy for choosing whether warm starts are enabled
#' during different resampling or simulation operations.
#'
#' @return A named list describing the warm start policy configuration.
#' @examples
#' get_warm_start_dispatch_policy()
#' @export
get_warm_start_dispatch_policy = function() {
  list(
    default = TRUE,
    jackknife = list(
      inference_class_overrides = c(
        "^InferenceSurvivalKKLWACoxPHIVWC$" = FALSE,
        "^InferenceSurvivalCoxPHRegr$" = FALSE,
        "^InferenceSurvivalStratCoxPHRegr$" = FALSE
      )
    ),
    non_param_boot = list(
      inference_class_overrides = c(
        "^InferenceCountNegBin$" = FALSE,
        "^InferenceSurvivalCoxPHRegr$" = FALSE,
        "^InferenceSurvivalStratCoxPHRegr$" = FALSE
      )
    ),
    bayesian_boot = list(
      inference_class_overrides = c(
        "^InferenceCountNegBin$" = FALSE,
        "^InferenceSurvivalCoxPHRegr$" = FALSE,
        "^InferenceSurvivalStratCoxPHRegr$" = FALSE
      )
    ),
    param_boot = list(
      inference_class_overrides = character(0)
    ),
    rand = list(
      inference_class_overrides = c(
        "^InferenceIncidKKCondLogitOneLik$" = FALSE,
        "^InferenceAllSimpleWilcox$" = FALSE
      )
    )
  )
}
edi_env$warm_start_dispatch_policy_config = get_warm_start_dispatch_policy()

#' Update the warm-start dispatch policy
#'
#' @param policy Either \code{NULL} or a named list of policy overrides.
#' @param reset If \code{TRUE}, restore the built-in default policy.
#' @return Invisible \code{NULL} or the current policy configuration.
#' @examples
#' set_warm_start_dispatch_policy(reset = TRUE)
#' @export
set_warm_start_dispatch_policy = function(policy = NULL, reset = FALSE) {
  checkmate::assertFlag(reset)
  if (isTRUE(reset)) {
    edi_env$warm_start_dispatch_policy_config = get_warm_start_dispatch_policy()
    return(invisible(edi_env$warm_start_dispatch_policy_config))
  }
  if (is.null(policy)) {
    return(invisible(edi_env$warm_start_dispatch_policy_config))
  }
  checkmate::assertList(policy, names = "named")
  edi_env$warm_start_dispatch_policy_config = utils::modifyList(edi_env$warm_start_dispatch_policy_config, policy)
  invisible(NULL)
}

edi_warm_start_dispatch_policy = function(inference_class, operation, n = NULL) {
  config = edi_env$warm_start_dispatch_policy_config
  inference_class = as.character(inference_class[[1]])
  operation = as.character(operation[[1]])
  n_val = suppressWarnings(as.integer(n))
  
  default_val = if (isTRUE(config$default)) TRUE else FALSE

  # Global disables: consistently negative across 3-4 n values or extreme single-n anomalies
  if (identical(operation, "bayesian_boot") &&
      (grepl("^InferenceContinKKGLMM$",                   inference_class, perl = TRUE) ||
       grepl("^InferenceSurvivalDepCensTransformRegr$",   inference_class, perl = TRUE) ||
       grepl("^InferenceOrdinalCloglogRegr$",             inference_class, perl = TRUE) ||
       grepl("^InferenceCountPoisson$",                   inference_class, perl = TRUE) ||
       grepl("^InferenceOrdinalPropOddsRegr$",            inference_class, perl = TRUE) ||
       grepl("^InferenceIncidKKNewcombeRiskDiff$",        inference_class, perl = TRUE) ||
       grepl("^InferenceOrdinalJonckheereTerpstraTest$",  inference_class, perl = TRUE) ||
       grepl("^InferenceOrdinalKKCLMMProbit$",            inference_class, perl = TRUE) ||
       grepl("^InferencePropFractionalLogit$",            inference_class, perl = TRUE) ||
       grepl("^InferencePropGCompMeanDiff$",              inference_class, perl = TRUE) ||
       grepl("^InferenceSurvivalWeibullRegr$",            inference_class, perl = TRUE))) {
    return(FALSE)
  }
  if (identical(operation, "bayesian_boot") &&
      grepl("^InferenceIncidBinomialIdentityRiskDiff$", inference_class, perl = TRUE)) {
    return(FALSE)
  }
  if (identical(operation, "jackknife") &&
      (grepl("^InferenceOrdinalContRatioRegr$",      inference_class, perl = TRUE) ||
       grepl("^InferenceOrdinalAdjCatLogitRegr$",    inference_class, perl = TRUE) ||
       grepl("^InferenceSurvivalRestrictedMeanDiff$", inference_class, perl = TRUE))) {
    return(FALSE)
  }
  if (identical(operation, "non_param_boot") &&
      grepl("^InferencePropZeroOneInflatedBetaRegr$", inference_class, perl = TRUE)) {
    return(FALSE)
  }
  if (identical(operation, "rand") &&
      grepl("^InferenceSurvivalKKLWACoxPHIVWC$", inference_class, perl = TRUE)) {
    return(FALSE)
  }

  if (!is.na(n_val) && n_val < 200L) {
    if (identical(operation, "non_param_boot") &&
        (grepl("^InferenceOrdinalContRatioRegr$",         inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKCLMMCauchit$",         inference_class, perl = TRUE) ||
         grepl("^InferencePropKKGEE$",                    inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKCondAdjCatLogitRegr$", inference_class, perl = TRUE) ||
         grepl("^InferenceCountQuasiPoisson$",            inference_class, perl = TRUE) ||
         grepl("^InferencePropKKQuantileRegrOneLik$",     inference_class, perl = TRUE) ||
         grepl("^InferenceCountHurdlePoisson$",           inference_class, perl = TRUE) ||
         grepl("^InferenceCountZeroInflatedPoisson$",     inference_class, perl = TRUE) ||
         grepl("^InferenceIncidBinomialIdentityRiskDiff$",inference_class, perl = TRUE) ||
         grepl("^InferenceIncidGCompRiskRatio$",          inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKGEE$",                   inference_class, perl = TRUE) ||
         grepl("^InferencePropBetaRegr$",                 inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "rand") &&
        (grepl("^InferencePropBetaRegr$",                 inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKCondAdjCatLogitRegr$", inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalKKClaytonCopulaOneLik$",inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalKKStratCoxPHOneLik$",   inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalKKWeibullFrailtyIVWC$", inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "jackknife") &&
        grepl("^InferenceSurvivalGehanWilcox$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
    if (identical(operation, "bayesian_boot") &&
        grepl("^InferenceContinQuantileRegr$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
  }

  if (!is.na(n_val) && n_val < 500L) {
    if (identical(operation, "rand") &&
        (grepl("^InferenceContinKKOLSIVWC$",       inference_class, perl = TRUE) ||
         grepl("^InferenceContinKKRobustRegrIVWC$", inference_class, perl = TRUE) ||
         grepl("^InferencePropKKQuantileRegrIVWC$", inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalContRatioRegr$",   inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalCoxPHRegr$",      inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKGEE$",             inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "non_param_boot") &&
        (grepl("^InferenceCountZeroInflatedNegBin$",     inference_class, perl = TRUE) ||
         grepl("^InferenceIncidLogBinomial$",            inference_class, perl = TRUE) ||
         grepl("^InferenceAllSimpleWilcox$",             inference_class, perl = TRUE) ||
         grepl("^InferenceContinKKGLMM$",               inference_class, perl = TRUE) ||
         grepl("^InferenceCountPoisson$",               inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalContRatioRegr$",       inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalDepCensTransformRegr$",inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "rand") &&
        grepl("^InferenceIncidLogBinomial$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
    if (identical(operation, "bayesian_boot") &&
        (grepl("^InferenceOrdinalKKCondAdjCatLogitRegr$", inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalGehanWilcox$",           inference_class, perl = TRUE) ||
         grepl("^InferenceCountKKGLMM$",                   inference_class, perl = TRUE) ||
         grepl("^InferenceCountHurdleNegBin$",             inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalContRatioRegr$",          inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalKKStratCoxPHOneLik$",    inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "jackknife") &&
        (grepl("^InferenceSurvivalLogRank$",              inference_class, perl = TRUE) ||
         grepl("^InferenceContinKKGLMM$",                 inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalCauchitRegr$",           inference_class, perl = TRUE) ||
         grepl("^InferencePropBetaRegr$",                 inference_class, perl = TRUE))) {
      return(FALSE)
    }
  }

  if (!is.na(n_val) && n_val >= 200L && n_val < 500L) {
    if (identical(operation, "rand") &&
        grepl("^InferenceContinKKQuantileRegrOneLik$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
    if (identical(operation, "jackknife") &&
        grepl("^InferenceIncidKKCondLogitPlusGLMMIVWC$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
  }

  if (!is.na(n_val) && n_val < 1000L) {
    if (identical(operation, "rand") &&
        (grepl("^InferenceContinKKGLMM$",         inference_class, perl = TRUE) ||
         grepl("^InferenceContinQuantileRegr$",    inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKModifiedPoisson$",inference_class, perl = TRUE) ||
         grepl("^InferencePropGCompMeanDiff$",     inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "bayesian_boot") &&
        grepl("^InferenceOrdinalKKCLMMCauchit$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
    if (identical(operation, "jackknife") &&
        grepl("^InferenceCountNegBin$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
    if (identical(operation, "non_param_boot") &&
        grepl("^InferenceCountHurdleNegBin$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
  }

  if (!is.na(n_val) && n_val >= 200L) {
    # At n>=200, KKHurdlePoissonOneLik rand warm starts cause the C++ GLMM to
    # fail convergence and fall back to slow glmmTMB (benchmark shows -75% at n=200)
    if (identical(operation, "rand") &&
        grepl("^InferenceCountKKHurdlePoissonOneLik$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
  }

  if (!is.na(n_val) && n_val >= 500L) {
    if (identical(operation, "bayesian_boot") &&
        (grepl("^InferenceSurvivalKKClaytonCopulaOneLik$",  inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKCondLogitPlusGLMMOneLik$", inference_class, perl = TRUE) ||
         grepl("^InferenceIncidModifiedPoisson$",           inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKCLMM$",                  inference_class, perl = TRUE) ||
         grepl("^InferencePropZeroOneInflatedBetaRegr$",    inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "rand") &&
        (grepl("^InferenceContinRobustRegr$",               inference_class, perl = TRUE) ||
         grepl("^InferenceIncidBinomialIdentityRiskDiff$",  inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKCondLogitPlusGLMMOneLik$", inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKCondLogitPlusGLMMIVWC$",   inference_class, perl = TRUE) ||
         grepl("^InferenceIncidGCompRiskRatio$",            inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKCLMM$",                  inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "non_param_boot") &&
        (grepl("^InferenceAllKKWilcoxIVWC$",       inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalCloglogRegr$",     inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKGLMM$",          inference_class, perl = TRUE) ||
         grepl("^InferencePropFractionalLogit$",    inference_class, perl = TRUE) ||
         grepl("^InferencePropKKQuantileRegrIVWC$", inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalKMDiff$",         inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "jackknife") &&
        (grepl("^InferenceContinQuantileRegr$",  inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalGCompMeanDiff$", inference_class, perl = TRUE))) {
      return(FALSE)
    }
  }

  if (!is.na(n_val) && n_val >= 1000L) {
    if (identical(operation, "rand") &&
        grepl("^InferenceCountPoisson$", inference_class, perl = TRUE)) {
      return(FALSE)
    }
    if (identical(operation, "jackknife") &&
        (grepl("^InferenceSurvivalKKClaytonCopulaIVWC$",  inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalOrderedProbitRegr$",      inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalDepCensTransformRegr$",  inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalGehanWilcox$",           inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalWeibullRegr$",           inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "bayesian_boot") &&
        (grepl("^InferenceIncidGCompRiskRatio$",          inference_class, perl = TRUE) ||
         grepl("^InferenceContinQuantileRegr$",           inference_class, perl = TRUE) ||
         grepl("^InferenceIncidRiskDiff$",                inference_class, perl = TRUE) ||
         grepl("^InferenceCountQuasiPoisson$",            inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKGCompRiskRatio$",        inference_class, perl = TRUE) ||
         grepl("^InferenceCountKKCondPoissonOneLik$",     inference_class, perl = TRUE) ||
         grepl("^InferenceCountKKGLMM$",                  inference_class, perl = TRUE) ||
         grepl("^InferenceCountKKHurdlePoissonOneLik$",   inference_class, perl = TRUE) ||
         grepl("^InferenceIncidExactBinomial$",           inference_class, perl = TRUE) ||
         grepl("^InferenceIncidProbitRegr$",              inference_class, perl = TRUE) ||
         grepl("^InferenceIncidenceExactZhang$",          inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKGLMM$",                inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalOrderedProbitRegr$",     inference_class, perl = TRUE))) {
      return(FALSE)
    }
    if (identical(operation, "non_param_boot") &&
        (grepl("^InferenceIncidKKCondLogitPlusGLMMOneLik$",inference_class, perl = TRUE) ||
         grepl("^InferenceIncidKKGCompRiskDiff$",          inference_class, perl = TRUE) ||
         grepl("^InferenceIncidRiskDiff$",                 inference_class, perl = TRUE) ||
         grepl("^InferenceOrdinalKKCLMM$",                 inference_class, perl = TRUE) ||
         grepl("^InferenceSurvivalRestrictedMeanDiff$",    inference_class, perl = TRUE))) {
      return(FALSE)
    }
  }

  if (operation %in% names(config)) {
    op_cfg = config[[operation]]
    if (is.list(op_cfg)) {
      overrides = op_cfg$inference_class_overrides
      if (!is.null(overrides) && length(overrides) > 0L) {
        for (pattern in names(overrides)) {
          if (is.na(pattern) || pattern == "") next
          if (grepl(pattern, inference_class, perl = TRUE)) {
            return(isTRUE(overrides[[pattern]]))
          }
        }
      }
    }
  }
  default_val
}
#' Update the parallel dispatch policy
#'
#' EDI uses an empirical, blocklist-first dispatch policy to decide when an
#' inference routine should be forced serial even if multiple cores are
#' available. This function lets the user update that default policy without
#' editing package internals.
#'
#' @details
#' The policy can be updated in two ways:
#' \itemize{
#'   \item Pass a named list to merge with the built-in default policy
#'   configuration. Supported top-level keys are \code{bootstrap} and
#'   \code{rand_ci}, and each key may contain
#'   \code{serial_inference_class_patterns} and
#'   \code{serial_response_types}.
#'   \item Pass a custom function with signature
#'   \code{function(inference_class, response_type, operation)} that returns a
#'   list with at least \code{force_serial} and \code{reason}.
#' }
#'
#' Use \code{reset = TRUE} to restore the built-in default policy. Do not expect
#' a universal "more cores is faster" rule.
#'
#' @param policy Either \code{NULL}, a named list of policy overrides, or a
#'   custom function. If \code{NULL} and \code{reset = FALSE}, the current
#' @examples
#' set_parallel_dispatch_policy(reset = TRUE)
#' @param reset If \code{TRUE}, restore the built-in default policy and remove
#'   any custom function override.
#'
#' @return Invisible \code{NULL} or the current policy configuration.
#'
#' @export
set_parallel_dispatch_policy = function(policy = NULL, reset = FALSE) {
  checkmate::assertFlag(reset)
  if (isTRUE(reset)) {
    edi_env$parallel_dispatch_policy_config = get_parallel_dispatch_policy()
    edi_env$parallel_dispatch_policy_override = NULL
    return(invisible(edi_env$parallel_dispatch_policy_config))
  }
  if (is.null(policy)) {
    return(invisible(edi_env$parallel_dispatch_policy_config))
  }
  if (is.function(policy)) {
    edi_env$parallel_dispatch_policy_override = policy
    return(invisible(NULL))
  }
  checkmate::assertList(policy, names = "named")
  current_config = edi_env$parallel_dispatch_policy_config
  for (nm in names(policy)) {
    if (!nm %in% names(current_config)) {
      stop("Unknown policy section: ", nm, call. = FALSE)
    }
    if (!is.list(policy[[nm]])) {
      stop("Policy section '", nm, "' must be a list.", call. = FALSE)
    }
    current_config[[nm]] = utils::modifyList(current_config[[nm]], policy[[nm]])
  }
  edi_env$parallel_dispatch_policy_config = current_config
  edi_env$parallel_dispatch_policy_override = NULL
  invisible(NULL)
}
#' Update the optimization dispatch policy
#'
#' @examples
#' set_optimization_dispatch_policy(reset = TRUE)
#' @param policy Either \code{NULL} or a named list of policy overrides.
#' @param reset If \code{TRUE}, restore the built-in default policy.
#'
#' @return Invisible \code{NULL} or the current policy configuration.
#' @export
set_optimization_dispatch_policy = function(policy = NULL, reset = FALSE) {
  checkmate::assertFlag(reset)
  if (isTRUE(reset)) {
    edi_env$optimization_dispatch_policy_config = get_optimization_dispatch_policy()
    return(invisible(edi_env$optimization_dispatch_policy_config))
  }
  if (is.null(policy)) {
    return(invisible(edi_env$optimization_dispatch_policy_config))
  }
  checkmate::assertList(policy, names = "named")
  edi_env$optimization_dispatch_policy_config = utils::modifyList(edi_env$optimization_dispatch_policy_config, policy)
  invisible(NULL)
}
# Internal helper to get the global fork cluster
get_global_fork_cluster = function() {
  edi_env$global_fork_cluster
}
# Internal helper to get the global mirai core count
get_global_mirai_cores = function() {
  edi_env$global_mirai_num_cores
}
get_num_cores = function() {
  if (!is.null(edi_env$num_cores_override)) return(edi_env$num_cores_override)
  cl = get_global_fork_cluster()
  if (!is.null(cl)) return(length(cl))
  mirai_cores = get_global_mirai_cores()
  if (!is.null(mirai_cores)) return(mirai_cores)
  1L
}
# Internal helper for empirical parallel dispatch policy.
# This is intentionally conservative and only forces serial execution for
# inference families that have shown repeat multicore regressions in the package
# benchmark suite. Everything else remains eligible for the usual warmup-based
# parallel heuristics. Do not expect a universal "more cores is faster" rule.
edi_parallel_dispatch_policy = function(inference_class, response_type, operation) {
  if (is.function(edi_env$parallel_dispatch_policy_override)) {
    return(edi_env$parallel_dispatch_policy_override(inference_class, response_type, operation))
  }
  inference_class = as.character(inference_class[[1]])
  response_type = as.character(response_type[[1]])
  operation = as.character(operation[[1]])
  config = edi_env$parallel_dispatch_policy_config
  matches_any = function(value, patterns) {
    if (is.null(patterns) || length(patterns) == 0L) return(FALSE)
    any(vapply(patterns, function(pattern) grepl(pattern, value, perl = TRUE), logical(1)))
  }
  reason = NULL
  if (identical(operation, "bootstrap")) {
    op_cfg = config$bootstrap
    if (matches_any(inference_class, op_cfg$serial_inference_class_patterns) ||
        matches_any(response_type, op_cfg$serial_response_types)) {
      reason = "bootstrap is forced serial by benchmark policy"
    }
  } else if (identical(operation, "rand_ci")) {
    op_cfg = config$rand_ci
    if (matches_any(inference_class, op_cfg$serial_inference_class_patterns) ||
        matches_any(response_type, op_cfg$serial_response_types)) {
      reason = "randomization confidence intervals are forced serial by benchmark policy"
    }
  }
  list(
    force_serial = !is.null(reason),
    reason = reason,
    inference_class = inference_class,
    response_type = response_type,
    operation = operation
  )
}
edi_bootstrap_dispatch_policy = function(inference_class, object = NULL) {
  config = edi_env$bootstrap_dispatch_policy_config
  inference_class = as.character(inference_class[[1]])
  overrides = config$inference_class_overrides
  design_overrides = config$design_class_overrides
  default_type = config$default_type
  if (is.null(default_type)) default_type = "bca"
  finalize_type = function(type) {
    type = tolower(type)
    if (identical(type, "bca") && !is.null(object)) {
      bca_unavailable = isTRUE(tryCatch(
        object$.__enclos_env__$private$jackknife_block_size_gt_one_unsupported(unit = "auto"),
        error = function(e) FALSE
      ))
      if (bca_unavailable) return("percentile")
    }
    type
  }
  if (!is.null(object) && !is.null(design_overrides) && length(design_overrides) > 0L) {
    des_obj = tryCatch(object$.__enclos_env__$private$des_obj, error = function(e) NULL)
    if (!is.null(des_obj)) {
      for (design_class in names(design_overrides)) {
        if (is(des_obj, design_class)) {
          design_map = design_overrides[[design_class]]
          for (pattern in names(design_map)) {
            if (is.na(pattern) || pattern == "") next
            if (grepl(pattern, inference_class, perl = TRUE)) {
              return(finalize_type(design_map[[pattern]]))
            }
          }
        }
      }
    }
  }
  if (!is.null(overrides) && length(overrides) > 0L) {
    for (pattern in names(overrides)) {
      if (is.na(pattern) || pattern == "") next
      if (grepl(pattern, inference_class, perl = TRUE)) {
        return(finalize_type(overrides[[pattern]]))
      }
    }
  }
  finalize_type(default_type)
}
edi_optimization_dispatch_policy = function(inference_class) {
  config = edi_env$optimization_dispatch_policy_config
  inference_class = as.character(inference_class[[1]])
  overrides = config$inference_class_overrides
  default_alg = config$default_alg
  if (is.null(default_alg)) default_alg = "newton_raphson"
  
  if (!is.null(overrides) && length(overrides) > 0L) {
    for (pattern in names(overrides)) {
      if (is.na(pattern) || pattern == "") next
      if (grepl(pattern, inference_class, perl = TRUE)) {
        return(tolower(overrides[[pattern]]))
      }
    }
  }
  tolower(default_alg)
}
# Internal helper
set_package_threads = function(num_cores) {
  # Ensure it's an integer
  num_cores = as.integer(num_cores)
  # Skip all expensive syscalls if threads are already set to this value.
  # Sys.setenv and the BLAS/OMP setters are relatively slow; calling them on
  # every serial replication causes a visible pause between reps.
  last = getOption(".edi_last_set_threads")
  if (identical(last, num_cores)) return(invisible(NULL))
  # R packages with global thread setters
	  if (check_package_installed("data.table")) {
	    data.table::setDTthreads(num_cores)
	  }
	  if (check_package_installed("fixest")) {
	    suppressWarnings(try(fixest::setFixest_nthreads(num_cores), silent = TRUE))
	  }
  # Environment variables for OpenMP and BLAS/LAPACK
  # This helps prevent thread explosion in child processes
  # that call multi-threaded native libraries.
  # Sys.setenv is relatively slow, so we only do it if needed.
  # We include as many as possible to catch various libraries.
  Sys.setenv(
    OMP_NUM_THREADS        = num_cores,
    MKL_NUM_THREADS        = num_cores,
    OPENBLAS_NUM_THREADS   = num_cores,
    GOTO_NUM_THREADS       = num_cores,
    BLAS_NUM_THREADS       = num_cores,
    LAPACK_NUM_THREADS     = num_cores,
    TAO_NUM_THREADS        = num_cores,
    VECLIB_MAXIMUM_THREADS = num_cores,
    NUMEXPR_NUM_THREADS    = num_cores,
    OMP_THREAD_LIMIT       = num_cores,
    OMP_DYNAMIC            = "FALSE",
    OMP_NESTED             = "FALSE"
  )
  
  # Direct C++ control for OpenMP and Eigen (more robust than env vars after startup)
  try(set_omp_num_threads_cpp(num_cores), silent = TRUE)
  # BLAS and OpenMP control via RhpcBLASctl (now a required dependency)
  try({
    RhpcBLASctl::blas_set_num_threads(num_cores)
    RhpcBLASctl::omp_set_num_threads(num_cores)
  }, silent = TRUE)
  
  # Also set R options for parallel/pbmcapply
  options(mc.cores = num_cores)
  
  options(".edi_last_set_threads" =   num_cores)
  invisible(NULL)
}

.create_match_dummies = function(m_vec) {
  m_vec = as.integer(m_vec)
  m_vec[is.na(m_vec)] = 0L
  max_m = max(m_vec, 0L)
  if (max_m == 0L) return(NULL)
  mm = stats::model.matrix(~ factor(m_vec) + 0)
  colnames(mm) = paste0("match_", 0:max_m)
  mm
}
