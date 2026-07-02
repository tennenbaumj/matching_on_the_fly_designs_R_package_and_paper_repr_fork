# simulations.R
#
# Simulation framework for comparing experimental designs and inference methods
# in the EDI (Experimental Design with Inference) package.
#
# Usage:
#   devtools::load_all("EDI")
#   source("simulations.R")
#
#   designs = list(
#     DesignSeqOneByOneKK21 = list(lambda = 0.5),
#     DesignSeqOneByOneUrn  = list(alpha = 2, beta = 2)
#   )
#   inf_cls = list(
#     InferenceContinOLS = list(max_resample_attempts = 25L),
#     InferenceContinKKOLSIVWC
#   )
#   sim = SimulationFramework$new(
#     response_type    = "continuous",
#     design_classes_and_params = designs,
#     inference_classes_and_params = inf_cls,
#     n = 100, p = 5, cond_exp_func_model = "linear", Nrep = 50, betaT = 1,
#     inference_types_and_params = list(asymp_pval = list(delta = 0))
#   )
#   sim$run()
#   SimulationFrameworkReport$new(sim)$summarize()
#' Generate Synthetic Simulation Covariates and Continuous Response
#'
#' @description A helper function to generate synthetic covariates and a latent continuous response
#' identical to the logic used within \code{SimulationFramework}.  Covariates may be
#' supplied directly via \code{X_mat} or drawn randomly via \code{cov_draw_method};
#' exactly one of the two must be non-\code{NULL}.
#'
#' @param n Integer. Sample size (number of rows).
#' @param p Integer. Number of covariates (number of columns).
#' @param cond_exp_func_model Character scalar. Either \code{"linear"} (latent response is a
#'   weighted linear combination of covariates) or \code{"nonlinear"} (Friedman 1991
#'   function applied to the first five covariates; requires \code{p >= 5}).
#' @param norm_sq_beta_vec Positive numeric scalar. The desired squared Euclidean norm
#'   of the coefficient vector, i.e. \code{sum(beta^2)}.  The coefficient vector (or
#'   the overall Friedman scale) is rescaled so that this quantity equals
#'   \code{norm_sq_beta_vec}.  Default \code{1}.
#' @param X_mat Numeric matrix of dimensions \code{n x p}, or \code{NULL} (default).
#'   When supplied, this matrix is used directly as the covariate matrix and
#'   \code{cov_draw_method} must be \code{NULL}.
#' @param cov_draw_method A function used to draw \code{n * p} i.i.d. covariate
#'   values, or \code{NULL}.  The function must accept the number of draws as its
#'   first positional argument followed by any named arguments in
#'   \code{cov_draw_method_args}.  Default \code{stats::rnorm}.
#'   Must be \code{NULL} when \code{X_mat} is supplied.
#' @param cov_draw_method_args Named list of additional arguments forwarded to
#'   \code{cov_draw_method} beyond the sample-size first argument.
#'   Default is \code{list(mean = 0, sd = 1)}.
#'
#' @return A list with two elements: \code{X} (a data frame of covariates) and
#'   \code{y_cont} (a numeric vector of the latent continuous response).
#'
#' @examples
#' generate_covariate_dataset(n = 10, p = 5)
#' @export
generate_covariate_dataset = function(n, p,
                                      cond_exp_func_model  = c("linear", "nonlinear"),
                                      norm_sq_beta_vec     = 1,
                                      X_mat                = NULL,
                                      cov_draw_method      = stats::rnorm,
                                      cov_draw_method_args = list(mean = 0, sd = 1)) {
  cond_exp_func_model = match.arg(cond_exp_func_model)
  user_supplied_X   = !is.null(X_mat)
  user_supplied_cov = !is.null(cov_draw_method)
  if (user_supplied_X && user_supplied_cov) {
    stop("generate_covariate_dataset: supply exactly one of 'X_mat' or 'cov_draw_method', not both.")
  }
  if (!user_supplied_X && !user_supplied_cov) {
    stop("generate_covariate_dataset: one of 'X_mat' or 'cov_draw_method' must be non-NULL.")
  }
  if (user_supplied_X) {
    if (!is.matrix(X_mat)) X_mat = as.matrix(X_mat)
  } else {
    X_mat = matrix(do.call(cov_draw_method, c(list(n * p), cov_draw_method_args)), nrow = n, ncol = p)
  }
  
  if (is.null(colnames(X_mat))) {
    colnames(X_mat) = paste0("x", seq_len(p))
  }
  X = data.table::as.data.table(X_mat)
  
  if (cond_exp_func_model == "linear") {
    beta_x = seq(1, -1, length.out = p)
    beta_x = beta_x * sqrt(norm_sq_beta_vec / sum(beta_x^2))
    y_cont = as.numeric(X_mat %*% beta_x)
  } else {
    if (p < 5L) stop("Friedman nonlinear cond_exp_func_model requires p >= 5.")
    # Friedman (1991) function — first five covariates, x in [-1,1]
    beta_friedman = c(10, 20, 10, 5)
    scale_factor  = sqrt(norm_sq_beta_vec / sum(beta_friedman^2))
    y_cont = scale_factor * (
      10 * sin(pi * X_mat[, 1L] * X_mat[, 2L]) +
      20 * (X_mat[, 3L] - 0.5)^2 +
      10 *  X_mat[, 4L] +
       5 *  X_mat[, 5L]
    )
  }
  list(X = X, y_cont = y_cont)
}
#' Transform continuous latent signal to the response type scale
#'
#' @description A helper function to transform a latent continuous signal to the scale
#' appropriate for a given \code{response_type}, identical to the logic used
#' within \code{SimulationFramework} but not used within \code{SimulationFramework}. 
#'
#' @param y_cont Numeric vector. The latent continuous response signal.
#' @param response_type Character scalar. One of \code{"continuous"}, \code{"incidence"},
#'   \code{"proportion"}, \code{"count"}, \code{"survival"}, \code{"ordinal"}.
#' @param n_ordinal_levels Positive integer. Number of ordinal categories when
#'   \code{response_type = "ordinal"}. Default \code{4L}.
#' @param proportion_epsilon Numeric scalar. Small value added to proportion to avoid 0 and 1. Default \code{1e-6}.
#' @param survival_min_time Numeric scalar. Minimum survival time and shift. Default \code{0.1}.
#' @param count_min_rate Integer scalar. Minimum baseline rate for count response. Default \code{0L}.
#' @param count_shift Numeric scalar. Constant added to counts after zero-centering. Default \code{0}.
#'
#' @examples
#' transform_cont_y_based_on_response_type(rnorm(10), 'incidence')
#' @return A numeric vector of transformed responses on the appropriate scale.
#'
#' @export
transform_cont_y_based_on_response_type = function(
  y_cont, 
  response_type, 
  n_ordinal_levels   = 4L,
  proportion_epsilon = 1e-6,
  survival_min_time  = 0.1,
  count_min_rate     = 0L,
  count_shift        = 0
) {
  y_sh = as.numeric(y_cont - mean(y_cont))
  switch(response_type,
    continuous = y_sh,
    incidence  = stats::plogis(y_sh),
    proportion = {
      y_tmp = y_cont - min(y_cont) + proportion_epsilon
      y_tmp / (max(y_tmp) + proportion_epsilon)
    },
    count      = pmax(count_min_rate, round(y_cont - min(y_cont) + count_shift)),
    survival   = pmax(survival_min_time, y_sh - min(y_sh) + survival_min_time),
    ordinal    = as.integer(cut(
      y_cont,
      breaks = unique(stats::quantile(y_cont, probs = seq(0, 1, length.out = n_ordinal_levels + 1L))),
      include.lowest = TRUE, labels = FALSE
    )),
    stop("Unknown response_type: ", response_type)
  )
}
# Internal helper for null-coalescing
`%||%` = function(a, b) if (!is.null(a)) a else b
# Walk an R6 class generator's inheritance chain to find the first class
# that defines initialize(), and return that function.
get_r6_init_fn = function(r6gen) {
  gen = r6gen
  while (!is.null(gen)) {
    fn = tryCatch(gen$public_methods$initialize, error = function(e) NULL)
    if (!is.null(fn)) return(fn)
    gen = tryCatch(gen$get_inherit(), error = function(e) NULL)
  }
  NULL
}
#' Simulation Framework for Experimental Designs and Inference Methods
#'
#' @description An R6 class for benchmarking experimental designs and inference methods by
#' Monte Carlo simulation. Each replication generates synthetic covariates and
#' responses, runs every requested \code{(design, inference)} pair, and records
#' point estimates, confidence intervals, and p-values. Raw and aggregated
#' results are available through \code{SimulationFrameworkReport}.
#'
#' @details
#' Covariates are drawn independently from \eqn{\mathrm{Uniform}(0, 1)}.
#' \itemize{
#'   \item \code{cond_exp_func_model = "linear"}: the base continuous signal is
#'     \eqn{y = X\beta} where \eqn{\beta} is evenly spaced from 1 to \eqn{-1}.
#'   \item \code{cond_exp_func_model = "nonlinear"}: the Friedman (1991) function
#'     \eqn{10\sin(\pi x_1 x_2) + 20(x_3-0.5)^2 + 10x_4 + 5x_5};
#'     requires \eqn{p \ge 5}.
#' }
#' The continuous base signal is transformed to the scale appropriate for
#' \code{response_type}.  Treatment effects are applied per-subject: additive
#' on the linear/logit/ordinal scale, log-multiplicative for count and survival.
#'
#' For each \code{(design, inference)} pair the framework runs whichever of the
#' following are supported by the inference class:
#' \itemize{
#'   \item \strong{asymptotic} (\code{InferenceAsymp} subclasses): Wald CI and
#'     p-value.
#'   \item \strong{bootstrap} (\code{InferenceNonParamBootstrap} subclasses): percentile CI
#'     and p-value.
#'   \item \strong{randomisation} (\code{InferenceRand} subclasses): p-value;
#'     additionally a test-inversion CI for \code{continuous},
#'     \code{proportion}, and \code{count} response types
#'     (\code{InferenceRandCI} subclasses).
#' }
#' Incompatible \code{(design, inference)} pairs (e.g.\ a KK-specific inference
#' class with a non-KK design) are silently skipped via \code{tryCatch}.
#'
#' Reported summary metrics include:
#' \itemize{
#'   \item \strong{MSE}: \eqn{\overline{(\hat\beta_T - \beta_T)^2}} over reps
#'     with a finite point estimate.
#'   \item \strong{coverage}: proportion of reps where \eqn{\beta_T} lies
#'     inside the CI (\code{NA} when no CI is available for that inference type).
#'   \item \strong{power}: proportion of p-values \eqn{< \alpha}; equals the
#'     empirical type-I error rate when \code{betaT = 0}.
#' }
#'
#' @examples
#' \donttest{
#' # Simple simulation with two designs and two inference methods
#' sim = SimulationFramework$new(
#'   response_type = "continuous",
#'   design_classes_and_params = list(
#'     DesignSeqOneByOneKK21 = list(lambda = 0.5),
#'     DesignSeqOneByOneBernoulli = list()
#'   ),
#'   inference_classes_and_params = list(
#'     InferenceContinOLS = list(),
#'     InferenceContinKKOLSIVWC = list()
#'   ),
#'   n = 100, p = 5, Nrep = 10, betaT = 1
#' )
#' sim$run()
#' report = SimulationFrameworkReport$new(sim)
#' report$summarize()
#' }
#' @export
SimulationFramework = R6::R6Class("SimulationFramework",
  lock_objects = FALSE,
  # ── public ─────────────────────────────────────────────────────────────────
  public = list(
    #' @description Create a new \code{SimulationFramework}.
    #'
    #' @param response_type \strong{(required)} Character scalar or vector.  The type of
    #'   outcome variable.  One of \code{"continuous"}, \code{"incidence"},
    #'   \code{"proportion"}, \code{"count"}, \code{"survival"},
    #'   \code{"ordinal"}.
    #'
    #' @param design_classes_and_params \code{NULL} (default) or a list
    #'   describing design classes and optional constructor parameters.
    #'   Unnamed R6 class generators use default parameters, for example
    #'   \code{list(DesignSeqOneByOneKK21, DesignFixedBernoulli)}.  Named
    #'   entries use the entry name as the design class and the value as the
    #'   parameter list, for example
    #'   \code{list(DesignSeqOneByOneUrn = list(alpha = 2, beta = 2))}.
    #'   Duplicate named entries are allowed for repeated designs with
    #'   different parameters.
    #'   Each generator must be constructable with only \code{response_type} and
    #'   \code{n} plus any extra params supplied in this list.
    #'   \code{NULL} uses the package's standard design set.
    #'   Designs requiring \code{strata_cols}, \code{cluster_col}, or
    #'   \code{factors} have sensible defaults auto-injected (first covariate
    #'   column; second for \code{cluster_col}; \code{list(treatment=2)} for
    #'   \code{factors}) when not supplied in the parameter list.
    #'   Example:
    #'   \preformatted{design_classes_and_params = list(
    #'   DesignSeqOneByOneKK21 = list(lambda = 0.5, t_0_pct = 0.1),
    #'   DesignSeqOneByOneUrn  = list(alpha  = 2,   beta    = 2),
    #'   DesignFixedBernoulli  # default params
    #' )}
    #'   Commonly useful design constructor parameters:
    #'   \describe{
    #'     \item{\code{lambda}}{Matching-weight decay for
    #'       \code{KK14} / \code{KK21} / \code{KK21stepwise}.}
    #'     \item{\code{t_0_pct}}{Burn-in fraction for
    #'       \code{KK14} / \code{KK21} / \code{KK21stepwise}.}
    #'     \item{\code{morrison}}{Logical; Morrison correction for \code{KK14}.}
    #'     \item{\code{alpha}, \code{beta}}{Shape parameters for
    #'       \code{DesignSeqOneByOneUrn}.}
    #'     \item{\code{preferred_num_bins_for_continuous_covariate}}{Bin count for
    #'       \code{DesignFixedBlocking} and \code{DesignFixedBlockedCluster}.}
    #'     \item{\code{B_target}}{Target number of blocks for \code{DesignFixedBlocking}.}
    #'   }
    #'
    #' @param inference_classes_and_params \code{NULL} (default) or a list
    #'   describing inference classes and optional constructor parameters.
    #'   Unnamed R6 class generators use default parameters, for example
    #'   \code{list(InferenceContinOLS, InferenceContinKKOLSIVWC)}.
    #'   Named entries use the entry name as the inference class and the value
    #'   as the constructor parameter list, for example
    #'   \code{list(InferenceContinOLS = list(max_resample_attempts = 25L))}.
    #'   Duplicate named entries are allowed for repeated inference classes with
    #'   different parameters. Supplied parameters must be accepted by the
    #'   inference class constructor.
    #'   \code{NULL} selects a curated set for the given
    #'   \code{response_type}: several universal classes that work with any
    #'   design, plus representative KK-specific classes (silently skipped for
    #'   non-KK designs at runtime).
    #'
    #' @param n Integer scalar or vector.  Sample size per simulation replication.  Default
    #'   \code{100}.
    #'
    #' @param p Integer scalar or vector.  Number of covariates.  Must be \eqn{\ge 5} when
    #'   \code{cond_exp_func_model = "nonlinear"}.  Default \code{5}.
    #'
    #' @param cond_exp_func_model Character scalar or vector.  How the latent continuous signal is
    #'   constructed before transformation to the \code{response_type} scale.
    #'   \describe{
    #'     \item{\code{"linear"}}{Linear combination \eqn{X\beta} with
    #'       coefficients evenly spaced from 1 to \eqn{-1}.}
    #'     \item{\code{"nonlinear"}}{Friedman (1991) function
    #'       \eqn{10\sin(\pi x_1 x_2)+20(x_3-0.5)^2+10x_4+5x_5};
    #'       requires \eqn{p \ge 5}.}
    #'   }
    #'   Default \code{"linear"}.
    #'
    #' @param Nrep Positive integer.  Number of Monte Carlo replications.
    #'   Default \code{100}.
    #'
    #' @param betaT Numeric scalar or vector.  True treatment effect added to treated
    #'   subjects' outcomes.  The scale is response-type specific: additive for
    #'   \code{continuous}, \code{proportion}, and \code{ordinal}; on the logit
    #'   scale for \code{incidence}; log-multiplicative for \code{count} and
    #'   \code{survival}.  Default \code{1}.  Set \code{betaT = 0} to check
    #'   type-I error.
    #'
    #' @param alpha Numeric in \eqn{(0,1)}.  Significance level used for all
    #'   confidence intervals and for computing power (\eqn{p < \alpha}).
    #'   Default \code{0.05}.
    #'
    #' @param B_boot Positive integer.  Bootstrap resamples per CI / p-value
    #'   call.  Default \code{201}.
    #'
    #' @param r_rand Positive integer.  Randomisation draws per rand p-value
    #'   call, and per bisection step of the rand CI.  Default \code{201}.
    #'
    #' @param pval_epsilon Numeric.  Bisection convergence tolerance for
    #'   randomisation-based CIs (\code{compute_rand_confidence_interval}).
    #'   Default \code{0.02}.
    #'
    #' @param sd_noise Numeric \eqn{> 0}. Standard deviation of independent
    #'   Gaussian noise added to each subject's outcome.  Default \code{1}.
    #'
    #' @param n_ordinal_levels Positive integer. Number of ordinal categories when
    #'   \code{response_type = "ordinal"}. Default \code{4L}.
    #'
    #' @param proportion_epsilon Numeric scalar. Small value added to proportion
    #'   base responses to avoid 0 and 1. Default \code{1e-6}.
    #'
    #' @param phi_proportion Positive numeric scalar. Precision parameter for
    #'   beta-distributed observed proportion outcomes. The beta mean is
    #'   \code{y_linear_model[i] + betaT * w[i]}. Default \code{100}.
    #'
    #' @param k_survival Positive numeric scalar. Scale parameter passed to
    #'   the Weibull draw for observed survival outcomes. Default \code{2}.
    #'
    #' @param incidence_clamp Numeric scalar in \eqn{(0, 0.5)}. Clamp applied
    #'   to the Bernoulli probability for observed incidence outcomes.
    #'   Default \code{1e-9}.
    #'
    #' @param proportion_clamp Numeric scalar in \eqn{(0, 0.5)}. Clamp applied
    #'   to the beta mean for observed proportion outcomes. Default \code{1e-9}.
    #'
    #' @param count_clamp Positive numeric scalar. Minimum Poisson mean for
    #'   observed count outcomes. Default \code{1e-9}.
    #'
    #' @param survival_clamp Positive numeric scalar. Minimum Weibull shape for
    #'   observed survival outcomes. Default \code{1e-9}.
    #'
    #' @param survival_min_time Numeric scalar. Minimum survival time and shift
    #'   for base responses. Default \code{0.1}.
    #'
    #' @param count_min_rate Integer scalar. Minimum baseline rate for count
    #'   responses. Default \code{0L}.
    #'
    #' @param count_shift Numeric scalar. Constant added to counts after
    #'   zero-centering for base responses. Default \code{0}.
    #'
    #' @param norm_sq_beta_vec Positive numeric scalar. The desired squared
    #'   Euclidean norm of the latent linear coefficient vector \eqn{\beta}.
    #'   The generated vector is scaled to match this norm. Default \code{1}.
    #'
    #' @param X_mat Numeric matrix of dimensions \code{n x p}, or \code{NULL} (default).
    #'   If provided, these fixed covariates are used for every replication.
    #'   In this case, \code{cov_draw_method} must be \code{NULL}.
    #'
    #' @param seed Integer or \code{NULL} (default).  Random seed for the
    #'   entire simulation run.
    #'
    #' @param cov_draw_method A function used to draw \code{n * p} i.i.d. covariate
    #'   values for every replication. The function must accept the total number
    #'   of values as its first argument, followed by arguments in
    #'   \code{cov_draw_method_args}.  Default \code{stats::rnorm}.
    #'   Must be \code{NULL} when \code{X_mat} is supplied.
    #'
    #' @param cov_draw_method_args Named list of additional arguments forwarded to
    #'   \code{cov_draw_method} beyond the sample-size first argument. Default is
    #'   \code{list(mean = 0, sd = 1)}.
    #'
    #' @param random_X_draws Logical. If \code{TRUE} (default), a new set of
    #'   covariates is drawn for every single replication. If \code{FALSE},
    #'   one set is drawn per \code{(n, p)} cell and shared across its
    #'   replications.
    #'
    #' @param prob_censoring Numeric in \eqn{[0,1]}.  Per-subject independent
    #'   censoring probability; applied only when
    #'   \code{response_type = "survival"}.  Default \code{0.25}.
    #'
    #' @param custom_replication_data_generator Optional function for custom
    #'   replication data. When supplied, it is called as
    #'   \code{fn(state, rep)} and must return a list containing at least
    #'   \code{X} and \code{y_linear_model}. Any additional fields in the
    #'   returned list (e.g. latent frailty draws) are passed forward as
    #'   \code{rep_data} to \code{custom_apply_treatment_and_noise} and the
    #'   function built by \code{make_estimand_fn}.
    #'
    #' @param custom_apply_treatment_and_noise Optional function for custom
    #'   response generation. Signature: \code{fn(y_linear_model, w, rep_data, state)}.
    #'   \code{w} is in \{-1, +1\} format; convert with \code{(w+1)/2} for \{0,1\}
    #'   semantics. \code{rep_data} is the full list returned by
    #'   \code{custom_replication_data_generator} (or \code{NULL} for the standard
    #'   path). Must return a list with components \code{y} and \code{dead}.
    #'   Three-argument functions \code{fn(y_linear_model, w, state)} are still
    #'   accepted for backwards compatibility.
    #'
    #' @param make_estimand_fn Optional factory function for a custom true
    #'   estimand. Signature: \code{fn(beta_T)}, returning a function with
    #'   signature \code{fn(y_linear_model, X, w, rep_data, state)}. Called once
    #'   per grid cell with that cell's \code{beta_T} so the returned estimand
    #'   function is always tied to the right effect size (important when
    #'   \code{betaT} is a vector of multiple values). The returned function is
    #'   invoked once per design class per replication \emph{after} the design
    #'   completes, so \code{w} (in \{-1, +1\} format) and \code{X} reflect the
    #'   realized assignment. Must return a numeric scalar. When supplied, its
    #'   return value is used as the ground truth for \emph{all} inference classes
    #'   (overriding the \code{is_mean_diff} gate). Three-argument functions
    #'   \code{fn(y_linear_model, state)} are still accepted for backwards
    #'   compatibility (they will not receive \code{X} or \code{w}). Default
    #'   \code{NULL} uses \code{beta_T} directly as the ground truth.
    #'
    #' @param dgp_params Optional named list of DGP configuration values (e.g.
    #'   \code{list(frailty_dist = "gamma", censoring_rate = 0.8)}). Injected into
    #'   \code{state} as \code{state\$dgp_params} and accessible in all three
    #'   custom-DGP hooks. Recommended over using closures to pass DGP parameters.
    #'
    #' @param custom_dgp Optional function for a fully custom DGP. Signature:
    #'   \code{fn(n, p, rep, state)} returning a list with components
    #'   \code{X} (data.frame, \code{n} x \code{p}),
    #'   \code{w} (integer vector in \{0, 1\}, length \code{n}),
    #'   \code{y} (numeric, length \code{n}),
    #'   \code{dead} (integer \{0,1\} or \code{NULL} for non-survival),
    #'   \code{true_estimand} (numeric scalar, optional).
    #'   When supplied, the design class acts as a data container only; it does
    #'   \emph{not} run its own randomization or matching. Requires a fixed design
    #'   class (not \code{DesignSeqOneByOne} variants). Cannot be combined with
    #'   \code{custom_replication_data_generator} or
    #'   \code{custom_apply_treatment_and_noise}.
    #'
    #' @param verbose Logical.  If \code{TRUE}, prints a message for every
    #'   replication and for every \code{(design, inference)} pair that is
    #'   skipped due to an error.  Default \code{TRUE}.
    #'
    #' @param keep_all_intermediate_data Logical. If \code{TRUE}, the framework
    #'   saves the instantiated design and inference objects for every replication.
    #'   These can be retrieved after the run using \code{$get_all_intermediate_data()}.
    #'   Warning: this can consume a lot of memory for many replications.
    #'   Default \code{FALSE}.
    #'
    #' @param turn_off_asserts_for_speed Logical. If \code{TRUE} (default),
    #'   all \pkg{checkmate} assertions across the package are globally
    #'   disabled during the simulation run to improve performance.
    #'
    #' @param num_cores Positive integer.  Number of worker processes for
    #'   parallel execution of Monte Carlo replications.  Note that when
    #'   \code{num_cores > 1}, parallelization *within* individual inference
    #'   routines (e.g. bootstrap, randomization) is automatically disabled
    #'   to prevent thread oversubscription.
    #'
    #'   \strong{Unix/Linux (recommended):} A \code{makeForkCluster} pool is
    #'   created once at \code{run()} start.  Workers inherit all pre-generated
    #'   design and SE caches via copy-on-write with zero serialization overhead.
    #'   Parallelism operates at the \emph{replication} level: \code{num_cores}
    #'   replications run simultaneously, each executing all DGP cells serially.
    #'   This eliminates the per-batch dispatch overhead that would arise from
    #'   cycling through cells within every replication, and keeps all cores
    #'   fully subscribed regardless of the number of DGP cells.  For best
    #'   performance, run on a Unix/Linux machine and set \code{num_cores} to
    #'   the number of physical cores available.
    #'
    #'   \strong{Non-Unix (Windows/macOS):} \pkg{mirai} is used when available
    #'   (parallelism at the DGP-cell level within each replication), otherwise
    #'   execution falls back to serial with a warning.  Default \code{1}.
    #'
    #' @param results_filename Character scalar. The filename for the results
    #'   file. Supported extensions are \code{.csv} and \code{.csv.bz2}.
    #'   Default \code{"simulation_framework_results.csv.bz2"}.
    #'
    #' @param continue_from_last_result_row Logical. If \code{TRUE} (default),
    #'   the framework loads existing results from \code{results_filename} and
    #'   skips previously completed replications.
    #'
    #' @param reuse_cache Logical. If \code{TRUE} (default), expensive
    #'   pre-generated design / SE cache objects are loaded from disk when
    #'   available. If \code{FALSE}, these cache objects are regenerated from
    #'   scratch, but each regenerated object is still saved to disk for later
    #'   restarts.
    #'
    #' @param stop_on_error Logical. If \code{TRUE} (default), any error raised
    #'   during a simulation path aborts the run immediately. If \code{FALSE},
    #'   the framework records the error, skips the failing path, and continues
    #'   with the remaining replications / design / inference combinations. Use
    #'   \code{$get_errors()} after \code{$run()} to inspect the captured
    #'   errors.
    #'
    #' @param save_to_disk_every_n_rep Positive integer. Results are flushed to
    #'   the on-disk staging file only once every this many replications, and
    #'   always after the final replication. Larger values reduce disk I/O
    #'   overhead at the cost of losing more progress if the run is interrupted.
    #'   Default \code{25L}.
    #'
    #' @param save_model_control_fits Logical. If \code{TRUE} (default), after the
    #'   design/SE cache is built, saves the per-subject model-implied potential
    #'   outcomes under treatment and control as CSV files in a subfolder named
    #'   \code{<stem>_response_values/} (where \code{<stem>} is \code{results_filename}
    #'   with its \code{.csv}/\code{.csv.bz2} extension stripped) next to
    #'   \code{results_filename}.  One file is written per unique
    #'   \code{(response_type, cond_exp_func_model, n, p, betaT)} cell.  Only
    #'   meaningful when \code{random_X_draws = FALSE}; silently skipped otherwise.
    #'   Column names depend on \code{response_type}:
    #'   \describe{
    #'     \item{\code{"continuous"}, \code{"survival"}}{columns \code{yt} and \code{yc}}
    #'     \item{\code{"incidence"}, \code{"proportion"}}{columns \code{pt} and \code{pc}}
    #'     \item{\code{"count"}}{columns \code{rt} and \code{rc}}
    #'   }
    #'   Default \code{TRUE}.
    #'
    #' @param inference_types_and_params \code{NULL} (default) or a named list
    #'   from inference type to a named list of arguments for that type's function
    #'   invocation.  The list names control which inference outputs are computed.
    #'   Valid names are \code{"asymp_ci"}, \code{"asymp_pval"},
    #'   \code{"exact_ci"}, \code{"exact_pval"}, \code{"boot_ci"},
    #'   \code{"boot_pval"}, \code{"rand_ci"}, and \code{"rand_pval"}.
    #'   Each value must be a named list whose names are accepted by the
    #'   corresponding inference function.  \code{NULL} runs all eight types with
    #'   default invocation arguments.
    #'   Example:
    #'   \preformatted{inference_types_and_params = list(
    #'   asymp_pval = list(delta = 0),
    #'   boot_ci    = list(B = 99, type = "perc"),
    #'   rand_pval  = list(r = 999, transform_responses = TRUE)
    #' )}
    #'   When no \code{*_ci} type is requested, \code{coverage} is omitted from
    #'   \code{SimulationFrameworkReport$summarize()}.  When no \code{*_pval}
    #'   type is requested, \code{power} is omitted.
    initialize = function(
      response_type,
      design_classes_and_params = NULL,
      inference_classes_and_params = NULL,
      n                     = 100L,
      p                     = 5L,
      cond_exp_func_model   = "linear",
      Nrep                  = 100L,
      betaT                 = 1,
      alpha                 = 0.05,
      B_boot                = 201L,
      r_rand                = 201L,
      pval_epsilon          = 0.02,
      sd_noise              = 1,
      n_ordinal_levels      = 4L,
      proportion_epsilon    = 1e-6,
      phi_proportion        = 100,
      k_survival            = 2,
      incidence_clamp       = 1e-9,
      proportion_clamp      = 1e-9,
      count_clamp           = 1e-9,
      survival_clamp        = 1e-9,
      survival_min_time     = 0.1,
      count_min_rate        = 0L,
      count_shift           = 0,
      norm_sq_beta_vec      = 1,
      X_mat                 = NULL,
      num_cores             = 1L,
      seed                  = NULL,
      cov_draw_method       = stats::rnorm,
      cov_draw_method_args  = list(mean = 0, sd = 1),
      random_X_draws        = TRUE,
      prob_censoring        = 0.25,
      custom_replication_data_generator = NULL,
      custom_apply_treatment_and_noise = NULL,
      make_estimand_fn = NULL,
      dgp_params = list(),
      custom_dgp = NULL,
      verbose                    = TRUE,
      keep_all_intermediate_data = FALSE,
      turn_off_asserts_for_speed = TRUE,
      inference_types_and_params = NULL,
      results_filename      = "simulation_framework_results.csv.bz2",
      continue_from_last_result_row = TRUE,
      reuse_cache           = TRUE,
      stop_on_error         = TRUE,
      save_to_disk_every_n_rep = 25L,
      save_model_control_fits = TRUE
    ) {
      valid_rt = c("continuous", "incidence", "proportion",
                   "count", "survival", "ordinal")
      if (any(!response_type %in% valid_rt))
        stop("response_type must be one of: ", paste(valid_rt, collapse = ", "))
      n_values = unique(as.integer(n))
      p_values = unique(as.integer(p))
      betaT_values = unique(as.numeric(betaT))
      cond_exp_func_model_values = unique(as.character(cond_exp_func_model))
      if (length(n_values) == 0L || any(!is.finite(n_values)) || any(n_values <= 1L))
        stop("n must contain finite integers greater than 1")
      if (length(p_values) == 0L || any(!is.finite(p_values)) || any(p_values < 1L))
        stop("p must contain finite positive integers")
      if (length(betaT_values) == 0L || any(!is.finite(betaT_values)))
        stop("betaT must contain finite numeric values")
      if (length(cond_exp_func_model_values) == 0L ||
          any(!cond_exp_func_model_values %in% c("linear", "nonlinear")))
        stop("cond_exp_func_model must contain only 'linear' and/or 'nonlinear'")
      if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)))
        stop("seed must be NULL or one finite numeric value")
      if (!is.null(X_mat) && (length(n_values) > 1L || length(p_values) > 1L))
        stop("X_mat can only be used when n and p are scalar")
      if (!isTRUE(random_X_draws) && is.null(seed))
        stop("random_X_draws = FALSE requires seed to be non-NULL")
      if (!is.character(results_filename) || length(results_filename) != 1L || is.na(results_filename))
        stop("results_filename must be a single non-missing character string")
      results_file_format = private$.results_file_format(results_filename)
      if (is.na(results_file_format)) {
        stop("results_filename must end in either '.csv' or '.csv.bz2'")
      }
	      if (!is.finite(phi_proportion) || phi_proportion <= 0)
	        stop("phi_proportion must be finite and > 0")
	      if (!is.finite(k_survival) || k_survival <= 0)
	        stop("k_survival must be finite and > 0")
	      if (!is.finite(incidence_clamp) || incidence_clamp <= 0 || incidence_clamp >= 0.5)
	        stop("incidence_clamp must be finite and in (0, 0.5)")
	      if (!is.finite(proportion_clamp) || proportion_clamp <= 0 || proportion_clamp >= 0.5)
	        stop("proportion_clamp must be finite and in (0, 0.5)")
	      if (!is.finite(count_clamp) || count_clamp <= 0)
	        stop("count_clamp must be finite and > 0")
	      if (!is.finite(survival_clamp) || survival_clamp <= 0)
	        stop("survival_clamp must be finite and > 0")
      valid_inf_types = c("asymp_ci", "asymp_pval", "exact_ci", "exact_pval",
                          "boot_ci",  "boot_pval",  "rand_ci",  "rand_pval")
      inf_type_spec = private$.parse_inference_types_and_params(
        inference_types_and_params,
        valid_inf_types
      )
      inf_types = names(inf_type_spec)
      private$response_type_values = unique(as.character(response_type))
      private$n_values         = n_values
      private$p_values         = p_values
      private$cond_exp_func_model_values = cond_exp_func_model_values
      private$Nrep             = as.integer(Nrep)
      private$betaT_values     = betaT_values
      private$alpha            = alpha
      private$B_boot           = as.integer(B_boot)
      private$r_rand           = as.integer(r_rand)
      private$pval_epsilon     = pval_epsilon
      private$sd_noise             = sd_noise
      private$n_ordinal_levels     = as.integer(n_ordinal_levels)
      private$proportion_epsilon    = proportion_epsilon
      private$phi_proportion        = phi_proportion
      private$k_survival            = k_survival
      private$incidence_clamp       = incidence_clamp
      private$proportion_clamp      = proportion_clamp
      private$count_clamp           = count_clamp
      private$survival_clamp        = survival_clamp
      private$survival_min_time     = survival_min_time
      private$count_min_rate        = as.integer(count_min_rate)
      private$count_shift           = count_shift
      private$norm_sq_beta_vec     = norm_sq_beta_vec
      private$X_mat                = X_mat
      private$num_cores            = as.integer(num_cores)
      private$seed                 = if (is.null(seed)) NULL else as.integer(seed)
      private$cov_draw_method      = if (is.null(X_mat)) cov_draw_method else NULL
      private$cov_draw_method_args = cov_draw_method_args
      private$random_X_draws       = random_X_draws
      private$prob_censoring       = prob_censoring
      private$custom_replication_data_generator = custom_replication_data_generator
      private$custom_apply_treatment_and_noise = custom_apply_treatment_and_noise
      private$make_estimand_fn = make_estimand_fn
      if (!is.null(custom_dgp)) {
        if (!is.null(custom_replication_data_generator) || !is.null(custom_apply_treatment_and_noise))
          stop("custom_dgp cannot be combined with custom_replication_data_generator or custom_apply_treatment_and_noise")
        if (!is.function(custom_dgp))
          stop("custom_dgp must be a function")
        if (all(betaT_values != 0))
          warning("custom_dgp is set; betaT is ignored as the true estimand comes from the DGP function")
      }
      if (!isTRUE(random_X_draws) && !is.null(custom_replication_data_generator))
        warning("custom_replication_data_generator overrides random_X_draws=FALSE; a new X will be drawn every replication")
      private$dgp_params = if (is.null(dgp_params)) list() else dgp_params
      private$custom_dgp = custom_dgp
      private$verbose                    = verbose
      private$turn_off_asserts_for_speed = turn_off_asserts_for_speed
      
      if (as.integer(num_cores) > 1L && keep_all_intermediate_data) {
        stop("Multithreading (num_cores > 1) is incompatible with 'keep_all_intermediate_data = TRUE'.")
      }
      
      private$keep_all_intermediate_data = keep_all_intermediate_data
      private$results_filename     = results_filename
      private$continue_from_last_result_row = continue_from_last_result_row
      private$reuse_cache          = isTRUE(reuse_cache)
      private$stop_on_error        = isTRUE(stop_on_error)
      checkmate::assertCount(save_to_disk_every_n_rep, positive = TRUE)
      private$save_to_disk_every_n_rep = as.integer(save_to_disk_every_n_rep)
      private$save_model_control_fits = isTRUE(save_model_control_fits)
      private$inf_types        = inf_types
      private$inference_type_params = inf_type_spec
      private$param_grid       = private$.build_param_grid(
        n_values,
        p_values,
        betaT_values,
        cond_exp_func_model_values,
        private$response_type_values
      )
      private$current_n        = private$param_grid$n[[1L]]
      private$current_p        = private$param_grid$p[[1L]]
      private$current_betaT    = private$param_grid$betaT[[1L]]
      private$current_cond_exp_func_model = private$param_grid$cond_exp_func_model[[1L]]
      private$current_response_type = private$param_grid$response_type[[1L]]
      design_spec = private$.parse_design_classes_and_params(
        design_classes_and_params,
        parent.frame()
      )
      private$design_classes = design_spec$classes
      private$design_params  = design_spec$params
      inference_spec = private$.parse_inference_classes_and_params(
        inference_classes_and_params,
        parent.frame()
      )
      private$inference_classes = inference_spec$classes
      private$inference_constructor_params = inference_spec$params
      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)
      private$design_labels    = private$.compute_design_labels()
      private$inference_labels = private$.compute_inference_labels()
      private$raw_results = list()
      private$error_log   = list()
      private$has_run     = FALSE
    },
    # ── run() ─────────────────────────────────────────────────────────────────
    #' @description Execute the simulation replications.
    #'
    #' @return The \code{SimulationFramework} object itself (invisibly).
    run = function() {
      if (!is.null(private$seed)) {
        had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (had_seed) {
          old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
          on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
        } else {
          on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
        }
        set.seed(private$seed)
      }
      # Disable assertions for the duration of the simulation for speed
      if (private$turn_off_asserts_for_speed){
        toggle_asserts(FALSE)
        on.exit(toggle_asserts(TRUE), add = TRUE)        
      }
      
      # ── Parallelism management ─────────────────────────────────────────────
      # Save state to restore on exit
      ns = asNamespace("EDI")
      prev_threads = getOption(".edi_last_set_threads")
      if (is.null(prev_threads)) prev_threads = 1L
      prev_global_cores = get_num_cores()
      prev_global_mirai_cores = get_global_mirai_cores()
      prev_override = ns$edi_env$num_cores_override
      
      num_cores = private$num_cores
      set_package_threads(num_cores)
      
      # If the simulation is running its own parallel loop, we force 
      # nested core budget to 1 to prevent N*M explosion.
      if (num_cores > 1L && prev_global_cores > 1L) {
        set_num_cores(1L)
      }
      
      on.exit({
        set_package_threads(prev_threads)
        if (get_num_cores() != prev_global_cores) {
           set_num_cores(prev_global_cores, force_mirai = !is.null(prev_global_mirai_cores))
        }
        assign("num_cores_override", prev_override, envir = ns$edi_env)
      }, add = TRUE)
      has_pregen_designs = any(vapply(private$design_classes, function(dc)
        !is.null(dc$public_methods$supports_batch_w_pregeneration), logical(1L)))
      has_cmh_inference = "InferenceIncidCMH" %in% private$inference_labels
      # force_serial was previously tied to has_pregen_designs to avoid JVM fork
      # issues; Java/GED is gone so parallelism is always safe.
      force_serial = FALSE
      use_fork_cluster  = !force_serial && isTRUE(num_cores > 1L) && .Platform$OS.type == "unix"
      use_mirai_backend = !force_serial && isTRUE(num_cores > 1L) && !use_fork_cluster
      # Rep/DGP bars always show 100% under fork parallelism because workers
      # return complete reps; suppress them to avoid misleading output.
      private$n_progress_lines = if (use_fork_cluster) 3L else 5L
      cl_cache = NULL  # fork cluster used only for cache building phase
      cl_rep   = NULL  # fork cluster for rep execution, created AFTER cell states are populated
      if (use_fork_cluster) {
        cl_cache = parallel::makeForkCluster(num_cores)
      }
      on.exit({
        if (!is.null(cl_cache)) try(parallel::stopCluster(cl_cache), silent = TRUE)
        if (!is.null(cl_rep))   try(parallel::stopCluster(cl_rep),   silent = TRUE)
      }, add = TRUE)
      if (use_mirai_backend) {
        if (!check_package_installed("mirai")) {
          use_mirai_backend = FALSE
          if (isTRUE(private$verbose)) {
            private$.message_stderr(
              "Warning: Parallelism (num_cores > 1) requires the 'mirai' package. Serial execution used for this run.\n"
            )
          }
        } else {
          private$.ensure_mirai_daemons(num_cores)
          on.exit({
            if (!is.null(prev_global_mirai_cores)) {
              private$.ensure_mirai_daemons(prev_global_mirai_cores)
            } else {
              tryCatch(mirai::daemons(0), error = function(e) invisible(NULL))
            }
          }, add = TRUE)
        }
      }
      # Handle cleanup/setup
      elapsed_f = private$.elapsed_file()
      if (!isTRUE(private$continue_from_last_result_row)) {
        if (file.exists(private$results_filename)) unlink(private$results_filename)
        private$.cleanup_results_staging_file()
        cache_dir = private$.simulation_cache_dir()
        if (dir.exists(cache_dir)) unlink(cache_dir, recursive = TRUE)
        if (file.exists(elapsed_f)) unlink(elapsed_f)
        private$prior_elapsed_secs = 0
      } else {
        prior = suppressWarnings(tryCatch(readRDS(elapsed_f), error = function(e) 0))
        private$prior_elapsed_secs = if (is.numeric(prior) && length(prior) >= 1L && is.finite(prior[[1L]])) prior[[1L]] else 0
      }
      
      private$simulation_start_time = as.numeric(Sys.time())
      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)
      n_met = length(private$inf_types)
      n_cells = nrow(private$param_grid)
      # Ensure staging file is ready if we are using bz2 and continuing
      if (isTRUE(private$continue_from_last_result_row)) {
        private$.ensure_staging_file_exists()
      }
      existing_results = private$.load_existing_results()
      n_existing = nrow(existing_results)
      
      if (isTRUE(private$verbose) && n_existing > 0L) {
        private$.message_stderr(sprintf("%s existing results loaded\n", format(n_existing, big.mark = ",")))
      }
      # Pre-allocate results list for data.table batches (one per chunk or serial replication)
      private$raw_results = vector("list", 1L + private$Nrep * n_cells)
      private$results_idx = 0L
      # Use C++ ResultKeyStore for O(1) key lookups without R string interning overhead
      init_result_key_store_cpp(n_existing + private$Nrep * n_cells)
      if (n_existing > 0L) {
        # Optimization: Store the entire data.table as the first element instead of row-by-row lists
        private$raw_results[[1L]] = existing_results
        private$results_idx = 1L
        
        # Extract columns once as vectors (fast pointer-based extraction in data.table)
        rt_v = existing_results$response_type
        cm_v = existing_results$cond_exp_func_model
        n_v  = existing_results$n
        p_v  = existing_results$p
        bt_v = existing_results$betaT
        rp_v = existing_results$rep
        ds_v = existing_results$design
        if_v = existing_results$inference
        it_v = existing_results$inference_type
        # Chunking allows showing progress and avoids one massive allocation.
        chunk_size_indexing = 200000L
        indices = seq(1L, n_existing, by = chunk_size_indexing)
        
        for (i in seq_along(indices)) {
          if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar("indexing existing results", (i - 1L) / length(indices))
          start = indices[i]
          end = min(n_existing, start + chunk_size_indexing - 1L)
          
          add_to_result_key_store_cpp(
            rt_v, cm_v, n_v, p_v, bt_v, rp_v, ds_v, if_v, it_v,
            start, end
          )
        }
        if (isTRUE(private$verbose)) {
          private$.draw_labeled_progress_bar("indexing existing results", 1)
          cat("\n", file = stderr())
        }
      }
      private$all_intermediate_data = vector("list", private$Nrep * n_cells)
      private$valid_combos          = list()
      private$seen_combo_keys      = character(0L)
      private$exact_warned_classes = character(0L)
      private$error_log            = list()
      private$pending_file_rows    = list()
      if (isTRUE(private$verbose)) {
        message(sprintf(
          "simulations: CEF_mod=%s  n=%s  p=%s  Nrep=%d  betaT=%s designs=%d inferences=%d num_cores=%d",
          private$.format_values(private$cond_exp_func_model_values),
          private$.format_values(private$n_values),
          private$.format_values(private$p_values),
          private$Nrep,
          private$.format_values(private$betaT_values), 
          n_des, 
          n_inf,
          num_cores
        ))
      }
      log_interval = max(1L, private$Nrep %/% 10L)
      shared_X_draws = list()
      if (!isTRUE(private$random_X_draws)) {
        np_grid = unique(private$param_grid[, .(n, p)])
        n_np = nrow(np_grid)
        
        # Optimization: convert X_mat once outside the loop
        X_mat_matrix = if (!is.null(private$X_mat)) as.matrix(private$X_mat) else NULL
        
        for (np_idx in seq_len(n_np)) {
          n_i = np_grid$n[[np_idx]]
          p_i = np_grid$p[[np_idx]]
          X_i = if (is.null(X_mat_matrix)) {
            m = matrix(
              do.call(private$cov_draw_method, c(list(n_i * p_i), private$cov_draw_method_args)),
              nrow = n_i,
              ncol = p_i
            )
            colnames(m) = paste0("x", seq_len(p_i))
            m
          } else {
            X_mat_matrix # Names should already be there or handled by generate_covariate_dataset
          }
          shared_X_draws[[paste(n_i, p_i, sep = "|")]] = X_i
        }
      }
      # Pre-calculate planned combos and cache cell metadata
      had_seed_plan = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed_plan) old_seed_plan = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      
      planned_combos_list = vector("list", n_cells)
      cell_reps_to_run    = vector("list", n_cells)
      validity_cache      = list()
      cell_tasks_count    = integer(n_cells)
      total_existing_planned_tasks = 0L
      
      for (cell_idx in seq_len(n_cells)) {
        n_i  = private$param_grid$n[[cell_idx]]
        p_i  = private$param_grid$p[[cell_idx]]
        dt_i = private$param_grid$cond_exp_func_model[[cell_idx]]
        bt_i = private$param_grid$betaT[[cell_idx]]
        rt_i = private$param_grid$response_type[[cell_idx]]
        # Show more granular progress during plan pre-calculation
        if (isTRUE(private$verbose)) {
          label = sprintf("pre-calculating plan (cell %d/%d: n=%d, p=%d, %s, %s)",
                          cell_idx, n_cells, n_i, p_i, dt_i, rt_i)
          private$.draw_labeled_progress_bar(label, (cell_idx - 1) / n_cells)
        }
        cache_key = paste(rt_i, n_i, p_i, dt_i, sep = "|")
        if (is.null(validity_cache[[cache_key]])) {
           private$current_n = n_i; private$current_p = p_i; private$current_betaT = bt_i; private$current_cond_exp_func_model = dt_i; private$current_response_type = rt_i
           rep_data = if (isTRUE(private$random_X_draws)) private$.generate_data() else private$.generate_data_from_X(shared_X_draws[[paste(n_i, p_i, sep = "|")]])
           validity_cache[[cache_key]] = private$.build_valid_combos_for_current_cell(rep_data)
        }
        
        cell_combos = lapply(validity_cache[[cache_key]], function(c) {
           c$betaT = bt_i; c$n = n_i; c$p = p_i; c$cond_exp_func_model = dt_i; c$response_type = rt_i; c
        })
        planned_combos_list[[cell_idx]] = cell_combos
        n_combos = length(cell_combos)
        cell_tasks_count[cell_idx] = n_combos
        
        # Identify missing replications for this specific cell
        reps_to_run = seq_len(private$Nrep)
        if (n_existing > 0L && n_combos > 0L) {
          # Use batches for safety with extremely large Nrep
          rep_batch_size = max(100L, 500000L %/% n_combos)
          rep_batches = split(seq_len(private$Nrep), (seq_len(private$Nrep) - 1) %/% rep_batch_size)
          
          c_des  = vapply(cell_combos, `[[`, "", "design")
          c_inf  = vapply(cell_combos, `[[`, "", "inference")
          c_type = vapply(cell_combos, `[[`, "", "inference_type")
          
          cell_req = rep(TRUE, private$Nrep)
          for (b_reps in rep_batches) {
            nb = length(b_reps)
            all_exists = check_in_result_key_store_cpp(
              rep(rt_i, nb * n_combos),
              rep(dt_i, nb * n_combos),
              rep(as.integer(n_i), nb * n_combos),
              rep(as.integer(p_i), nb * n_combos),
              rep(as.numeric(bt_i), nb * n_combos),
              rep(as.integer(b_reps), each = n_combos),
              rep(c_des, nb),
              rep(c_inf, nb),
              rep(c_type, nb)
            )
            total_existing_planned_tasks = total_existing_planned_tasks + sum(all_exists)
            exists_mat = matrix(all_exists, nrow = n_combos)
            cell_req[b_reps] = colSums(exists_mat) < n_combos
          }
          if (!private$keep_all_intermediate_data) {
             reps_to_run = which(cell_req)
          }
        }
        cell_reps_to_run[[cell_idx]] = reps_to_run
      }
      if (isTRUE(private$verbose)) {
        private$.draw_labeled_progress_bar("pre-calculating simulation plan", 1)
        cat("\n", file = stderr())
      }
      if (had_seed_plan) assign(".Random.seed", old_seed_plan, envir = .GlobalEnv)
      
      planned_combos = unlist(planned_combos_list, recursive = FALSE)
      private$valid_combos = planned_combos
      private$progress_total = length(planned_combos) * private$Nrep
      
      private$progress_count = total_existing_planned_tasks
      private$initial_progress_count = total_existing_planned_tasks
      private$progress_log_interval = max(1L, private$progress_total %/% 20L)
      private$progress_bar = NULL; private$use_progress_bar = FALSE
      private$total_cells = n_cells; private$last_progress_draw_time = 0
      private$current_task_label = "Des/Inf"

      if (isTRUE(private$verbose)) private$.print_plan_summary(planned_combos_list)
      # ── Save per-subject model-implied response values (if requested) ────────
      if (isTRUE(private$save_model_control_fits) && !isTRUE(private$random_X_draws)) {
        rv_dir = {
          path = private$.results_output_path()
          stem = sub("\\.csv(\\.bz2)?$", "", basename(path), ignore.case = TRUE)
          file.path(dirname(path), paste0(stem, "_response_values"))
        }
        if (!dir.exists(rv_dir)) dir.create(rv_dir, recursive = TRUE, showWarnings = FALSE)
        cell_grid = unique(private$param_grid[, .(response_type, cond_exp_func_model, n, p, betaT)])
        for (gi in seq_len(nrow(cell_grid))) {
          rt_g  = cell_grid$response_type[[gi]]
          dt_g  = cell_grid$cond_exp_func_model[[gi]]
          n_g   = cell_grid$n[[gi]]
          p_g   = cell_grid$p[[gi]]
          bt_g  = cell_grid$betaT[[gi]]
          X_g   = shared_X_draws[[paste(n_g, p_g, sep = "|")]]
          if (is.null(X_g)) next
          dat_g = generate_covariate_dataset(
            n = n_g, p = p_g,
            cond_exp_func_model  = dt_g,
            norm_sq_beta_vec     = private$norm_sq_beta_vec,
            X_mat                = X_g,
            cov_draw_method      = NULL,
            cov_draw_method_args = private$cov_draw_method_args
          )
          ylm_g = as.numeric(dat_g$y_cont - mean(dat_g$y_cont))
          rv_dt = switch(rt_g,
            continuous = data.table::data.table(
              yt = ylm_g + bt_g,
              yc = ylm_g
            ),
            survival = data.table::data.table(
              yt = exp(ylm_g + bt_g),
              yc = exp(ylm_g)
            ),
            incidence = data.table::data.table(
              pt = stats::plogis(ylm_g + bt_g),
              pc = stats::plogis(ylm_g)
            ),
            proportion = data.table::data.table(
              pt = stats::plogis(ylm_g + bt_g),
              pc = stats::plogis(ylm_g)
            ),
            count = data.table::data.table(
              rt = exp(ylm_g + bt_g),
              rc = exp(ylm_g)
            ),
            NULL
          )
          if (is.null(rv_dt)) next
          bt_str   = gsub("[^A-Za-z0-9_-]", "_", as.character(bt_g))
          rv_fname = sprintf("rv_%s_%s_n%d_p%d_betaT%s.csv", rt_g, dt_g, n_g, p_g, bt_str)
          data.table::fwrite(rv_dt, file.path(rv_dir, rv_fname))
        }
      }
      # Early exit if everything is already done (before drawing the main progress bar)
      if (private$progress_count >= private$progress_total) {
        if (isTRUE(private$verbose)) message("Simulation already complete.")
        private$.finish_run()
        return(invisible(self))
      }
      # Pre-build all per-cell state objects and run-required vectors so the
      # outer rep loop can iterate cells without rebuilding them each rep.
      all_cell_states       = vector("list", n_cells)
      all_run_required_v    = vector("list", n_cells)

      for (cell_idx in seq_len(n_cells)) {
        n_ci  = private$param_grid$n[[cell_idx]]
        p_ci  = private$param_grid$p[[cell_idx]]
        bt_ci = private$param_grid$betaT[[cell_idx]]
        dt_ci = private$param_grid$cond_exp_func_model[[cell_idx]]
        rt_ci = private$param_grid$response_type[[cell_idx]]

        reps_to_run_ci = cell_reps_to_run[[cell_idx]]
        run_req_ci = rep(FALSE, private$Nrep)
        run_req_ci[reps_to_run_ci] = TRUE
        all_run_required_v[[cell_idx]] = run_req_ci

        cell_key_prefix_ci = paste(rt_ci, dt_ci, n_ci, p_ci, bt_ci, sep = "|")
        cell_combos_here = planned_combos_list[[cell_idx]]
        cell_valid_combo_keys_ci = if (length(cell_combos_here) > 0L) {
          vapply(cell_combos_here, function(c)
            paste(c$design, c$inference, c$inference_type, sep = "|"), character(1L))
        } else {
          character(0L)
        }
        cs = list(
          n = n_ci, p = p_ci, betaT = bt_ci,
          cond_exp_func_model = dt_ci, norm_sq_beta_vec = private$norm_sq_beta_vec,
          response_type = rt_ci, random_X_draws = private$random_X_draws,
          shared_X = if (!isTRUE(private$random_X_draws)) shared_X_draws[[paste(n_ci, p_ci, sep = "|")]] else NULL,
          X_mat = private$X_mat, cov_draw_method = private$cov_draw_method, cov_draw_method_args = private$cov_draw_method_args,
          custom_replication_data_generator = private$custom_replication_data_generator,
          custom_apply_treatment_and_noise = private$custom_apply_treatment_and_noise,
          make_estimand_fn = private$make_estimand_fn,
          dgp_params = private$dgp_params,
          custom_dgp = private$custom_dgp,
          design_classes = private$design_classes, design_labels = private$design_labels, design_params = private$design_params,
          inference_classes = private$inference_classes, inference_labels = private$inference_labels, inference_ctor_params = private$inference_constructor_params,
          inf_types = private$inf_types, alpha = private$alpha, B_boot = private$B_boot, r_rand = private$r_rand,
          pval_epsilon = private$pval_epsilon, sd_noise = private$sd_noise, n_ordinal_levels = private$n_ordinal_levels,
          phi_proportion = private$phi_proportion, k_survival = private$k_survival,
          incidence_clamp = private$incidence_clamp, proportion_clamp = private$proportion_clamp, count_clamp = private$count_clamp,
          survival_clamp = private$survival_clamp, survival_min_time = private$survival_min_time, count_min_rate = private$count_min_rate,
          count_shift = private$count_shift, prob_censoring = private$prob_censoring,
          inference_type_params = private$inference_type_params,
          cell_key_prefix = cell_key_prefix_ci,
          valid_combo_keys = cell_valid_combo_keys_ci,
          stop_on_error = private$stop_on_error
        )
        all_cell_states[[cell_idx]] = cs
      }
      private$n_active_reps_total       = if (length(all_run_required_v) > 0L) sum(Reduce(`|`, all_run_required_v)) else 0L
      private$n_total_active_work_units = if (length(all_run_required_v) > 0L) sum(vapply(all_run_required_v, sum, 0L)) else 0L
      use_parallel_workers = !force_serial && num_cores > 1L
      num_cores_to_use = if (use_parallel_workers) num_cores else 1L

      # Pre-generate ALL w vectors for designs whose matching/clustering structure
      # is expensive but deterministic given X (e.g. BinaryMatch, OptimalBlocks),
      # and for non-blocking designs used with CMH inference (K = Nrep so the
      # SE estimate uses the full population of draws, no extra parameter needed).
      if ((has_pregen_designs || has_cmh_inference) && !isTRUE(private$random_X_draws)) {
        cache_jobs = list()
        for (ci in seq_len(n_cells)) {
          cs = all_cell_states[[ci]]
          if (is.null(cs$shared_X)) next
          active_reps_ci = which(all_run_required_v[[ci]])
          if (length(active_reps_ci) == 0L) next
          cell_combos_ci = planned_combos_list[[ci]]
          for (di in seq_along(private$design_classes)) {
            dc = private$design_classes[[di]]
            dl = private$design_labels[[di]]
            dp = private$design_params[[di]]
            is_pregen = !is.null(dc$public_methods$supports_batch_w_pregeneration)
            if (!is_pregen && !has_cmh_inference) next
            des_combos = Filter(function(co) co$design == dl, cell_combos_ci)
            if (length(des_combos) == 0L) next
            # Use active_reps_ci directly — avoids 312 expensive check_in_result_key_store_cpp
            # batches (18M C++ lookups) that caused a multi-second pause at startup.
            # active_reps_ci already filters to reps where any combo for this cell is
            # incomplete, which is a valid superset of per-design reps_needing.
            reps_needing_dl = active_reps_ci
            if (length(reps_needing_dl) == 0L) next
            # Stripped cs: only what _run_simulation_cache_job needs.
            # Avoids serializing cs$design_classes / cs$inference_classes (R6 generators)
            # when jobs are dispatched to fork cluster workers via clusterApply.
            cs_cache = list(n = cs$n, response_type = cs$response_type, shared_X = cs$shared_X)
            cache_file = private$.simulation_cache_file(cs, dl, dp, "design_w")
            cache_jobs[[length(cache_jobs) + 1L]] = list(
              cell_idx = ci,
              design_idx = di,
              cs = cs_cache,
              design_gen_class = dc$classname,
              design_label = dl,
              design_params = dp,
              cache_type = "design_w",
              is_pregen = is_pregen,
              reps_needing = reps_needing_dl,
              cache_file = cache_file,
              ready = isTRUE(private$reuse_cache) && file.exists(cache_file)
            )
          }
        }
        n_cache_jobs = length(cache_jobs)
        # Distribute the thread budget across workers: if there are fewer cache jobs than
        # cores, each worker gets the spare cores for its drawing step and BLAS calls.
        # When jobs >= cores every worker gets 1 to avoid N*M thread explosion.
        if (n_cache_jobs > 0L && use_parallel_workers) {
          cores_per_cache_job = max(1L, num_cores_to_use %/% n_cache_jobs)
          cache_jobs = lapply(cache_jobs, function(j) { j$num_cores_for_worker = cores_per_cache_job; j })
        }
        attach_cache_job = function(job) {
          cs = all_cell_states[[job$cell_idx]]
          obj = private$.load_simulation_cache_object(
            cs,
            private$design_labels[[job$design_idx]],
            private$design_params[[job$design_idx]],
            job$cache_type,
            reps_needing = job$reps_needing,
            restore_rng = FALSE,
            cache_file = job$cache_file
          )
          if (is.null(obj)) return(FALSE)
          dl = private$design_labels[[job$design_idx]]
          if (identical(job$cache_type, "design_w")) {
            if (is.list(obj) && isTRUE(obj$blocking_design)) {
              # sentinel for non-pregen blocking designs — CMH uses closed-form SE, no w matrix needed
            } else {
              if (is.null(cs$design_w_cache)) cs$design_w_cache = list()
              cs$design_w_cache[[dl]] = obj
            }
          }
          all_cell_states[[job$cell_idx]] <<- cs
          TRUE
        }
        mark_cache_job_done = function(job) {
          attached = attach_cache_job(job)
          cache_jobs_done <<- cache_jobs_done + 1L
          if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar(
            "caching design/SE data",
            cache_jobs_done / max(1L, n_cache_jobs)
          )
          attached
        }
        ready_cache_jobs = cache_jobs[vapply(cache_jobs, function(job) isTRUE(job$ready), logical(1L))]
        pending_cache_jobs = cache_jobs[!vapply(cache_jobs, function(job) isTRUE(job$ready), logical(1L))]
        if (isTRUE(private$verbose) && length(ready_cache_jobs) > 0L) {
          message(sprintf("Loading %d cached design/SE objects from disk...", length(ready_cache_jobs)))
        }
        invisible(lapply(ready_cache_jobs, attach_cache_job))
        cache_jobs_done = length(ready_cache_jobs)
        if (isTRUE(private$verbose) && n_cache_jobs > 0L) {
          private$.draw_labeled_progress_bar(
            "caching design/SE data",
            cache_jobs_done / max(1L, n_cache_jobs)
          )
        }

        if (length(pending_cache_jobs) > 0L) {
          RUN_CACHE_JOB_DETACHED = private$.run_simulation_cache_job
          environment(RUN_CACHE_JOB_DETACHED) = asNamespace("EDI")
          cache_seed = private$seed
          cache_job_seed = function(job_idx) {
            if (is.null(cache_seed)) NULL else as.integer(cache_seed + 1000003L + job_idx)
          }
          # Build a minimal dispatch environment so run_cache_job_local and its
          # callers don't capture the full run() call frame when serialized by clusterApply.
          cache_dispatch_env = new.env(parent = asNamespace("EDI"))
          cache_dispatch_env$RUN_FN   = RUN_CACHE_JOB_DETACHED
          cache_dispatch_env$SEED_VAL = private$seed
          run_cache_job_local = function(job, job_idx) {
            seed = if (is.null(SEED_VAL)) NULL else as.integer(SEED_VAL + 1000003L + job_idx)
            RUN_FN(job = job, seed = seed)
          }
          environment(run_cache_job_local) = cache_dispatch_env
          # Worker wrapper with baseenv() closure — avoids serializing the run() frame.
          fork_cache_worker_fn = function(fn) tryCatch(fn(), error = function(e) e)
          environment(fork_cache_worker_fn) = baseenv()
          if (use_parallel_workers && use_mirai_backend) {
            chunk_size = num_cores_to_use
            cache_chunks = split(seq_along(pending_cache_jobs), (seq_along(pending_cache_jobs) - 1L) %/% chunk_size)
            for (chunk_idxs in cache_chunks) {
              jobs = lapply(chunk_idxs, function(jj) {
                job_j = pending_cache_jobs[[jj]]
                mirai::mirai({
                  RUN_CACHE_JOB_DETACHED(
                    job = job_j,
                    seed = seed_j
                  )
                },
                RUN_CACHE_JOB_DETACHED = RUN_CACHE_JOB_DETACHED,
                job_j = job_j,
                seed_j = cache_job_seed(jj))
              })
              names(jobs) = as.character(chunk_idxs)
              completed = 0L
              while (completed < length(chunk_idxs)) {
                ready = names(jobs)[!vapply(jobs, mirai::unresolved, logical(1L))]
                if (length(ready) == 0L) { Sys.sleep(0.1); next }
                for (nm in ready) {
                  job_res = tryCatch(jobs[[nm]][], error = function(e) e)
                  job_idx = as.integer(nm)
                  if (!inherits(job_res, "error") &&
                      !inherits(job_res, "miraiError") &&
                      !inherits(job_res, "errorValue")) {
                    mark_cache_job_done(pending_cache_jobs[[job_idx]])
                  } else {
                    cache_jobs_done = cache_jobs_done + 1L
                    if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar(
                      "caching design/SE data",
                      cache_jobs_done / max(1L, n_cache_jobs)
                    )
                  }
                  completed = completed + 1L
                }
                jobs[ready] = NULL
              }
            }
          } else if (use_parallel_workers) {
            # fork cluster path (Unix persistent workers)
            chunk_size = num_cores_to_use
            cache_chunks = split(seq_along(pending_cache_jobs), (seq_along(pending_cache_jobs) - 1L) %/% chunk_size)
            for (chunk_idxs in cache_chunks) {
              chunk_fns = lapply(chunk_idxs, function(jj) {
                # Clean closure env: only job + index + the detached fn.
                # Prevents clusterApply from serializing the entire run() call frame.
                e = new.env(parent = asNamespace("EDI"))
                e$job_j    = pending_cache_jobs[[jj]]
                e$jj_j     = jj
                e$local_fn = run_cache_job_local
                f = function() local_fn(job_j, jj_j)
                environment(f) = e
                f
              })
              results = parallel::clusterApply(cl_cache, chunk_fns, fork_cache_worker_fn)
              for (k in seq_along(chunk_idxs)) {
                jj = chunk_idxs[[k]]
                job_res = results[[k]]
                if (!inherits(job_res, "error") && !is.null(job_res)) {
                  mark_cache_job_done(pending_cache_jobs[[jj]])
                } else {
                  cache_jobs_done = cache_jobs_done + 1L
                  if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar(
                    "caching design/SE data",
                    cache_jobs_done / max(1L, n_cache_jobs)
                  )
                }
              }
            }
          } else {
            for (jj in seq_along(pending_cache_jobs)) {
              tryCatch(run_cache_job_local(pending_cache_jobs[[jj]], jj), error = function(e) NULL)
              mark_cache_job_done(pending_cache_jobs[[jj]])
            }
          }
        }
        if (isTRUE(private$verbose)) {
          private$.draw_labeled_progress_bar("caching design/SE data", 1)
          cat("\n", file = stderr())
        }
      }

      # ── Create fork cluster for rep execution (AFTER cell states + caches are populated) ──
      # Workers inherit all_cell_states via copy-on-write at fork time, so per-rep dispatch
      # items are tiny (just ci + rep_i + rep_seed) — no large matrix serialization.
      if (use_fork_cluster) {
        if (!is.null(cl_cache)) {
          try(parallel::stopCluster(cl_cache), silent = TRUE)
          cl_cache = NULL
        }
        RUN_REP_DETACHED_G = private$.run_single_replication_in_worker
        environment(RUN_REP_DETACHED_G) = asNamespace("EDI")
        assign(".edi_sim_cell_states",  all_cell_states,    envir = .GlobalEnv)
        assign(".edi_sim_run_fn",       RUN_REP_DETACHED_G, envir = .GlobalEnv)
        assign(".edi_sim_run_required", all_run_required_v, envir = .GlobalEnv)
        on.exit(suppressWarnings(
          rm(list = intersect(
               c(".edi_sim_cell_states", ".edi_sim_run_fn", ".edi_sim_run_required"),
               ls(envir = .GlobalEnv, all.names = TRUE)
             ), envir = .GlobalEnv, inherits = FALSE)
        ), add = TRUE)
        cl_rep = parallel::makeForkCluster(num_cores)
      }

      # ── Set up and draw initial progress bar (after pregen, before main loop) ──
      if (isTRUE(private$verbose)) {
        private$use_progress_bar = TRUE
        first_active_rep = NA_integer_
        for (r_init in seq_len(private$Nrep)) {
          if (any(vapply(cell_reps_to_run, function(v) r_init %in% v, FALSE))) {
            first_active_rep = r_init; break
          }
        }
        if (!is.na(first_active_rep)) {
          private$current_rep_idx  = first_active_rep
          first_active_cell_init   = which(vapply(cell_reps_to_run, function(v) first_active_rep %in% v, FALSE))[[1L]]
          private$current_cell_idx = first_active_cell_init
          private$tasks_per_rep    = cell_tasks_count[[first_active_cell_init]]
          private$current_task_in_rep_idx = 0L
        } else {
          private$current_rep_idx  = private$Nrep
          private$current_cell_idx = n_cells
          private$tasks_per_rep    = if (length(cell_tasks_count) > 0L) cell_tasks_count[[length(cell_tasks_count)]] else 0L
          private$current_task_in_rep_idx = private$tasks_per_rep
        }
        on.exit({
           if (!is.null(private$progress_bar_drawn)) {
              cat("\n", file = stderr())
              private$progress_bar_drawn = NULL
           }
        }, add = TRUE)
        private$.draw_progress()
      }

      # ── Main simulation loop: outer = rep, inner = DGP cell ──────────────────
      private$rep_elapsed_times        = numeric(private$Nrep)
      private$rep_elapsed_idx          = 0L
      private$session_work_units_done  = 0L
      private$session_start_time       = as.numeric(Sys.time())
      private$rep_start_capture        = private$session_start_time
      # inter_rep_start tracks the wall-clock moment after the previous active rep's
      # post-processing (flush, draw_progress) so that dead space is included in the
      # NEXT rep's elapsed measurement and the ETA is not biased downward.
      inter_rep_start = private$session_start_time
      use_parallel_workers = !force_serial && num_cores > 1L
      num_cores_to_use = if (use_parallel_workers) num_cores else 1L

      # Helper closures for the parallel paths — declared once outside the loops.
      # Fork cluster worker wrapper: baseenv() closure stops serialization chain.
      # Workers look up cell state + run function from globals inherited at fork time
      # (assigned to .GlobalEnv before makeForkCluster — copy-on-write, no serialization).
      # Fork workers run all DGP cells for one replication serially, returning a
      # named list {ci → worker_out}.  Each cell is wrapped in its own tryCatch so
      # a single failing cell does not kill the other cells in the same rep.
      fork_rep_worker_fn = function(item) {
        tryCatch({
          if (!is.null(item$rep_seed)) set.seed(item$rep_seed)
          states  = get(".edi_sim_cell_states",  envir = globalenv())
          run_req = get(".edi_sim_run_required", envir = globalenv())
          run_fn  = get(".edi_sim_run_fn",       envir = globalenv())
          active_ci = which(vapply(run_req, `[[`, FALSE, item$rep_i))
          if (length(active_ci) == 0L) return(list())
          out = vector("list", length(active_ci))
          names(out) = as.character(active_ci)
          for (k in seq_along(active_ci)) {
            ci = active_ci[[k]]
            out[[k]] = tryCatch(
              run_fn(item$rep_i, states[[ci]], progress_cb = NULL, is_forked = TRUE),
              error = function(e) e
            )
          }
          out
        }, error = function(e) e)
      }
      environment(fork_rep_worker_fn) = baseenv()
      is_mirai_failed = function(x) {
        is.null(x) || inherits(x, "error") || inherits(x, "miraiError") || inherits(x, "errorValue")
      }
      is_fork_failed = function(x) is.null(x) || inherits(x, "error")
      mirai_err_msg = function(x) {
        if (is.null(x)) return("Worker returned NULL output.")
        if (inherits(x, "miraiError") || inherits(x, "errorValue")) {
          msg = attr(x, "message")
          if (!is.null(msg)) return(as.character(msg)[1L])
          return(as.character(x)[1L])
        }
        conditionMessage(x)
      }
      seed_val = private$seed

      # ── Fork path: rep-level parallelism ─────────────────────────────────────
      # num_cores reps are dispatched simultaneously; each worker runs all DGP
      # cells for its rep serially.  Workers inherit all_cell_states and
      # all_run_required_v via copy-on-write at fork time — no serialization.
      if (use_fork_cluster) {
        private$current_task_label = "Rep"
        rep_chunks = split(seq_len(private$Nrep),
                           (seq_len(private$Nrep) - 1L) %/% num_cores_to_use)
        for (rep_chunk in rep_chunks) {
          chunk_items = lapply(rep_chunk, function(rep_i) {
            rep_seed = if (!is.null(seed_val)) seed_val + rep_i else NULL
            list(rep_i = rep_i, rep_seed = rep_seed)
          })
          chunk_iter_start   = as.numeric(Sys.time())
          chunk_results_list = parallel::clusterApply(cl_rep, chunk_items, fork_rep_worker_fn)

          # Preliminary per-rep elapsed from parallel-only wall time.
          # Will be corrected after the inner loop to include post-processing overhead.
          n_active_in_chunk = sum(vapply(rep_chunk, function(r)
            any(vapply(all_run_required_v, `[[`, FALSE, r)), logical(1L)))
          chunk_parallel_elapsed = as.numeric(Sys.time()) - chunk_iter_start
          elapsed_per_rep = if (n_active_in_chunk > 0L)
            chunk_parallel_elapsed / n_active_in_chunk else 0
          # Reset rep_start_capture so elapsed_in_rep ticks up correctly during the inner loop.
          private$rep_start_capture = as.numeric(Sys.time())

          for (k in seq_along(rep_chunk)) {
            rep             = rep_chunk[[k]]
            worker_out_list = chunk_results_list[[k]]
            active_cells_rep = which(vapply(all_run_required_v, `[[`, FALSE, rep))

            if (length(active_cells_rep) == 0L) {
              private$current_rep_idx = rep
              private$current_cell_idx = n_cells
              private$current_task_in_rep_idx = if (n_cells > 0L) cell_tasks_count[[n_cells]] else 0L
              private$tasks_per_rep = private$current_task_in_rep_idx
              if (rep %% 100L == 0L || rep == private$Nrep) private$.draw_progress()
              next
            }

            if (is_fork_failed(worker_out_list)) {
              failed_error = private$.make_error_record(
                stage = "worker_execution", rep = rep,
                design = NA_character_, design_params = NULL,
                inference = NA_character_, inference_params = NULL,
                inference_type = NA_character_, inference_type_params = NULL,
                message = if (is.null(worker_out_list)) "Worker returned NULL output."
                          else conditionMessage(worker_out_list),
                metadata = list(rep = rep, backend = "fork_cluster")
              )
              private$.append_errors(list(failed_error))
              if (isTRUE(private$stop_on_error)) private$.abort_from_error_record(failed_error)
              private$progress_count = private$progress_count + length(active_cells_rep)
              private$.draw_progress()
            } else {
              cell_dts   = list()
              total_skip = 0L
              for (ci_nm in names(worker_out_list)) {
                cell_out = worker_out_list[[ci_nm]]
                if (is_fork_failed(cell_out)) {
                  ci = as.integer(ci_nm)
                  cell_err = private$.make_error_record(
                    stage = "worker_execution", rep = rep,
                    design = NA_character_, design_params = NULL,
                    inference = NA_character_, inference_params = NULL,
                    inference_type = NA_character_, inference_type_params = NULL,
                    message = if (is.null(cell_out)) "Worker returned NULL output."
                              else conditionMessage(cell_out),
                    metadata = list(cell_index = ci, rep = rep, backend = "fork_cluster")
                  )
                  private$.append_errors(list(cell_err))
                  if (isTRUE(private$stop_on_error)) private$.abort_from_error_record(cell_err)
                  private$progress_count = private$progress_count + 1L
                } else {
                  private$.append_errors(cell_out$errors)
                  if (!is.null(cell_out$fatal_error))
                    private$.abort_from_error_record(cell_out$fatal_error)
                  if (!is.null(cell_out$results_dt) && nrow(cell_out$results_dt) > 0L)
                    cell_dts[[length(cell_dts) + 1L]] = cell_out$results_dt
                  total_skip = total_skip + cell_out$skipped_count
                }
              }
              if (length(cell_dts) > 0L) {
                all_rep_dt = data.table::rbindlist(cell_dts, use.names = TRUE, fill = TRUE)
                private$.record_batch(all_rep_dt, total_skip)
              } else if (total_skip > 0L) {
                private$.record_batch(NULL, total_skip)
              }
              private$session_work_units_done =
                private$session_work_units_done + length(active_cells_rep)
            }

            private$current_rep_idx = rep
            private$current_cell_idx = n_cells
            private$current_task_in_rep_idx = if (n_cells > 0L) cell_tasks_count[[n_cells]] else 0L
            private$tasks_per_rep = private$current_task_in_rep_idx
            private$.draw_progress()
            if (rep %% private$save_to_disk_every_n_rep == 0L) private$.flush_pending_to_disk()
            private$rep_elapsed_idx = private$rep_elapsed_idx + 1L
            private$rep_elapsed_times[[private$rep_elapsed_idx]] = elapsed_per_rep
          }
          # Back-fill rep_elapsed_times with actual elapsed including post-processing.
          if (n_active_in_chunk > 0L) {
            elapsed_per_rep_actual = (as.numeric(Sys.time()) - chunk_iter_start) / n_active_in_chunk
            idx_end = private$rep_elapsed_idx
            private$rep_elapsed_times[seq(idx_end - n_active_in_chunk + 1L, idx_end)] = elapsed_per_rep_actual
          }
        }
      }

      # ── Mirai / serial path ───────────────────────────────────────────────────
      if (!use_fork_cluster) for (rep in seq_len(private$Nrep)) {
        # Determine which cells need work for this rep.
        active_cells = which(vapply(all_run_required_v, `[[`, FALSE, rep))

        # Fast-forward progress for fully-skipped reps (all cells already done).
        if (length(active_cells) == 0L) {
          private$current_rep_idx = rep
          private$current_cell_idx = n_cells
          private$current_task_in_rep_idx = if (n_cells > 0L) cell_tasks_count[[n_cells]] else 0L
          private$tasks_per_rep = private$current_task_in_rep_idx
          if (rep %% 100L == 0L || rep == private$Nrep) private$.draw_progress()
          next
        }

        private$rep_start_capture = inter_rep_start
        private$current_rep_idx = rep
        private$current_cell_idx = active_cells[[1L]]

        if (use_mirai_backend) {
            private$current_task_label = "Chunk Cells"
            chunk_size = num_cores_to_use
            cell_chunks = split(active_cells, (seq_along(active_cells) - 1L) %/% chunk_size)
            for (chunk_cells in cell_chunks) {
              private$current_cell_idx = chunk_cells[[1L]]
              private$tasks_per_rep    = length(chunk_cells)
              private$current_task_in_rep_idx = 0L
              private$last_progress_draw_time = 0
              private$.draw_progress()

              RUN_REP_DETACHED = private$.run_single_replication_in_worker
              environment(RUN_REP_DETACHED) = asNamespace("EDI")

              jobs = lapply(chunk_cells, function(ci) {
                cs_i = all_cell_states[[ci]]
                RUN_REP_ci = function(rep_i) RUN_REP_DETACHED(rep_i, cs_i, progress_cb = NULL, is_forked = FALSE)
                rep_seed = if (!is.null(seed_val)) seed_val + rep else NULL
                mirai::mirai({
                  if (!is.null(rep_seed)) set.seed(rep_seed)
                  RUN_REP_ci(rep_i)
                }, RUN_REP_ci = RUN_REP_ci, rep_i = rep, rep_seed = rep_seed)
              })
              names(jobs) = as.character(chunk_cells)
              chunk_results = vector("list", length(chunk_cells))
              names(chunk_results) = as.character(chunk_cells)
              completed_in_chunk = 0L
              while (completed_in_chunk < length(chunk_cells)) {
                ready = names(jobs)[!vapply(jobs, mirai::unresolved, logical(1L))]
                if (length(ready) == 0L) { Sys.sleep(0.1); next }
                for (nm in ready) {
                  if (is.null(chunk_results[[nm]])) {
                    job_res = tryCatch(jobs[[nm]][], error = function(e) e)
                    chunk_results[[nm]] = job_res
                    if (!is_mirai_failed(job_res) && !is.null(job_res$fatal_error)) {
                      pending_jobs = jobs[setdiff(names(jobs), ready)]
                      if (length(pending_jobs) > 0L)
                        invisible(lapply(pending_jobs, function(job) try(mirai::stop_mirai(job), silent = TRUE)))
                    }
                    completed_in_chunk = completed_in_chunk + 1L
                    private$session_work_units_done = private$session_work_units_done + 1L
                  }
                }
                jobs[ready] = NULL
                private$current_cell_idx = chunk_cells[[min(length(chunk_cells), completed_in_chunk + 1L)]]
                private$current_task_in_rep_idx = completed_in_chunk
                private$.draw_progress()
              }
              chunk_results_ok = chunk_results[!vapply(chunk_results, is_mirai_failed, logical(1L))]
              if (length(chunk_results_ok) > 0L) {
                all_chunk_dt = data.table::rbindlist(lapply(chunk_results_ok, function(w) w$results_dt), use.names = TRUE, fill = TRUE)
                all_chunk_skipped = sum(vapply(chunk_results_ok, function(w) w$skipped_count, 0L))
                private$.append_errors(unlist(lapply(chunk_results_ok, `[[`, "errors"), recursive = FALSE))
                private$current_cell_idx = chunk_cells[[length(chunk_cells)]]
                private$current_task_in_rep_idx = length(chunk_cells)
                private$.record_batch(all_chunk_dt, all_chunk_skipped)
              }
              fatal_errors = unlist(lapply(chunk_results_ok, function(x) {
                if (is.null(x$fatal_error)) list() else list(x$fatal_error)
              }), recursive = FALSE)
              if (length(fatal_errors) > 0L) private$.abort_from_error_record(fatal_errors[[1L]])
              failed_indices = which(vapply(chunk_results, is_mirai_failed, logical(1L)))
              if (length(failed_indices) > 0L) {
                failed_errors = lapply(failed_indices, function(idx) {
                  ci = chunk_cells[[idx]]
                  private$.make_error_record(
                    stage = "worker_execution", rep = rep,
                    design = NA_character_, design_params = NULL,
                    inference = NA_character_, inference_params = NULL,
                    inference_type = NA_character_, inference_type_params = NULL,
                    message = mirai_err_msg(chunk_results[[idx]]),
                    metadata = list(cell_index = ci, rep = rep, backend = "mirai")
                  )
                })
                private$.append_errors(failed_errors)
                if (isTRUE(private$stop_on_error)) private$.abort_from_error_record(failed_errors[[1L]])
                private$progress_count = private$progress_count + length(failed_indices)
                private$.draw_progress()
              }
            }
        } else {
          # ── Serial path ───────────────────────────────────────────────────────
          private$current_task_label = "Des/Inf"
          for (cell_idx in active_cells) {
            cell_state = all_cell_states[[cell_idx]]
            private$current_cell_idx  = cell_idx
            private$current_n         = cell_state$n
            private$current_p         = cell_state$p
            private$current_betaT     = cell_state$betaT
            private$current_cond_exp_func_model = cell_state$cond_exp_func_model
            private$current_response_type       = cell_state$response_type
            private$tasks_per_rep     = cell_tasks_count[[cell_idx]]
            private$current_task_in_rep_idx = 0L
            private$last_progress_draw_time = 0
            private$.draw_progress()
            worker_out = private$.run_single_replication_in_worker(rep, cell_state, progress_cb = private$.advance_progress, is_forked = FALSE)
            if (private$keep_all_intermediate_data) {
              rep_data = if (isTRUE(private$random_X_draws)) private$.generate_data() else private$.generate_data_from_X(shared_X_draws[[paste(private$current_n, private$current_p, sep = "|")]])
              X = rep_data$X; y_linear_model = rep_data$y_linear_model
              true_mean_diff_ate = private$compute_true_mean_diff_ate(y_linear_model)
              rep_slot = (cell_idx - 1L) * private$Nrep + rep
              for (di in seq_along(private$design_classes)) {
                design_gen  = private$design_classes[[di]]; design_name = private$design_labels[[di]]
                design_extra = if (!is.null(private$design_params)) private$design_params[[di]] else list()
                des_obj = private$.build_design(design_gen, X, y_linear_model, design_extra, rep_data = rep_data)
                private$all_intermediate_data[[rep_slot]]$designs[[design_name]] = des_obj
                if (is.null(des_obj)) next
                for (ii in seq_along(private$inference_classes)) {
                  inf_gen = private$inference_classes[[ii]]; inf_name = private$inference_labels[[ii]]
                  inf_ctor_extra = private$inference_constructor_params[[ii]]
                  inf_obj = do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra))
                  private$all_intermediate_data[[rep_slot]]$inferences[[design_name]][[inf_name]] = inf_obj
                }
              }
              private$all_intermediate_data[[rep_slot]]$y_linear_model = y_linear_model
            }
            private$session_work_units_done = private$session_work_units_done + 1L
            if (!is(worker_out, "try-error") && !is.null(worker_out)) {
              private$.append_errors(worker_out$errors)
              if (!is.null(worker_out$fatal_error)) private$.abort_from_error_record(worker_out$fatal_error)
              n_res = if (is.null(worker_out$results_dt)) 0L else nrow(worker_out$results_dt)
              private$current_task_in_rep_idx = worker_out$skipped_count + n_res
              private$.record_batch(worker_out$results_dt, worker_out$skipped_count)
            } else {
              worker_error = private$.make_error_record(
                stage = "worker_execution", rep = rep,
                design = NA_character_, design_params = NULL,
                inference = NA_character_, inference_params = NULL,
                inference_type = NA_character_, inference_type_params = NULL,
                message = if (is.null(worker_out)) "Worker returned NULL output." else as.character(worker_out),
                metadata = list(cell_index = cell_idx, num_cores = num_cores_to_use)
              )
              private$.append_errors(list(worker_error))
              if (isTRUE(private$stop_on_error)) private$.abort_from_error_record(worker_error)
              private$current_task_in_rep_idx = private$tasks_per_rep
              private$progress_count = private$progress_count + cell_tasks_count[[cell_idx]]
              private$.draw_progress()
            }
          }
        }

        # Record elapsed time for this rep and update ETA.
        # flush and draw_progress happen AFTER recording so they're included in the
        # NEXT rep's inter_rep_start measurement rather than being dead space.
        private$current_rep_idx = rep
        private$current_cell_idx = n_cells
        private$current_task_in_rep_idx = if (n_cells > 0L) cell_tasks_count[[n_cells]] else 0L
        private$tasks_per_rep = private$current_task_in_rep_idx
        private$.draw_progress()
        if (rep %% private$save_to_disk_every_n_rep == 0L) private$.flush_pending_to_disk()
        inter_rep_now = as.numeric(Sys.time())
        private$rep_elapsed_idx = private$rep_elapsed_idx + 1L
        private$rep_elapsed_times[[private$rep_elapsed_idx]] = inter_rep_now - inter_rep_start
        inter_rep_start = inter_rep_now
      }
      private$.finish_run()
      invisible(self)
    },
    #' @description Retrieve the stored intermediate data (design and inference objects)
    #' for every replication. Only available if \code{keep_all_intermediate_data = TRUE}
    #' was passed to the constructor.
    #'
    #' @return A nested list containing the intermediate data for each replication,
    #'   or \code{NULL} if not recorded.
    get_all_intermediate_data = function() {
      if (!private$has_run) stop("Call $run() first.")
      if (!private$keep_all_intermediate_data) return(NULL)
      private$all_intermediate_data
    },
    # ── clear_all_intermediate_data_and_gc() ─────────────────────────────────
    #' @description Release all stored intermediate data and invoke the garbage collector.
    #' Useful after inspecting intermediate results to free memory before
    #' further processing.  Sets the internal store to \code{NULL} and calls
    #' \code{gc()}.
    #'
    #' @return The \code{SimulationFramework} object itself (invisibly).
    clear_all_intermediate_data_and_gc = function() {
      private$all_intermediate_data = NULL
      gc()
      invisible(self)
    }
  ),
  # ── private ────────────────────────────────────────────────────────────────
  private = list(
    .finish_run = function() {
      private$has_run = TRUE
      private$.flush_pending_to_disk()
      if (private$.results_file_format(private$results_filename) == "csv.bz2" &&
          file.exists(private$.results_staging_filename()) &&
          private$results_idx > 0L) {
        private$.sync_results_bz2_from_staging()
      }
      private$.cleanup_results_staging_file()
      total_elapsed = private$prior_elapsed_secs + (as.numeric(Sys.time()) - private$simulation_start_time)
      tryCatch(saveRDS(total_elapsed, private$.elapsed_file()), error = function(e) NULL)
      clear_result_key_store_cpp() # Free memory
      if (isTRUE(private$use_progress_bar)) {
        private$current_cell_idx = private$total_cells
        private$current_rep_idx = private$Nrep
        private$current_task_in_rep_idx = private$tasks_per_rep
        private$.draw_simulation_progress_bars()
        cat("\n", file = stderr())
        private$progress_bar_drawn = NULL
      }
    },
    response_type_values = NULL,
    current_response_type = NULL,
    n_values         = NULL,
    p_values         = NULL,
    cond_exp_func_model_values = NULL,
    Nrep             = NULL,
    betaT_values     = NULL,
    alpha            = NULL,
    B_boot           = NULL,
    r_rand           = NULL,
    pval_epsilon     = NULL,
    sd_noise             = NULL,
	    n_ordinal_levels     = NULL,
	    proportion_epsilon   = NULL,
	    phi_proportion       = NULL,
	    k_survival           = NULL,
	    incidence_clamp      = NULL,
	    proportion_clamp     = NULL,
	    count_clamp          = NULL,
	    survival_clamp       = NULL,
	    survival_min_time    = NULL,
    count_min_rate       = NULL,
    count_shift          = NULL,
    norm_sq_beta_vec     = NULL,
    X_mat                = NULL,
    num_cores            = NULL,
    seed                 = NULL,
    cov_draw_method      = NULL,
    cov_draw_method_args = NULL,
    random_X_draws       = NULL,
    prob_censoring       = NULL,
    custom_replication_data_generator = NULL,
    custom_apply_treatment_and_noise = NULL,
    make_estimand_fn = NULL,
    dgp_params = NULL,
    custom_dgp = NULL,
    verbose          = NULL,
    results_filename = NULL,
    simulation_start_time = NULL,
    initial_progress_count = NULL,
    continue_from_last_result_row = NULL,
    reuse_cache = TRUE,
    stop_on_error = TRUE,
    save_to_disk_every_n_rep = 50L,
    save_model_control_fits = TRUE,
    pending_file_rows = NULL,
    design_params    = NULL,
    inference_constructor_params = NULL,
    inference_type_params = NULL,
    inf_types        = NULL,
    design_classes   = NULL,
    inference_classes = NULL,
    design_labels    = NULL,
    inference_labels = NULL,
    turn_off_asserts_for_speed = NULL,
    raw_results               = NULL,
    error_log                 = NULL,
    results_idx               = 0L,
    all_intermediate_data     = NULL,    keep_all_intermediate_data = FALSE,
    has_run                   = FALSE,
    exact_warned_classes      = NULL,
    valid_combos              = NULL,
    seen_combo_keys  = NULL,
    seen_result_keys = NULL,
    total_cells              = 0L,
    current_cell_idx         = 0L,
    last_progress_draw_time  = 0,
    current_task_label       = "Des/Inf",
    current_rep_idx          = 0L,
    current_task_in_rep_idx  = 0L,
    tasks_per_rep            = 0L,
    progress_total           = 0L,
    progress_count           = 0L,
    rep_elapsed_times           = numeric(0),
    rep_elapsed_idx             = 0L,
    prior_elapsed_secs          = 0,
    session_start_time          = NULL,
    rep_start_capture           = NULL,
    n_active_reps_total         = 0L,
    n_total_active_work_units   = 0L,
    session_work_units_done     = 0L,
    progress_bar             = NULL,
    use_progress_bar         = FALSE,
    progress_log_interval    = 0L,
    progress_bar_drawn       = NULL,
    n_progress_lines         = 5L,
    .ensure_mirai_daemons = function(n) {
      s = tryCatch(mirai::status(), error = function(e) list(connections = 0L))
      n_running = if (is.numeric(s$connections) && length(s$connections) == 1L) as.integer(s$connections) else 0L
      if (n_running != as.integer(n)) mirai::daemons(as.integer(n))
      # everywhere() in mirai 2.x is asynchronous; collect before submitting real tasks.
      setup_tasks = tryCatch(
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
          invisible(NULL)
        }),
        error = function(e) NULL
      )
      if (!is.null(setup_tasks)) {
        tryCatch(mirai::collect_mirai(setup_tasks), error = function(e) invisible(NULL))
      }
      invisible(NULL)
    },
    # ── Design spec parsing ───────────────────────────────────────────────────
    .parse_design_classes_and_params = function(spec, eval_env) {
      if (is.null(spec)) {
        classes = private$.default_design_classes()
        return(list(classes = classes, params = lapply(classes, function(...) list())))
      }
      if (!is.list(spec))
        stop("design_classes_and_params must be NULL or a list")
      nm = names(spec)
      if (is.null(nm)) nm = rep("", length(spec))
      classes = vector("list", length(spec))
      params  = vector("list", length(spec))
      for (i in seq_along(spec)) {
        entry_name = nm[[i]]
        entry      = spec[[i]]
        if (inherits(entry, "R6ClassGenerator")) {
          classes[[i]] = entry
          params[[i]]  = list()
          next
        }
        if (!nzchar(entry_name)) {
          stop(
            "design_classes_and_params[[", i, "]] must be an R6 class generator ",
            "or a named parameter list whose name is the design class"
          )
        }
        cls = private$.resolve_design_class(entry_name, eval_env)
        if (!inherits(cls, "R6ClassGenerator"))
          stop("design class '", entry_name, "' is not an R6 class generator")
        if (is.null(entry)) {
          entry = list()
        } else if (!is.list(entry)) {
          stop("design_classes_and_params[['", entry_name, "']] must be a list of parameters")
        }
        classes[[i]] = cls
        params[[i]]  = entry
      }
      list(classes = classes, params = params)
    },
    .resolve_design_class = function(class_name, eval_env) {
      if (exists(class_name, envir = eval_env, inherits = TRUE))
        return(get(class_name, envir = eval_env, inherits = TRUE))
      ns = asNamespace("EDI")
      if (exists(class_name, envir = ns, inherits = FALSE))
        return(get(class_name, envir = ns, inherits = FALSE))
      stop("could not find design class '", class_name, "'")
    },
    # ── Inference spec parsing ────────────────────────────────────────────────
    .parse_inference_classes_and_params = function(spec, eval_env) {
      if (is.null(spec)) {
        classes = private$.default_inference_classes()
        return(list(classes = classes, params = lapply(classes, function(...) list())))
      }
      if (!is.list(spec))
        stop("inference_classes_and_params must be NULL or a list")
      nm = names(spec)
      if (is.null(nm)) nm = rep("", length(spec))
      classes = vector("list", length(spec))
      params  = vector("list", length(spec))
      for (i in seq_along(spec)) {
        entry_name = nm[[i]]
        entry      = spec[[i]]
        if (inherits(entry, "R6ClassGenerator")) {
          classes[[i]] = entry
          params[[i]]  = list()
          next
        }
        if (!nzchar(entry_name)) {
          stop(
            "inference_classes_and_params[[", i, "]] must be an R6 class generator ",
            "or a named parameter list whose name is the inference class"
          )
        }
        cls = private$.resolve_inference_class(entry_name, eval_env)
        if (!inherits(cls, "R6ClassGenerator"))
          stop("inference class '", entry_name, "' is not an R6 class generator")
        if (is.null(entry)) {
          entry = list()
        } else if (!is.list(entry)) {
          stop("inference_classes_and_params[['", entry_name, "']] must be a list of parameters")
        }
        private$.validate_r6_init_args(cls, entry, "inference_classes_and_params")
        classes[[i]] = cls
        params[[i]]  = entry
      }
      list(classes = classes, params = params)
    },
    .resolve_inference_class = function(class_name, eval_env) {
      if (exists(class_name, envir = eval_env, inherits = TRUE))
        return(get(class_name, envir = eval_env, inherits = TRUE))
      ns = asNamespace("EDI")
      if (exists(class_name, envir = ns, inherits = FALSE))
        return(get(class_name, envir = ns, inherits = FALSE))
      stop("could not find inference class '", class_name, "'")
    },
    .validate_r6_init_args = function(r6gen, args, arg_name) {
      if (length(args) == 0L) return(invisible(TRUE))
      if (is.null(names(args)) || any(!nzchar(names(args)))) {
        stop(arg_name, " for ", r6gen$classname,
             " must be a named list of constructor arguments")
      }
      init_fn = get_r6_init_fn(r6gen)
      if (is.null(init_fn)) {
        stop(r6gen$classname, " has no discoverable initialize() constructor; ",
             arg_name, " parameters cannot be validated")
      }
      fn_formals = names(formals(init_fn))
      if ("..." %in% fn_formals) return(invisible(TRUE))
      bad = setdiff(names(args), fn_formals)
      if (length(bad)) {
        stop(
          arg_name, " for ", r6gen$classname,
          " contains constructor argument(s) not accepted by initialize(): ",
          paste(bad, collapse = ", ")
        )
      }
      invisible(TRUE)
    },
    # ── Inference type parsing and invocation args ────────────────────────────
    .parse_inference_types_and_params = function(spec, valid_inf_types) {
      if (is.null(spec)) {
        return(stats::setNames(lapply(valid_inf_types, function(...) list()),
                               valid_inf_types))
      }
      if (!is.list(spec) || is.null(names(spec)) || any(!nzchar(names(spec)))) {
        stop("inference_types_and_params must be NULL or a named list")
      }
      bad = setdiff(names(spec), valid_inf_types)
      if (length(bad)) {
        stop("Invalid inference_types_and_params names: ", paste(bad, collapse = ", "),
             ".  Valid values: ", paste(valid_inf_types, collapse = ", "))
      }
      spec = spec[!duplicated(names(spec))]
      for (inf_type in names(spec)) {
        if (is.null(spec[[inf_type]])) {
          spec[[inf_type]] = list()
        } else if (!is.list(spec[[inf_type]])) {
          stop("inference_types_and_params[['", inf_type, "']] must be a named list")
        } else if (length(spec[[inf_type]]) > 0L &&
                   (is.null(names(spec[[inf_type]])) || any(!nzchar(names(spec[[inf_type]]))))) {
          stop("inference_types_and_params[['", inf_type, "']] must be a named list")
        }
      }
      spec
    },
    .has_inf_type = function(inf_type) {
      inf_type %in% private$inf_types
    },
    .any_inf_type = function(inf_types) {
      any(inf_types %in% private$inf_types)
    },
    .inf_type_method_name = function(inf_type) {
      switch(inf_type,
        asymp_ci   = "compute_asymp_confidence_interval",
        asymp_pval = "compute_asymp_two_sided_pval",
        exact_ci   = "compute_exact_confidence_interval",
        exact_pval = "compute_exact_two_sided_pval_for_treatment_effect",
        boot_ci    = "compute_bootstrap_confidence_interval",
        boot_pval  = "compute_bootstrap_two_sided_pval",
        rand_ci    = "compute_rand_confidence_interval",
        rand_pval  = "compute_rand_two_sided_pval",
        stop("Unknown inference type: ", inf_type)
      )
    },
    .args_for_inf_type = function(inf_obj, inf_type, defaults = list()) {
      user_args = private$inference_type_params[[inf_type]]
      if (is.null(user_args)) user_args = list()
      method_name = private$.inf_type_method_name(inf_type)
      private$.validate_method_args(inf_obj, method_name, user_args, inf_type)
      modifyList(defaults, user_args)
    },
    .validate_method_args = function(inf_obj, method_name, args, inf_type) {
      if (length(args) == 0L) return(invisible(TRUE))
      fn = tryCatch(inf_obj[[method_name]], error = function(e) NULL)
      if (!is.function(fn)) {
        stop("Cannot validate parameters for ", inf_type, ": function ",
             method_name, "() is not available on ", class(inf_obj)[1L])
      }
      fn_formals = names(formals(fn))
      if ("..." %in% fn_formals) return(invisible(TRUE))
      bad = setdiff(names(args), fn_formals)
      if (length(bad)) {
        stop("inference_types_and_params[['", inf_type, "']] contains argument(s) ",
             "not accepted by ", method_name, "(): ", paste(bad, collapse = ", "))
      }
      invisible(TRUE)
    },
    .params_for_inference_type_to_str = function(inference_type) {
      ps = private$.params_to_str(private$inference_type_params[[inference_type]])
      if (nchar(ps) > 0L) paste0(inference_type, "(", ps, ")") else ""
    },
    # ── R6 formals helpers ────────────────────────────────────────────────────
    # Filter an arg list to only those accepted by r6gen's initialize().
    # If initialize accepts '...' or cannot be found, all args pass through.
    .filter_by_formals = function(r6gen, args) {
      if (length(args) == 0L) return(args)
      init_fn = get_r6_init_fn(r6gen)
      if (is.null(init_fn)) return(args)
      fn_formals = names(formals(init_fn))
      if ("..." %in% fn_formals) return(args)
      args[names(args) %in% fn_formals]
    },
    .has_private_method_on_object = function(obj, method_name) {
      exists(method_name, envir = obj$.__enclos_env__$private, inherits = FALSE)
    },
    # Serialize a named params list to "k=v, ..." string.
    .params_to_str = function(p) {
      if (is.null(p) || length(p) == 0L) return("")
      kv = mapply(function(k, v) paste0(k, "=", paste(deparse(v), collapse = "")),
                  names(p), p, SIMPLIFY = TRUE)
      paste(kv, collapse = ", ")
    },
    .build_param_grid = function(n_values, p_values, betaT_values, cond_exp_func_model_values, response_type_values) {
      grid = data.table::as.data.table(expand.grid(
        response_type = response_type_values,
        cond_exp_func_model = cond_exp_func_model_values,
        n = n_values,
        p = p_values,
        betaT = betaT_values,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      ))
      grid = grid[!(cond_exp_func_model == "nonlinear" & p < 5L)]
      if (nrow(grid) == 0L)
        stop("No valid simulation cells remain after filtering cond_exp_func_model / p combinations")
      grid
    },
    .format_values = function(x) {
      if (length(x) == 1L) as.character(x) else paste0("c(", paste(x, collapse = ", "), ")")
    },
    .load_existing_results = function() {
      empty_dt = data.table::data.table(
        response_type = character(),
        cond_exp_func_model = character(),
        n = integer(),
        p = integer(),
        betaT = numeric(),
        rep = integer(),
        design = character(),
        inference = character(),
        inference_type = character(),
        estimate = numeric(),
        ci_lo = numeric(),
        ci_hi = numeric(),
        pval = numeric(),
        true_estimand = numeric()
      )
      if (!isTRUE(private$continue_from_last_result_row) || !file.exists(private$results_filename))
        return(private$.load_existing_results_from_staging_or_empty(empty_dt))
      if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar("loading previous results", 0)
      dt = private$.read_results_file(show_progress = FALSE)
      if (isTRUE(private$verbose)) {
        private$.draw_labeled_progress_bar("loading previous results", 1)
        cat("\n", file = stderr())
      }
      if (!"response_type" %in% names(dt))
        dt[, response_type := private$response_type_values[[1L]]]
      dt = dt[response_type %in% private$response_type_values]
      for (nm in names(empty_dt)) {
        if (!nm %in% names(dt))
          dt[, (nm) := empty_dt[[nm]]]
      }
      dt[, response_type := as.character(response_type)]
      dt[, cond_exp_func_model := as.character(cond_exp_func_model)]
      dt[, design := as.character(design)]
      dt[, inference := as.character(inference)]
      dt[, inference_type := as.character(inference_type)]
      dt[, n := as.integer(n)]
      dt[, p := as.integer(p)]
      dt[, rep := as.integer(rep)]
      dt[, betaT := as.numeric(betaT)]
      dt[, estimate := as.numeric(estimate)]
      dt[, ci_lo := as.numeric(ci_lo)]
      dt[, ci_hi := as.numeric(ci_hi)]
      dt[, pval := as.numeric(pval)]
      dt[, true_estimand := as.numeric(true_estimand)]
      dt[, names(empty_dt), with = FALSE]
    },
    .load_existing_results_from_staging_or_empty = function(empty_dt) {
      if (private$.results_file_format(private$results_filename) != "csv.bz2")
        return(empty_dt)
      staging_filename = private$.results_staging_filename()
      if (!file.exists(staging_filename))
        return(empty_dt)
      data.table::fread(staging_filename, showProgress = FALSE)
    },
    .results_file_format = function(filename) {
      if (grepl("\\.csv\\.bz2$", filename, ignore.case = TRUE)) return("csv.bz2")
      if (grepl("\\.csv$", filename, ignore.case = TRUE)) return("csv")
      NA_character_
    },
    .results_staging_filename = function() {
      path = private$.results_output_path()
      file.path(
        dirname(path),
        paste0(
          sub("\\.csv\\.bz2$", "", basename(path), ignore.case = TRUE),
          "__staging.csv"
        )
      )
    },
    .results_output_path = function() {
      # Absolute path forms: Unix ('/...'), Windows drive letter ('C:/...',
      # 'C:\...'), or Windows UNC ('\\server\share\...').
      is_absolute = grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", private$results_filename)
      if (is_absolute) {
        private$results_filename
      } else {
        file.path(getwd(), private$results_filename)
      }
    },
    .simulation_cache_dir = function() {
      path = private$.results_output_path()
      stem = sub("\\.csv(\\.bz2)?$", "", basename(path), ignore.case = TRUE)
      file.path(dirname(path), paste0(stem, "__cache"))
    },
    .elapsed_file = function() {
      path = private$.results_output_path()
      stem = sub("\\.csv(\\.bz2)?$", "", basename(path), ignore.case = TRUE)
      file.path(dirname(path), paste0(stem, "__elapsed.rds"))
    },
    .safe_cache_component = function(x) {
      x = gsub("[^A-Za-z0-9_.=-]+", "_", as.character(x))
      x = gsub("^_+|_+$", "", x)
      if (!nzchar(x)) "cache" else x
    },
    .hash_object = function(x) {
      digest::digest(x, algo = "md5")
    },
    .simulation_cache_file = function(cs, design_label, design_params, cache_type) {
      cache_key = private$.hash_object(list(
        cache_format_version = 1L,
        cache_type = cache_type,
        response_type = cs$response_type,
        cond_exp_func_model = cs$cond_exp_func_model,
        n = as.integer(cs$n),
        p = as.integer(cs$p),
        betaT = as.numeric(cs$betaT),
        Nrep = as.integer(private$Nrep),
        seed = private$seed,
        random_X_draws = private$random_X_draws,
        shared_X = cs$shared_X,
        X_mat = cs$X_mat,
        cov_draw_method_args = cs$cov_draw_method_args,
        norm_sq_beta_vec = cs$norm_sq_beta_vec,
        design_label = design_label,
        design_params = design_params
      ))
      filename = paste(
        private$.safe_cache_component(cache_type),
        private$.safe_cache_component(cs$response_type),
        private$.safe_cache_component(cs$cond_exp_func_model),
        paste0("n", cs$n),
        paste0("p", cs$p),
        private$.safe_cache_component(design_label),
        cache_key,
        sep = "__"
      )
      file.path(private$.simulation_cache_dir(), paste0(filename, ".rds"))
    },
    .load_simulation_cache_object = function(cs, design_label, design_params,
                                             cache_type, reps_needing,
                                             restore_rng = FALSE,
                                             cache_file = NULL) {
      if (is.null(cache_file)) {
        cache_file = private$.simulation_cache_file(cs, design_label, design_params, cache_type)
      }
      if (!file.exists(cache_file)) return(NULL)
      cache_record = tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.null(cache_record)) return(NULL)
      obj = if (is.list(cache_record) &&
                identical(cache_record$cache_format_version, 1L) &&
                !is.null(cache_record$value)) {
        cache_record$value
      } else {
        cache_record
      }
      if (identical(cache_type, "design_w")) {
        if (is.list(obj) && isTRUE(obj$blocking_design)) return(obj)  # sentinel
        if (!is.list(obj) || is.null(obj$ws) || is.null(obj$rep_to_col)) return(NULL)
        rep_names = names(obj$rep_to_col)
        if (is.null(rep_names) ||
            any(!as.character(reps_needing) %in% rep_names)) return(NULL)
        if (!is.matrix(obj$ws) || nrow(obj$ws) != as.integer(cs$n)) return(NULL)
      } else {
        return(NULL)
      }
      if (isTRUE(restore_rng) &&
          is.list(cache_record) &&
          identical(cache_record$cache_format_version, 1L) &&
          !is.null(cache_record$rng_after)) {
        assign(".Random.seed", cache_record$rng_after, envir = .GlobalEnv)
      }
      obj
    },
    .save_simulation_cache_object = function(obj, cs, design_label, design_params, cache_type) {
      cache_dir = private$.simulation_cache_dir()
      if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      cache_file = private$.simulation_cache_file(cs, design_label, design_params, cache_type)
      tmp = tempfile(
        paste0(".", basename(cache_file), "."),
        tmpdir = dirname(cache_file),
        fileext = ".tmp"
      )
      on.exit(unlink(tmp, force = TRUE), add = TRUE)
      cache_record = list(
        cache_format_version = 1L,
        value = obj,
        rng_after = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else NULL
      )
      saveRDS(cache_record, tmp, version = 2)
      if (!file.rename(tmp, cache_file)) {
        file.copy(tmp, cache_file, overwrite = TRUE)
        unlink(tmp, force = TRUE)
      }
      invisible(cache_file)
    },
    .run_simulation_cache_job = function(job, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      ns = asNamespace("EDI")
      prev_nc = ns$edi_env$num_cores_override
      on.exit(assign("num_cores_override", prev_nc, envir = ns$edi_env), add = TRUE)
      # Use the thread budget computed by the caller; default to 1 to stay safe when
      # multiple workers run simultaneously (prevents N*M thread explosion).
      worker_nc = if (!is.null(job$num_cores_for_worker)) as.integer(job$num_cores_for_worker) else 1L
      assign("num_cores_override", worker_nc, envir = ns$edi_env)

      save_cache_record = function(obj, cache_file) {
        cache_dir = dirname(cache_file)
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        tmp = tempfile(
          paste0(".", basename(cache_file), "."),
          tmpdir = cache_dir,
          fileext = ".tmp"
        )
        on.exit(unlink(tmp, force = TRUE), add = TRUE)
        cache_record = list(
          cache_format_version = 1L,
          value = obj,
          rng_after = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
          } else NULL
        )
        saveRDS(cache_record, tmp, version = 2)
        if (!file.rename(tmp, cache_file)) {
          file.copy(tmp, cache_file, overwrite = TRUE)
          unlink(tmp, force = TRUE)
        }
        invisible(cache_file)
      }

      cs = job$cs
      design_gen = asNamespace("EDI")[[job$design_gen_class]]
      d_pg = do.call(
        design_gen$new,
        c(list(response_type = cs$response_type, n = cs$n), job$design_params)
      )
      d_pg$add_all_subjects_to_experiment(as.data.frame(cs$shared_X))

      if (identical(job$cache_type, "design_w")) {
        # Non-pregen blocking designs use closed-form CMH SE — no w-pregeneration needed.
        if (!isTRUE(job$is_pregen) && d_pg$is_blocking_design()) {
          save_cache_record(list(blocking_design = TRUE), job$cache_file)
          return(list(status = "skipped_blocking_design", cache_file = job$cache_file))
        }
        w_mat = d_pg$draw_ws_according_to_design(length(job$reps_needing))
        storage.mode(w_mat) = "integer"
        obj = list(
          ws = w_mat,
          rep_to_col = setNames(seq_along(job$reps_needing), as.character(job$reps_needing))
        )
        save_cache_record(obj, job$cache_file)
        return(list(status = "created", cache_file = job$cache_file))
      }

      list(status = "skipped_unknown_cache_type", cache_file = job$cache_file)
    },
    .reps_needing_design_cache = function(active_reps, des_combos, cs, n_existing) {
      reps_needing = active_reps
      if (n_existing > 0L) {
        n_dc = length(des_combos)
        nb   = length(active_reps)
        c_des  = vapply(des_combos, `[[`, "", "design")
        c_inf  = vapply(des_combos, `[[`, "", "inference")
        c_type = vapply(des_combos, `[[`, "", "inference_type")
        all_ex = check_in_result_key_store_cpp(
          rep(cs$response_type,       nb * n_dc),
          rep(cs$cond_exp_func_model, nb * n_dc),
          rep(as.integer(cs$n),       nb * n_dc),
          rep(as.integer(cs$p),       nb * n_dc),
          rep(as.numeric(cs$betaT),   nb * n_dc),
          rep(as.integer(active_reps), each = n_dc),
          rep(c_des,  nb),
          rep(c_inf,  nb),
          rep(c_type, nb)
        )
        exists_mat = matrix(all_ex, nrow = n_dc)
        reps_needing = active_reps[colSums(exists_mat) < n_dc]
      }
      reps_needing
    },
    .copy_binary_stream = function(from, to, chunk_size = 1024L * 1024L) {
      repeat {
        bytes = readBin(from, what = "raw", n = chunk_size)
        if (length(bytes) == 0L) break
        writeBin(bytes, to)
      }
      invisible(NULL)
    },
    .read_results_file = function(show_progress = FALSE) {
      format = private$.results_file_format(private$results_filename)
      if (identical(format, "csv")) {
        return(data.table::fread(private$results_filename, showProgress = show_progress))
      }
      if (!identical(format, "csv.bz2")) {
        stop("Unsupported results file format: ", private$results_filename)
      }
      staging_filename = private$.results_staging_filename()
      if (file.exists(staging_filename)) {
        return(data.table::fread(staging_filename, showProgress = show_progress))
      }
      extracted_file = tempfile("simulation_framework_results_", fileext = ".csv")
      input_con = bzfile(private$.results_output_path(), open = "rb")
      output_con = file(extracted_file, open = "wb")
      on.exit(try(close(input_con), silent = TRUE), add = TRUE)
      on.exit(try(close(output_con), silent = TRUE), add = TRUE)
      on.exit(unlink(extracted_file, force = TRUE), add = TRUE)
      private$.copy_binary_stream(input_con, output_con)
      close(output_con)
      close(input_con)
      data.table::fread(extracted_file, showProgress = show_progress)
    },
    .result_key_for_values = function(response_type, cond_exp_func_model, n, p, betaT, rep, design, inference, inference_type) {
      paste(
        response_type,
        cond_exp_func_model,
        n,
        p,
        betaT,
        rep,
        design,
        inference,
        inference_type,
        sep = "|"
      )
    },
    .result_key = function(rep, design, inference, inference_type) {
      private$.result_key_for_values(
        private$current_response_type,
        private$current_cond_exp_func_model,
        private$current_n,
        private$current_p,
        private$current_betaT,
        rep,
        design,
        inference,
        inference_type
      )
    },
    .result_key_from_row = function(row) {
      if (!is.list(row)) {
        stop(".result_key_from_row expected a list, but got: ", typeof(row))
      }
      private$.result_key_for_values(
        row[["response_type"]],
        row[["cond_exp_func_model"]],
        row[["n"]],
        row[["p"]],
        row[["betaT"]],
        row[["rep"]],
        row[["design"]],
        row[["inference"]],
        row[["inference_type"]]
      )
    },
    .result_metadata_dt = function(rep, design, inference, inference_type) {
      data.table::data.table(
        response_type = private$current_response_type,
        cond_exp_func_model = private$current_cond_exp_func_model,
        n = as.integer(private$current_n),
        p = as.integer(private$current_p),
        betaT = as.numeric(private$current_betaT),
        rep = as.integer(rep),
        design = design,
        inference = inference,
        inference_type = inference_type
      )
    },
    .make_error_record = function(stage, rep, design, design_params, inference,
                                  inference_params, inference_type,
                                  inference_type_params, message,
                                  metadata = NULL) {
      list(
        response_type = private$current_response_type,
        cond_exp_func_model = private$current_cond_exp_func_model,
        n = as.integer(private$current_n),
        p = as.integer(private$current_p),
        betaT = as.numeric(private$current_betaT),
        rep = if (is.null(rep) || !length(rep)) NA_integer_ else as.integer(rep),
        design = design %||% NA_character_,
        design_params = design_params,
        inference = inference %||% NA_character_,
        inference_params = inference_params,
        inference_type = inference_type %||% NA_character_,
        inference_type_params = inference_type_params,
        stage = stage,
        error_message = as.character(message)[1L],
        metadata = metadata %||% list(),
        timestamp = as.character(Sys.time())
      )
    },
    .append_errors = function(errors) {
      if (length(errors) == 0L) return(invisible(NULL))
      private$error_log = c(private$error_log, errors)
      invisible(NULL)
    },
    .format_error_record = function(err) {
      path_parts = c(err$design, err$inference, err$inference_type)
      path_parts = path_parts[is.finite(nchar(path_parts)) & nzchar(path_parts) & !is.na(path_parts)]
      path = if (length(path_parts) == 0L) "<no path>" else paste(path_parts, collapse = " -> ")
      sprintf(
        paste0(
          "SimulationFramework stopped on error.\n",
          "stage: %s\n",
          "cell: response_type=%s, cond_exp_func_model=%s, n=%s, p=%s, betaT=%s\n",
          "rep: %s\n",
          "path: %s\n",
          "message: %s"
        ),
        err$stage,
        err$response_type,
        err$cond_exp_func_model,
        err$n,
        err$p,
        err$betaT,
        if (is.na(err$rep)) "NA" else err$rep,
        path,
        err$error_message
      )
    },
    .abort_from_error_record = function(err) {
      msg = private$.format_error_record(err)
      tryCatch(
        writeLines(c(as.character(Sys.time()), msg, ""), file.path(tempdir(), "edi_sim_crash.log")),
        error = function(e) invisible(NULL)
      )
      stop(msg, call. = FALSE)
    },
    .log_skip = function(rep, design, inference, inference_type) {
      invisible(NULL)
    },
    .valid_inference_types = function(inf_obj) {
      valid_inference_types = character(0L)
      if (is(inf_obj, "InferenceAsymp")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, c("asymp_ci", "asymp_pval"))
        )
      }
      if (is(inf_obj, "InferenceExact")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, c("exact_ci", "exact_pval"))
        )
      }
      if (is(inf_obj, "InferenceNonParamBootstrap")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, c("boot_ci", "boot_pval"))
        )
      }
      if (is(inf_obj, "InferenceRand")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, "rand_pval")
        )
        if (is(inf_obj, "InferenceRandCI") &&
            private$current_response_type %in% c("continuous", "proportion", "count")) {
          valid_inference_types = c(
            valid_inference_types,
            intersect(private$inf_types, "rand_ci")
          )
        }
      }
      valid_inference_types
    },
    .build_valid_combos_for_current_cell = function(rep_data) {
      X = rep_data$X
      y_linear_model = rep_data$y_linear_model
      sim_mode = compute_simulation_mode(
        private$custom_dgp, private$custom_replication_data_generator,
        private$custom_apply_treatment_and_noise, private$make_estimand_fn
      )
      combos = list()
      for (di in seq_along(private$design_classes)) {
        design_gen   = private$design_classes[[di]]
        design_name  = private$design_labels[[di]]
        design_extra = if (!is.null(private$design_params)) private$design_params[[di]] else list()
        des_obj = tryCatch(
          private$.build_design(design_gen, X, y_linear_model, design_extra, skip_assignment = TRUE),
          error = function(e) NULL
        )
        if (is.null(des_obj)) next
        for (ii in seq_along(private$inference_classes)) {
          inf_gen  = private$inference_classes[[ii]]
          inf_name = private$inference_labels[[ii]]
          inf_ctor_extra = private$inference_constructor_params[[ii]]
          toggle_asserts(TRUE)
          inf_obj = tryCatch(
            do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra)),
            error = function(e) NULL
          )
          toggle_asserts(FALSE)
          if (is.null(inf_obj)) next
          valid_inference_types = private$.valid_inference_types(inf_obj)
          if (length(valid_inference_types) == 0L) next
          for (it in valid_inference_types) {
            # Validate user-supplied params for this inference type early
            private$.args_for_inf_type(inf_obj, it)
            
            combos[[length(combos) + 1L]] = list(
              response_type = private$current_response_type,
              cond_exp_func_model = private$current_cond_exp_func_model,
              n = private$current_n,
              p = private$current_p,
              betaT = private$current_betaT,
              design = design_name,
              inference = inf_name,
              inference_type = it,
              simulation_mode = sim_mode
            )
          }
        }
      }
      combos
    },
    .run_single_replication_in_worker = function(rep_i, state, progress_cb = NULL, is_forked = FALSE) {
      # This runs in a worker process. It must be self-contained.
      # 1. Cap threads and nested parallelism to avoid N*M oversubscription.
      # Use loadNamespace to ensure EDI is loaded and functions are accessible,
      # because closure environment assignment is not reliably preserved across
      # mirai serialization.
      ns_edi = loadNamespace("EDI")
      ns_edi$set_package_threads(1L)
      if (is_forked) {
        ns_edi$unset_num_cores() # Clear any inherited clusters
      }
      
      error_records = list()
      make_error = function(stage, design = NA_character_, design_params = NULL,
                            inference = NA_character_, inference_params = NULL,
                            inference_type = NA_character_, inference_type_params = NULL,
                            message, metadata = NULL) {
        list(
          response_type = state$response_type,
          cond_exp_func_model = state$cond_exp_func_model,
          n = as.integer(state$n),
          p = as.integer(state$p),
          betaT = as.numeric(state$betaT),
          rep = as.integer(rep_i),
          design = design %||% NA_character_,
          design_params = design_params,
          inference = inference %||% NA_character_,
          inference_params = inference_params,
          inference_type = inference_type %||% NA_character_,
          inference_type_params = inference_type_params,
          stage = stage,
          error_message = as.character(message)[1L],
          metadata = metadata %||% list(),
          timestamp = as.character(Sys.time())
        )
      }
      handle_error = function(err) {
        error_records[[length(error_records) + 1L]] <<- err
        if (isTRUE(state$stop_on_error)) {
          return(list(
            results_dt = if (length(results) > 0L) data.table::rbindlist(results) else NULL,
            skipped_count = skipped_count,
            errors = error_records,
            fatal_error = err
          ))
        }
        NULL
      }
      state_for_rep = function() {
        state$rep = rep_i
        state
      }
      validate_custom_replication_data = function(data) {
        if (!is.list(data))
          stop("custom_replication_data_generator must return a list")
        if (is.null(data$X) || is.null(data$y_linear_model))
          stop("custom_replication_data_generator must return 'X' and 'y_linear_model'")
        X = data$X
        if (!is.data.frame(X)) X = as.data.frame(X)
        y_linear_model = as.numeric(data$y_linear_model)
        if (nrow(X) != state$n)
          stop("custom replication data returned X with ", nrow(X), " rows; expected ", state$n)
        if (length(y_linear_model) != state$n)
          stop("custom replication data returned y_linear_model of length ", length(y_linear_model), "; expected ", state$n)
        data$X = X
        data$y_linear_model = y_linear_model
        data
      }
      apply_treatment_and_noise = function(y_linear_model, w, rep_data = NULL) {
        if (!is.null(state$custom_apply_treatment_and_noise)) {
          fn = state$custom_apply_treatment_and_noise
          fn_args = names(formals(fn))
          out = if ("rep_data" %in% fn_args || "..." %in% fn_args) {
            fn(y_linear_model, w, rep_data, state_for_rep())
          } else {
            fn(y_linear_model, w, state_for_rep())
          }
          if (!is.list(out) || is.null(out$y) || is.null(out$dead)) {
            stop("custom_apply_treatment_and_noise must return a list with 'y' and 'dead'")
          }
          return(list(y = out$y, dead = out$dead))
        }
        apply_treatment_and_noise_cpp(
          y_linear_model, w,
          state$response_type, state$betaT,
          state$sd_noise, state$prob_censoring,
          state$n_ordinal_levels,
          phi_proportion = state$phi_proportion,
          k_survival = state$k_survival,
          incidence_clamp = state$incidence_clamp,
          proportion_clamp = state$proportion_clamp,
          count_clamp = state$count_clamp,
          survival_clamp = state$survival_clamp
        )
      }
      # 1. Generate data
      data = tryCatch(
        if (!is.null(state$custom_replication_data_generator)) {
          validate_custom_replication_data(state$custom_replication_data_generator(state_for_rep(), rep_i))
        } else {
          generate_covariate_dataset(
            n                    = state$n,
            p                    = state$p,
            cond_exp_func_model  = state$cond_exp_func_model,
            norm_sq_beta_vec     = state$norm_sq_beta_vec,
            X_mat                = state$shared_X %||% state$X_mat,
            cov_draw_method      = if (is.null(state$shared_X) && is.null(state$X_mat)) state$cov_draw_method else NULL,
            cov_draw_method_args = state$cov_draw_method_args
          )
        },
        error = function(e) {
          fatal = handle_error(make_error(
            stage = "data_generation",
            message = conditionMessage(e),
            metadata = list(condition_class = class(e))
          ))
          if (!is.null(fatal)) return(fatal)
          NULL
        }
      )
      if (is.list(data) && !is.null(data$fatal_error)) return(data)
      if (is.null(data)) {
        return(list(results_dt = NULL, skipped_count = 0L, errors = error_records, fatal_error = NULL))
      }
      y_linear_model = if (!is.null(state$custom_replication_data_generator)) {
        as.numeric(data$y_linear_model)
      } else {
        as.numeric(data$y_cont - mean(data$y_cont))
      }
      X = data$X
      rep_data = data  # full CRDG output; passed to custom hooks so extra fields (e.g. frailty draws) survive

      # True ATE calculation logic (extracted from compute_true_mean_diff_ate)
      clamp = function(x, lo, hi) {
        pmin(hi, pmax(lo, x))
      }
      # Default per-rep estimand (make_estimand_fn override applied per-design inside the loop)
      true_mean_diff_ate = switch(state$response_type,
        continuous = state$betaT,
        incidence = {
          p_t = clamp(stats::plogis(y_linear_model + state$betaT), state$incidence_clamp, 1 - state$incidence_clamp)
          p_c = clamp(stats::plogis(y_linear_model), state$incidence_clamp, 1 - state$incidence_clamp)
          mean(p_t - p_c)
        },
        proportion = {
          p_t = clamp(stats::plogis(y_linear_model + state$betaT), state$proportion_clamp, 1 - state$proportion_clamp)
          p_c = clamp(stats::plogis(y_linear_model), state$proportion_clamp, 1 - state$proportion_clamp)
          mean(p_t - p_c)
        },
        count = {
          r_t = clamp(exp(y_linear_model + state$betaT), state$count_clamp, Inf)
          r_c = clamp(exp(y_linear_model), state$count_clamp, Inf)
          mean(r_t - r_c)
        },
        NA_real_
      )
      results = list()
      result_keys = character()
      skipped_count = 0L
      # 2. Design and Inference loop
      for (di in seq_along(state$design_classes)) {
        design_gen   = state$design_classes[[di]]
        design_name  = state$design_labels[[di]]
        design_extra = if (!is.null(state$design_params)) state$design_params[[di]] else list()
        # Auto-inject covariate-dependent params (mirrors .build_design)
        init_fn_w = get_r6_init_fn(design_gen)
        if (!is.null(init_fn_w)) {
          fn_formals_w = names(formals(init_fn_w))
          x_names_w    = names(X)
        if ("strata_cols" %in% fn_formals_w &&
            !"strata_cols" %in% names(design_extra) &&
            !identical(design_gen$classname, "DesignFixedBlocking"))
          design_extra$strata_cols = x_names_w[1L]
          if ("cluster_col" %in% fn_formals_w && !"cluster_col" %in% names(design_extra))
            design_extra$cluster_col = x_names_w[min(2L, length(x_names_w))]
          if ("factors"     %in% fn_formals_w && !"factors"     %in% names(design_extra))
            design_extra$factors = list(treatment = 2L)
        }
        # Build design (extracted from .build_design)
        obs_out = NULL  # reset per-design; populated in Mode 3 branch below
        des_obj = tryCatch({
          d = do.call(design_gen$new, c(list(response_type = state$response_type, n = state$n), design_extra))
          if (!is.null(state$custom_dgp)) {
            # Mode 3: load complete dataset from observational DGP; design is a data container only
            if (grepl("SeqOneByOne", design_gen$classname, fixed = FALSE))
              stop("custom_dgp requires a fixed design class; '",
                   design_gen$classname, "' is a sequential design")
            dgp_fn = state$custom_dgp
            dgp_fn_args = names(formals(dgp_fn))
            obs_out = if ("state" %in% dgp_fn_args || "..." %in% dgp_fn_args) {
              dgp_fn(state$n, state$p, rep_i, state_for_rep())
            } else {
              dgp_fn(state$n, state$p, rep_i)
            }
            if (!is.list(obs_out) || is.null(obs_out$X) || is.null(obs_out$w) || is.null(obs_out$y))
              stop("custom_dgp must return a list with 'X', 'w', and 'y'")
            if (!is.data.frame(obs_out$X)) obs_out$X = as.data.frame(obs_out$X)
            if (nrow(obs_out$X) != state$n)
              stop("custom_dgp returned X with ", nrow(obs_out$X), " rows; expected ", state$n)
            d$add_all_subjects_to_experiment(obs_out$X)
            d$assign_w_to_all_subjects(2L * as.integer(obs_out$w) - 1L)
            dead_obs = obs_out$dead %||% rep(1L, state$n)
            d$add_all_subject_responses(obs_out$y, dead_obs)
          } else if (inherits(d, "DesignSeqOneByOne")) {
            for (t in seq_len(state$n)) {
              w_t = d$add_one_subject_to_experiment_and_assign(X[t, , drop = FALSE])
              out = apply_treatment_and_noise(y_linear_model[t], w_t, rep_data)
              d$add_one_subject_response(t, out$y, out$dead)
            }
          } else {
            # Check for precomp_w BEFORE building the design so we can skip the
            # expensive add_all_subjects_to_experiment (ILP/matching) when cached.
            precomp_w = NULL
            if (!is.null(state$design_w_cache) && !is.null(state$design_w_cache[[design_name]])) {
              cache = state$design_w_cache[[design_name]]
              col_j = cache$rep_to_col[as.character(rep_i)]
              if (!is.na(col_j)) precomp_w = cache$ws[, col_j]
            }
            if (!is.null(precomp_w)) {
              # Fast path: bypass ILP/matching — populate minimum state directly.
              priv = d$.__enclos_env__$private
              priv$Xraw = data.table::as.data.table(X)
              priv$Ximp = data.table::copy(priv$Xraw)
              priv$X    = as.matrix(X)
              priv$t    = as.integer(state$n)
              # Block / strata IDs needed by CMH, ExtendedRobins, etc.
              if (inherits(d, "DesignFixedBinaryMatch") &&
                  exists("ensure_matching_structure_computed", envir = priv, inherits = FALSE)) {
                priv$ensure_matching_structure_computed()
              } else if (inherits(d, "DesignFixedOptimalBlocks") &&
                         exists("get_or_compute_block_ids", envir = priv, inherits = FALSE)) {
                priv$m = as.integer(priv$get_or_compute_block_ids())
              } else if (!is.null(priv$strata_cols) && length(priv$strata_cols) > 0L &&
                         exists("get_strata_keys", envir = priv, inherits = FALSE)) {
                strata_keys = priv$get_strata_keys()
                if (length(strata_keys) == state$n) priv$m = match(strata_keys, unique(strata_keys))
              }
            } else {
              d$add_all_subjects_to_experiment(X)
            }
            # Inject the full Nrep W matrix for CMH SE estimation (K = Nrep).
            # Covers both pregen designs and non-blocking designs with CMH inference.
            if (!is.null(state$design_w_cache) && !is.null(state$design_w_cache[[design_name]])) {
              d$inject_cmh_se_w_mat(state$design_w_cache[[design_name]]$ws)
            }
            d$assign_w_to_all_subjects(precomp_w)
            w = d$get_w()
            out = apply_treatment_and_noise(y_linear_model, w, rep_data)
            d$add_all_subject_responses(out$y, out$dead)
          }
          d
        }, error = function(e) {
          fatal = handle_error(make_error(
            stage = "design_build",
            design = design_name,
            design_params = design_extra,
            message = conditionMessage(e),
            metadata = list(
              design_class = design_gen$classname,
              condition_class = class(e)
            )
          ))
          if (!is.null(fatal)) return(fatal)
          NULL
        })
        if (is.list(des_obj) && !is.null(des_obj$fatal_error)) return(des_obj)

        if (is.null(des_obj)) next
        # Per-design true estimand: make_estimand_fn (Mode 2) or custom_dgp$true_estimand (Mode 3)
        # These override the pre-rep default and apply to ALL inference class types.
        if (!is.null(state$make_estimand_fn)) {
          fn_cte = state$make_estimand_fn(state$betaT)
          fn_cte_args = names(formals(fn_cte))
          true_mean_diff_ate = as.numeric(
            if ("X" %in% fn_cte_args || "w" %in% fn_cte_args || "rep_data" %in% fn_cte_args || "..." %in% fn_cte_args) {
              fn_cte(y_linear_model, X, des_obj$get_w(), rep_data, state_for_rep())
            } else {
              fn_cte(y_linear_model, state_for_rep())
            }
          )[1L]
        } else if (!is.null(obs_out) && !is.null(obs_out$true_estimand)) {
          true_mean_diff_ate = as.numeric(obs_out$true_estimand)[1L]
        }
        for (ii in seq_along(state$inference_classes)) {
          inf_gen  = state$inference_classes[[ii]]
          inf_name = state$inference_labels[[ii]]
          inf_ctor_extra = state$inference_ctor_params[[ii]]
          # Skip combos that plan-validation already rejected so a re-fired
          # assertion (e.g. "requires a blocking design") never becomes a
          # fatal stop_on_error in the worker.
          if (!is.null(state$valid_combo_keys) && length(state$valid_combo_keys) > 0L) {
            any_valid_type = any(startsWith(state$valid_combo_keys,
                                            paste0(design_name, "|", inf_name, "|")))
            if (!any_valid_type) next
          }
          inf_obj = tryCatch({
            do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra))
          }, error = function(e) {
            fatal = handle_error(make_error(
              stage = "inference_initialize",
              design = design_name,
              design_params = design_extra,
              inference = inf_name,
              inference_params = inf_ctor_extra,
              message = conditionMessage(e),
              metadata = list(
                design_class = design_gen$classname,
                inference_class = inf_gen$classname,
                condition_class = class(e)
              )
            ))
            if (!is.null(fatal)) return(fatal)
            NULL
          })
          if (is.list(inf_obj) && !is.null(inf_obj$fatal_error)) return(inf_obj)
          
          if (is.null(inf_obj)) next
          is_mean_diff = is(inf_obj, "InferenceAllSimpleMeanDiff") ||
                         is(inf_obj, "InferenceIncidWald") ||
                         is(inf_obj, "InferenceIncidCMH") ||
                         is(inf_obj, "InferenceIncidExtendedRobins") ||
                         is(inf_obj, "InferenceIncidRiskDiff") ||
                         is(inf_obj, "InferenceIncidMiettinenNurminenRiskDiff") ||
                         is(inf_obj, "InferenceIncidNewcombeRiskDiff") ||
                         is(inf_obj, "InferenceIncidKKNewcombeRiskDiff") ||
                         is(inf_obj, "InferenceIncidKKGCompRiskDiff") ||
                         is(inf_obj, "InferencePropGCompMeanDiff")
          te = if (!is.null(state$make_estimand_fn) || !is.null(state$custom_dgp)) {
            true_mean_diff_ate  # custom estimand always wins regardless of inference class type
          } else if (is_mean_diff) {
            true_mean_diff_ate
          } else {
            state$betaT
          }
          # .valid_inference_types logic (extracted)
          valid_inference_types = character(0)
          if (is(inf_obj, "InferenceAsymp")) 
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, c("asymp_ci", "asymp_pval")))
          if (is(inf_obj, "InferenceExact"))
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, c("exact_ci", "exact_pval")))
          if (is(inf_obj, "InferenceNonParamBootstrap"))
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, c("boot_ci", "boot_pval")))
          if (is(inf_obj, "InferenceRand")) {
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, "rand_pval"))
            if (is(inf_obj, "InferenceRandCI") && state$response_type %in% c("continuous", "proportion", "count"))
              valid_inference_types = c(valid_inference_types, intersect(state$inf_types, "rand_ci"))
          }
          # Result key and pending check
          pending_inference_types = valid_inference_types[!check_in_result_key_store_cpp(
            rep(state$response_type, length(valid_inference_types)),
            rep(state$cond_exp_func_model, length(valid_inference_types)),
            rep(state$n, length(valid_inference_types)),
            rep(state$p, length(valid_inference_types)),
            rep(state$betaT, length(valid_inference_types)),
            rep(rep_i, length(valid_inference_types)),
            rep(design_name, length(valid_inference_types)),
            rep(inf_name, length(valid_inference_types)),
            valid_inference_types
          )]
          
          skipped_count = skipped_count + (length(valid_inference_types) - length(pending_inference_types))
          if (length(pending_inference_types) == 0L) {
            # Advance progress for skipped types
            if (!is.null(progress_cb)) {
               for (k in seq_along(valid_inference_types)) progress_cb()
            }
            next
          }
          est = tryCatch({
            v = inf_obj$compute_estimate()
            if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v)[1L]
          }, error = function(e) {
            fatal = handle_error(make_error(
              stage = "estimate",
              design = design_name,
              design_params = design_extra,
              inference = inf_name,
              inference_params = inf_ctor_extra,
              message = conditionMessage(e),
              metadata = list(
                design_class = design_gen$classname,
                inference_class = inf_gen$classname,
                pending_inference_types = pending_inference_types,
                condition_class = class(e)
              )
            ))
            if (!is.null(fatal)) return(fatal)
            NA_real_
          })
          if (is.list(est) && !is.null(est$fatal_error)) return(est)
          # Helper to merge user-supplied params with defaults
          get_args = function(type, defaults = list()) {
            user_args = state$inference_type_params[[type]]
            if (is.null(user_args)) user_args = list()
            modifyList(defaults, user_args)
          }
          # helper for recording in worker
          local_record = function(type, ci, pval) {
            ci2 = if (length(ci) >= 2L) as.numeric(ci[1:2]) else c(NA_real_, NA_real_)
            if (all(is.finite(ci2)) && ci2[1L] > ci2[2L]) ci2 = rev(ci2)
            sim_mode = compute_simulation_mode(
              state$custom_dgp, state$custom_replication_data_generator,
              state$custom_apply_treatment_and_noise, state$make_estimand_fn
            )
            results[[length(results) + 1L]] <<- list(
              response_type = state$response_type,
              rep           = rep_i,
              cond_exp_func_model = state$cond_exp_func_model,
              n             = as.integer(state$n),
              p             = as.integer(state$p),
              betaT         = as.numeric(state$betaT),
              design        = design_name,
              inference     = inf_name,
              inference_type = type,
              estimate      = if (is.null(est) || !is.finite(est)) NA_real_ else as.numeric(est),
              ci_lo          = ci2[1L],
              ci_hi          = ci2[2L],
              pval          = if (is.null(pval) || length(pval) == 0L || !is.finite(pval[1L]))
                                NA_real_ else as.numeric(pval[1L]),
              true_estimand = as.numeric(te),
              simulation_mode = sim_mode
            )
            if (!is.null(progress_cb)) progress_cb()
          }
          # Advance progress for already-present results (cached)
          if (!is.null(progress_cb) && length(valid_inference_types) > length(pending_inference_types)) {
             already_done = length(valid_inference_types) - length(pending_inference_types)
             for (k in seq_len(already_done)) progress_cb()
          }
          # Inference execution logic
          if (is(inf_obj, "InferenceAsymp") && any(c("asymp_ci", "asymp_pval") %in% pending_inference_types)) {
            if ("asymp_pval" %in% pending_inference_types) {
              args = get_args("asymp_pval")
              pval_a = tryCatch(do.call(inf_obj$compute_asymp_two_sided_pval, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "asymp_pval",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_asymp_two_sided_pval",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                NA_real_
              })
              if (is.list(pval_a) && !is.null(pval_a$fatal_error)) return(pval_a)
              local_record("asymp_pval", c(NA_real_, NA_real_), pval_a)
            }
            if ("asymp_ci" %in% pending_inference_types) {
              args = get_args("asymp_ci", list(alpha = state$alpha))
              ci_a = tryCatch(do.call(inf_obj$compute_asymp_confidence_interval, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "asymp_ci",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_asymp_confidence_interval",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                c(NA_real_, NA_real_)
              })
              if (is.list(ci_a) && !is.null(ci_a$fatal_error)) return(ci_a)
              local_record("asymp_ci", ci_a, NA_real_)
            }
          }
          
          if (is(inf_obj, "InferenceExact") && any(c("exact_ci", "exact_pval") %in% pending_inference_types)) {
            if ("exact_pval" %in% pending_inference_types) {
              args = get_args("exact_pval")
              pval_e = tryCatch(do.call(inf_obj$compute_exact_two_sided_pval_for_treatment_effect, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "exact_pval",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_exact_two_sided_pval_for_treatment_effect",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                NA_real_
              })
              if (is.list(pval_e) && !is.null(pval_e$fatal_error)) return(pval_e)
              local_record("exact_pval", c(NA_real_, NA_real_), pval_e)
            }
            if ("exact_ci" %in% pending_inference_types) {
              args = get_args("exact_ci", list(alpha = state$alpha))
              ci_e = tryCatch(do.call(inf_obj$compute_exact_confidence_interval, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "exact_ci",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_exact_confidence_interval",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                c(NA_real_, NA_real_)
              })
              if (is.list(ci_e) && !is.null(ci_e$fatal_error)) return(ci_e)
              local_record("exact_ci", ci_e, NA_real_)
            }
          }
          
          if (is(inf_obj, "InferenceNonParamBootstrap") && any(c("boot_ci", "boot_pval") %in% pending_inference_types)) {
            if ("boot_pval" %in% pending_inference_types) {
              args = get_args("boot_pval", list(B = state$B_boot, na.rm = TRUE))
              pval_b = tryCatch(do.call(inf_obj$compute_bootstrap_two_sided_pval, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "boot_pval",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_bootstrap_two_sided_pval",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                NA_real_
              })
              if (is.list(pval_b) && !is.null(pval_b$fatal_error)) return(pval_b)
              local_record("boot_pval", c(NA_real_, NA_real_), pval_b)
            }
            if ("boot_ci" %in% pending_inference_types) {
              args = get_args("boot_ci", list(B = state$B_boot, alpha = state$alpha, na.rm = TRUE, show_progress = FALSE))
              ci_b = tryCatch(do.call(inf_obj$compute_bootstrap_confidence_interval, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "boot_ci",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_bootstrap_confidence_interval",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                c(NA_real_, NA_real_)
              })
              if (is.list(ci_b) && !is.null(ci_b$fatal_error)) return(ci_b)
              local_record("boot_ci", ci_b, NA_real_)
            }
          }
          
          if (is(inf_obj, "InferenceRand") && any(c("rand_ci", "rand_pval") %in% pending_inference_types)) {
            if ("rand_pval" %in% pending_inference_types) {
              args = get_args("rand_pval", list(r = state$r_rand, na.rm = TRUE, show_progress = FALSE))
              pval_r = tryCatch(do.call(inf_obj$compute_rand_two_sided_pval, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "rand_pval",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_rand_two_sided_pval",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                NA_real_
              })
              if (is.list(pval_r) && !is.null(pval_r$fatal_error)) return(pval_r)
              local_record("rand_pval", c(NA_real_, NA_real_), pval_r)
            }
            if ("rand_ci" %in% pending_inference_types && is(inf_obj, "InferenceRandCI") && state$response_type %in% c("continuous", "proportion", "count")) {
              args = get_args("rand_ci", list(r = state$r_rand, alpha = state$alpha, pval_epsilon = state$pval_epsilon, show_progress = FALSE))
              ci_r = tryCatch(do.call(inf_obj$compute_rand_confidence_interval, args), error = function(e) {
                fatal = handle_error(make_error(
                  stage = "inference_call",
                  design = design_name,
                  design_params = design_extra,
                  inference = inf_name,
                  inference_params = inf_ctor_extra,
                  inference_type = "rand_ci",
                  inference_type_params = args,
                  message = conditionMessage(e),
                  metadata = list(
                    method = "compute_rand_confidence_interval",
                    design_class = design_gen$classname,
                    inference_class = inf_gen$classname,
                    condition_class = class(e)
                  )
                ))
                if (!is.null(fatal)) return(fatal)
                c(NA_real_, NA_real_)
              })
              if (is.list(ci_r) && !is.null(ci_r$fatal_error)) return(ci_r)
              local_record("rand_ci", ci_r, NA_real_)
            }
          }
        }
      }
      # Return results as a data.table for efficiency in master loop
      list(
        results_dt = if (length(results) > 0L) data.table::rbindlist(results) else NULL,
        skipped_count = skipped_count,
        errors = error_records,
        fatal_error = NULL
      )
    },
    .advance_progress = function() {
      private$current_task_in_rep_idx = private$current_task_in_rep_idx + 1L
      if (!isTRUE(private$verbose)) return(invisible(NULL))
      
      if (isTRUE(private$use_progress_bar)) {
        private$.draw_simulation_progress_bars()
      }
      invisible(NULL)
    },
    .print_plan_summary = function(planned_combos_list) {
      n_cells = length(planned_combos_list)
      cat("Simulation Plan Summary:\n", file = stderr())
      
      # Group by response_type to keep it concise
      rt_summaries = list()
      for (cell_idx in seq_len(n_cells)) {
        rt = private$param_grid$response_type[[cell_idx]]
        combos = planned_combos_list[[cell_idx]]
        design_names = unique(vapply(combos, `[[`, "", "design"))
        design_names = gsub("Design", "", design_names)
        inference_names = unique(vapply(combos, `[[`, "", "inference"))
        inference_names = gsub("Inference", "", inference_names)
        inference_names = gsub("^(Contin|Count|Incidence|Incid|Prop|Survival|Ordinal|All)", "", inference_names)
        if (is.null(rt_summaries[[rt]])) {
          rt_summaries[[rt]] = list(
            designs = design_names,
            inferences = inference_names,
            n_tasks = length(combos)
          )
        } else {
          rt_summaries[[rt]]$designs = unique(c(rt_summaries[[rt]]$designs, design_names))
          rt_summaries[[rt]]$inferences = unique(c(rt_summaries[[rt]]$inferences, inference_names))
          rt_summaries[[rt]]$n_tasks = rt_summaries[[rt]]$n_tasks + length(combos)
        }
      }
      response_types = names(rt_summaries)
      if (length(response_types) == 1L) {
        s = rt_summaries[[response_types[[1L]]]]
        cat(sprintf("  Designs (%d): %s\n", length(s$designs), paste(s$designs, collapse = ", ")), file = stderr())
        cat(sprintf("  Inferences (%d): %s\n", length(s$inferences), paste(s$inferences, collapse = ", ")), file = stderr())
      } else {
        for (rt in response_types) {
          s = rt_summaries[[rt]]
          cat(sprintf("  - Response Type: %s\n", rt), file = stderr())
          cat(sprintf("    Designs (%d): %s\n", length(s$designs), paste(s$designs, collapse = ", ")), file = stderr())
          cat(sprintf("    Inferences (%d): %s\n", length(s$inferences), paste(s$inferences, collapse = ", ")), file = stderr())
        }
      }
      cat("\n", file = stderr())
    },
    .draw_progress = function() {
      if (!isTRUE(private$verbose)) return(invisible(NULL))
      if (isTRUE(private$use_progress_bar)) {
        private$.draw_simulation_progress_bars()
      } else if (private$progress_log_interval > 0L &&
                 (private$progress_count %% private$progress_log_interval == 0L ||
                  private$progress_count == private$progress_total)) {
        message(sprintf("Completed %d / %d runs", private$progress_count, private$progress_total))
      }
    },
    .draw_labeled_progress_bar = function(label, prop) {
      width = getOption("width", 80L)
      if (is.null(width) || width < 80L) width = 80L
      
      bar_width = width - nchar(label) - 10L
      if (bar_width < 10L) bar_width = 10L
      
      make_bar = function(p, b_width) {
        pct_str = sprintf(" %3d%% ", floor(p * 100))
        n_pct = nchar(pct_str)
        fill = floor(p * b_width)
        full_bar = paste0(strrep("=", fill), strrep(" ", b_width - fill))
        if (b_width >= n_pct) {
          start_pos = (b_width - n_pct) %/% 2 + 1
          substr(full_bar, start_pos, start_pos + n_pct - 1) = pct_str
        }
        sprintf("[%s]", full_bar)
      }
      
      content = sprintf("%s %s", label, make_bar(prop, bar_width))
      content = substr(content, 1, width - 1L)
      n_pad = (width - 1L) - nchar(content)
      if (n_pad > 0L) content = paste0(content, strrep(" ", n_pad))
      cat(paste0("\r", content), file = stderr())
      if (exists("flush.console")) utils::flush.console()
    },
    .message_stderr = function(msg) {
      if (isTRUE(private$progress_bar_drawn)) {
        cat(sprintf("\033[%dA\r\033[J", private$n_progress_lines), file = stderr())
        private$progress_bar_drawn = NULL
      }
      cat(msg, file = stderr())
      if (exists("flush.console")) utils::flush.console()
    },
    .draw_simulation_progress_bars = function() {
      now = as.numeric(Sys.time())
      # Throttle: only redraw if 100ms have passed OR we are at 100%
      is_done = private$current_rep_idx == private$Nrep &&
                private$current_cell_idx == private$total_cells &&
                private$current_task_in_rep_idx == private$tasks_per_rep
      if (!is_done && (now - private$last_progress_draw_time) < 0.1) return(invisible(NULL))
      private$last_progress_draw_time = now
      # Force a narrow width to prevent wrapping which breaks ANSI sequences.
      width = 60L
      task_total_display = max(private$tasks_per_rep, private$current_task_in_rep_idx)
      # Proportions (clamped to [0,1]).
      # Outer loop is rep; inner is cell; innermost is task.
      task_in_progress_prop = max(0, min(1, if (task_total_display > 0) private$current_task_in_rep_idx / task_total_display else 0))
      cell_in_progress_prop = max(0, min(1, if (private$total_cells > 0) (max(0, private$current_cell_idx - 1) + task_in_progress_prop) / private$total_cells else 0))
      rep_in_progress_prop  = max(0, min(1, if (private$Nrep > 0) (max(0, private$current_rep_idx - 1) + cell_in_progress_prop) / private$Nrep else 0))
      # overall_prop accounts for work done before this session started.
      # Force to 1.0 when all reps/cells/tasks are complete to avoid showing < 100%
      # due to deduplication reducing progress_count increments.
      overall_prop = if (is_done) 1.0 else if (private$progress_total > 0L) private$progress_count / private$progress_total else rep_in_progress_prop

      # ETA: rep-based once ≥1 rep has finished; task-throughput rough estimate before then.
      .fmt_secs = function(secs) {
        d = floor(secs / 86400); secs = secs %% 86400
        h = floor(secs / 3600);  secs = secs %% 3600
        m = floor(secs / 60);    s = round(secs %% 60)
        parts = character()
        if (d > 0) parts = c(parts, paste0(d, "d"))
        if (h > 0) parts = c(parts, paste0(h, "h"))
        if (m > 0) parts = c(parts, paste0(m, "m"))
        parts = c(parts, paste0(s, "s"))
        paste(parts, collapse = " ")
      }
      n_done = private$rep_elapsed_idx
      eta_str = if (overall_prop >= 0.9999) {
        "Status: Completed."
      } else if (n_done > 0L && private$n_active_reps_total > 0L) {
        # Rep-based ETA: rolling mean of last 5 reps, minus time already spent in
        # the current in-flight rep so the countdown ticks every 100 ms.
        secs_per_active_rep  = mean(private$rep_elapsed_times[seq(max(1L, n_done - 4L), n_done)])
        n_remaining          = max(0L, private$n_active_reps_total - n_done)
        elapsed_in_rep       = if (!is.null(private$rep_start_capture)) max(0, now - private$rep_start_capture) else 0
        elapsed_in_rep       = min(elapsed_in_rep, secs_per_active_rep)
        remaining            = max(0, secs_per_active_rep * n_remaining - elapsed_in_rep)
        paste0("Time Left: ", .fmt_secs(remaining))
      } else {
        # Cell-throughput fallback — counts worker cell completions (not result rows),
        # so it is unaffected by deduplication that suppresses progress_count advances
        # when resuming a nearly-finished run.
        # Require at least one full parallel chunk before measuring throughput:
        # mid-chunk, only the fastest workers have reported, so the apparent rate
        # (wu_done / session_elapsed) is num_cores-fold lower than steady-state.
        wu_done = private$session_work_units_done
        chunk_threshold = max(1L, min(private$num_cores, private$n_total_active_work_units))
        session_elapsed = if (!is.null(private$session_start_time)) now - private$session_start_time else 0
        if (wu_done >= chunk_threshold && session_elapsed > 0.5) {
          remaining_wu = max(0L, private$n_total_active_work_units - wu_done)
          remaining = session_elapsed * remaining_wu / wu_done
          paste0("~Time Left: ", .fmt_secs(remaining))
        } else {
          "Status: Estimating..."
        }
      }
      make_bar_line = function(label, prop, b_width, digits = 0, label_width = 25) {
        padded_label = sprintf("%-*s", label_width, substr(label, 1, label_width))
        label_len = nchar(padded_label)
        bar_available = b_width - label_len - 2L
        if (bar_available < 10L) return(padded_label)
        fill = floor(prop * bar_available)
        fill = max(0, min(bar_available, fill))
        if (digits == 0) {
          pct_str = sprintf(" %d%% ", floor(prop * 100))
        } else {
          pct_str = sprintf(" %.1f%% ", prop * 100)
        }
        n_pct = nchar(pct_str)
        full_bar = paste0(strrep("=", fill), strrep(" ", bar_available - fill))
        if (bar_available >= n_pct) {
           start_pos = (bar_available - n_pct) %/% 2 + 1
           substr(full_bar, start_pos, start_pos + n_pct - 1) = pct_str
        }
        sprintf("%s[%s]", padded_label, full_bar)
      }
      total_elapsed = private$prior_elapsed_secs + (now - private$simulation_start_time)
      line5 = paste0("Time elapsed: ", .fmt_secs(total_elapsed))
      line1 = make_bar_line(sprintf("Rep %d/%d Overall", private$current_rep_idx, private$Nrep), overall_prop, width, digits = 1)
      # Use \r and \033[2K on EACH line for maximum robustness.
      # Under fork parallelism workers return whole reps at once, so DGP/Task bars
      # always read 100% and are suppressed; only the Rep bar is meaningful.
      if (private$n_progress_lines == 5L) {
        line2 = make_bar_line(sprintf("DGP %d/%d", private$current_cell_idx, private$total_cells), cell_in_progress_prop, width)
        line3 = make_bar_line(sprintf("%s %d/%d", private$current_task_label, private$current_task_in_rep_idx, task_total_display), task_in_progress_prop, width)
        output_block = paste0(
          "\r\033[2K", eta_str, "\n",
          "\r\033[2K", line1, "\n",
          "\r\033[2K", line2, "\n",
          "\r\033[2K", line3, "\n",
          "\r\033[2K", line5, "\n"
        )
      } else {
        output_block = paste0(
          "\r\033[2K", eta_str, "\n",
          "\r\033[2K", line1, "\n",
          "\r\033[2K", line5, "\n"
        )
      }
      if (is.null(private$progress_bar_drawn)) {
         cat(output_block, sep = "", file = stderr())
         private$progress_bar_drawn = TRUE
      } else {
         cat(sprintf("\033[%dA", private$n_progress_lines), output_block, sep = "", file = stderr())
      }
      if (exists("flush.console")) utils::flush.console()
    },
    .append_result_row_to_file = function(row) {
      format = private$.results_file_format(private$results_filename)
      if (identical(format, "csv")) {
        file_exists = file.exists(private$results_filename)
        data.table::fwrite(row, private$results_filename, append = file_exists, col.names = !file_exists)
        return(invisible(NULL))
      }
      if (!identical(format, "csv.bz2")) {
        stop("Unsupported results file format: ", private$results_filename)
      }
      staging_filename = private$.results_staging_filename()
      staging_exists = file.exists(staging_filename)
      data.table::fwrite(row, staging_filename, append = staging_exists, col.names = !staging_exists)
      invisible(NULL)
    },
    .sync_results_bz2_from_staging = function(staging_filename = private$.results_staging_filename()) {
      private$.message_stderr("Compressing results into a bz2 file...\n")
      if (!file.exists(staging_filename)) {
        stop("Cannot update compressed results because staging CSV is missing: ", staging_filename)
      }
      results_path = private$.results_output_path()
      results_dir = dirname(results_path)
      if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
      }
      compressed_tmpfile = tempfile(
        pattern = paste0(sub("\\.csv\\.bz2$", "", basename(private$results_filename), ignore.case = TRUE), "_"),
        tmpdir = results_dir,
        fileext = ".csv.bz2"
      )
      on.exit(if (file.exists(compressed_tmpfile)) unlink(compressed_tmpfile), add = TRUE)
      input_con = file(staging_filename, open = "rb")
      output_con = bzfile(compressed_tmpfile, open = "wb")
      on.exit(try(close(input_con), silent = TRUE), add = TRUE)
      on.exit(try(close(output_con), silent = TRUE), add = TRUE)
      private$.copy_binary_stream(input_con, output_con)
      close(output_con)
      close(input_con)
      if (file.exists(results_path))
        unlink(results_path)
      if (!file.rename(compressed_tmpfile, results_path)) {
        stop("Failed to move temporary compressed results into place: ", private$results_filename)
      }
      invisible(NULL)
    },
    .cleanup_results_staging_file = function() {
      if (private$.results_file_format(private$results_filename) != "csv.bz2")
        return(invisible(NULL))
      staging_filename = private$.results_staging_filename()
      if (file.exists(staging_filename))
        unlink(staging_filename)
      invisible(NULL)
    },
    .ensure_staging_file_exists = function() {
      if (private$.results_file_format(private$results_filename) != "csv.bz2")
        return(invisible(NULL))
      
      staging_filename = private$.results_staging_filename()
      if (file.exists(staging_filename))
        return(invisible(NULL))
      
      results_path = private$.results_output_path()
      if (!file.exists(results_path))
        return(invisible(NULL))
      
      if (isTRUE(private$verbose)) message("Extracting existing results to staging file...")
      input_con = bzfile(results_path, open = "rb")
      output_con = file(staging_filename, open = "wb")
      on.exit(try(close(input_con), silent = TRUE), add = TRUE)
      on.exit(try(close(output_con), silent = TRUE), add = TRUE)
      private$.copy_binary_stream(input_con, output_con)
      invisible(NULL)
    },
    # Build unique per-instance design labels: "ClassName (params)" with
    # " [k]" suffix when two instances would otherwise share the same label.
    .compute_design_labels = function() {
      labels = vapply(seq_along(private$design_classes), function(di) {
        cls = private$design_classes[[di]]$classname
        ps  = private$.params_to_str(
          if (!is.null(private$design_params)) private$design_params[[di]] else NULL)
        if (nchar(ps) > 0L) paste0(cls, " (", ps, ")") else cls
      }, "")
      # Append [k] for any group of labels that are still identical
      for (lbl in unique(labels[duplicated(labels)])) {
        idx = which(labels == lbl)
        for (k in seq_along(idx)) labels[idx[k]] = paste0(lbl, " [", k, "]")
      }
      labels
    },
    .compute_inference_labels = function() {
      labels = vapply(seq_along(private$inference_classes), function(ii) {
        cls = private$inference_classes[[ii]]$classname
        ps  = private$.params_to_str(private$inference_constructor_params[[ii]])
        if (nchar(ps) > 0L) paste0(cls, " (", ps, ")") else cls
      }, "")
      for (lbl in unique(labels[duplicated(labels)])) {
        idx = which(labels == lbl)
        for (k in seq_along(idx)) labels[idx[k]] = paste0(lbl, " [", k, "]")
      }
      labels
    },
    .state_for_current_cell = function(rep = NULL) {
      list(
        response_type = private$current_response_type,
        cond_exp_func_model = private$current_cond_exp_func_model,
        n = private$current_n,
        p = private$current_p,
        betaT = private$current_betaT,
        alpha = private$alpha,
        B_boot = private$B_boot,
        r_rand = private$r_rand,
        pval_epsilon = private$pval_epsilon,
        sd_noise = private$sd_noise,
        n_ordinal_levels = private$n_ordinal_levels,
        proportion_epsilon = private$proportion_epsilon,
        phi_proportion = private$phi_proportion,
        k_survival = private$k_survival,
        incidence_clamp = private$incidence_clamp,
        proportion_clamp = private$proportion_clamp,
        count_clamp = private$count_clamp,
        survival_clamp = private$survival_clamp,
        survival_min_time = private$survival_min_time,
        count_min_rate = private$count_min_rate,
        count_shift = private$count_shift,
        norm_sq_beta_vec = private$norm_sq_beta_vec,
        prob_censoring = private$prob_censoring,
        rep = rep
      )
    },
    .validate_custom_replication_data = function(data) {
      if (!is.list(data))
        stop("custom_replication_data_generator must return a list")
      if (is.null(data$X) || is.null(data$y_linear_model))
        stop("custom_replication_data_generator must return 'X' and 'y_linear_model'")
      X = data$X
      if (!is.data.frame(X)) X = as.data.frame(X)
      y_linear_model = as.numeric(data$y_linear_model)
      if (nrow(X) != private$current_n)
        stop("custom replication data returned X with ", nrow(X), " rows; expected ", private$current_n)
      if (length(y_linear_model) != private$current_n)
        stop("custom replication data returned y_linear_model of length ", length(y_linear_model), "; expected ", private$current_n)
      data$X = X
      data$y_linear_model = y_linear_model
      data
    },
    .apply_treatment_and_noise = function(y_linear_model, w, rep_data = NULL) {
      if (!is.null(private$custom_apply_treatment_and_noise)) {
        fn = private$custom_apply_treatment_and_noise
        fna = names(formals(fn))
        out = if ("rep_data" %in% fna || "..." %in% fna) {
          fn(y_linear_model, w, rep_data, private$.state_for_current_cell())
        } else {
          fn(y_linear_model, w, private$.state_for_current_cell())
        }
        if (!is.list(out) || is.null(out$y) || is.null(out$dead)) {
          stop("custom_apply_treatment_and_noise must return a list with 'y' and 'dead'")
        }
        return(list(y = out$y, dead = out$dead))
      }
      apply_treatment_and_noise_cpp(
        y_linear_model, w,
        private$current_response_type, private$current_betaT,
        private$sd_noise, private$prob_censoring,
        private$n_ordinal_levels,
        phi_proportion = private$phi_proportion,
        k_survival = private$k_survival,
        incidence_clamp = private$incidence_clamp,
        proportion_clamp = private$proportion_clamp,
        count_clamp = private$count_clamp,
        survival_clamp = private$survival_clamp
      )
    },
    # ── Data generation ───────────────────────────────────────────────────────
    .generate_data = function() {
      if (!is.null(private$custom_replication_data_generator)) {
        return(private$.validate_custom_replication_data(
          private$custom_replication_data_generator(private$.state_for_current_cell(), NULL)
        ))
      }
      data = generate_covariate_dataset(
        n                    = private$current_n,
        p                    = private$current_p,
        cond_exp_func_model  = private$current_cond_exp_func_model,
        norm_sq_beta_vec     = private$norm_sq_beta_vec,
        X_mat                = private$X_mat,
        cov_draw_method      = private$cov_draw_method,
        cov_draw_method_args = private$cov_draw_method_args
      )
      data$y_linear_model = as.numeric(scale(data$y_cont))
      data$y_cont = NULL # SimulationFramework doesn't need the raw cont y anymore
      data
    },
    .generate_data_from_X = function(X_mat) {
      if (!is.null(private$custom_replication_data_generator)) {
        return(private$.validate_custom_replication_data(
          private$custom_replication_data_generator(private$.state_for_current_cell(), NULL)
        ))
      }
      data = generate_covariate_dataset(
        n                    = private$current_n,
        p                    = private$current_p,
        cond_exp_func_model  = private$current_cond_exp_func_model,
        norm_sq_beta_vec     = private$norm_sq_beta_vec,
        X_mat                = X_mat,
        cov_draw_method      = NULL,
        cov_draw_method_args = private$cov_draw_method_args
      )
      data$y_linear_model = as.numeric(scale(data$y_cont))
      data$y_cont = NULL # SimulationFramework doesn't need the raw cont y anymore
      data
    },
    compute_true_mean_diff_ate = function(y_linear_model, X = NULL, w = NULL, rep_data = NULL) {
      if (!is.null(private$make_estimand_fn)) {
        fn  = private$make_estimand_fn(private$current_betaT)
        fna = names(formals(fn))
        return(as.numeric(
          if ("X" %in% fna || "w" %in% fna || "rep_data" %in% fna || "..." %in% fna) {
            fn(y_linear_model, X, w, rep_data, private$.state_for_current_cell())
          } else {
            fn(y_linear_model, private$.state_for_current_cell())
          }
        )[1L])
      }
      eta_c = y_linear_model
      eta_t = y_linear_model + private$current_betaT
      clamp = function(x, lo, hi) {
        pmin(hi, pmax(lo, x))
      }
      switch(private$current_response_type,
        continuous = private$current_betaT,
        incidence = {
          p_t = clamp(stats::plogis(eta_t), private$incidence_clamp, 1 - private$incidence_clamp)
          p_c = clamp(stats::plogis(eta_c), private$incidence_clamp, 1 - private$incidence_clamp)
          mean(p_t - p_c)
        },
        proportion = {
          mu_t = clamp(stats::plogis(eta_t), private$proportion_clamp, 1 - private$proportion_clamp)
          mu_c = clamp(stats::plogis(eta_c), private$proportion_clamp, 1 - private$proportion_clamp)
          mean(mu_t - mu_c)
        },
        count = {
          mu_t = pmax(private$count_clamp, exp(eta_t))
          mu_c = pmax(private$count_clamp, exp(eta_c))
          mean(mu_t - mu_c)
        },
        survival = {
          shape_t = pmax(private$survival_clamp, exp(eta_t))
          shape_c = pmax(private$survival_clamp, exp(eta_c))
          mean_t = private$k_survival * gamma(1 + 1 / shape_t)
          mean_c = private$k_survival * gamma(1 + 1 / shape_c)
          (1 - private$prob_censoring / 2) * mean(mean_t - mean_c)
        },
        ordinal = private$compute_true_ordinal_mean_diff(eta_c, eta_t),
        private$current_betaT
      )
    },
    compute_true_ordinal_mean_diff = function(eta_c, eta_t) {
      expected_ordinal = function(eta) {
        K = private$n_ordinal_levels
        if (private$sd_noise <= 0) {
          rounded_eta = sign(eta) * floor(abs(eta) + 0.5)
          return(pmin(K, pmax(1, rounded_eta)))
        }
        sigma = private$sd_noise
        probs = matrix(0, nrow = length(eta), ncol = K)
        probs[, 1L] = stats::pnorm((1.5 - eta) / sigma)
        if (K > 2L) {
          for (k in 2L:(K - 1L)) {
            lo = (k - 0.5 - eta) / sigma
            hi = (k + 0.5 - eta) / sigma
            probs[, k] = stats::pnorm(hi) - stats::pnorm(lo)
          }
        }
        if (K > 1L) {
          probs[, K] = 1 - stats::pnorm((K - 0.5 - eta) / sigma)
        }
        as.numeric(probs %*% seq_len(K))
      }
      mean(expected_ordinal(eta_t) - expected_ordinal(eta_c))
    },
    # Instantiate design and run the full experiment (assign + observe all n).
    .build_design = function(design_gen, X, y_linear_model, design_extra, skip_assignment = FALSE, rep_data = NULL) {
      n       = private$current_n
      # Auto-inject required args that depend on the covariate matrix when the
      # user has not already supplied them via design_classes_and_params.
      init_fn = get_r6_init_fn(design_gen)
      if (!is.null(init_fn)) {
        fn_formals = names(formals(init_fn))
        x_names    = names(X)
        if ("strata_cols" %in% fn_formals &&
            !"strata_cols" %in% names(design_extra) &&
            !identical(design_gen$classname, "DesignFixedBlocking"))
          design_extra$strata_cols = x_names[1L]
        if ("cluster_col" %in% fn_formals && !"cluster_col" %in% names(design_extra))
          design_extra$cluster_col = x_names[min(2L, length(x_names))]
        if ("factors"     %in% fn_formals && !"factors"     %in% names(design_extra))
          design_extra$factors = list(treatment = 2L)
      }
      des_obj = do.call(design_gen$new, c(
        list(response_type = private$current_response_type, n = n),
        design_extra
      ))
      if (skip_assignment) {
        # Bypass heavy assignment logic during validation phase.
        # Populate the minimum state needed for inference constructors that
        # validate against completed-design metadata such as block IDs or
        # binary-match structure.
        priv = des_obj$.__enclos_env__$private
        priv$Xraw = data.table::as.data.table(X)
        priv$Ximp = data.table::copy(priv$Xraw)
        priv$X = as.matrix(X)
        priv$w = rep(c(0L, 1L), length.out = n)
        priv$y = rep(0, n)
        priv$y_original = priv$y
        priv$dead = rep(1L, n)  # 1 = uncensored; 0 would trigger "uncensored responses" asserts
        priv$t = n
        # Some fixed designs derive their blocking structure lazily from the
        # covariates. Build that structure here so validation-time assertions
        # (e.g. CMH / Extended Robins) see the same metadata as a fully
        # realized design.
        if (inherits(des_obj, "DesignFixedBinaryMatch") &&
            private$.has_private_method_on_object(des_obj, "ensure_matching_structure_computed")) {
          priv$ensure_matching_structure_computed()
        } else if (inherits(des_obj, "DesignFixedOptimalBlocks") &&
                   private$.has_private_method_on_object(des_obj, "get_or_compute_block_ids")) {
          priv$m = as.integer(priv$get_or_compute_block_ids())
        } else if (!is.null(priv$strata_cols) &&
                   length(priv$strata_cols) > 0L &&
                   private$.has_private_method_on_object(des_obj, "get_strata_keys")) {
          strata_keys = priv$get_strata_keys()
          if (length(strata_keys) == n) {
            priv$m = match(strata_keys, unique(strata_keys))
          }
        }
        return(des_obj)
      }
      if (inherits(des_obj, "DesignSeqOneByOne")) {
        # Sequential: assignment depends on prior responses so w is obtained
        # one subject at a time.  Call Rcpp with length-1 vectors per subject.
        for (t in seq_len(n)) {
          w_t = des_obj$add_one_subject_to_experiment_and_assign(
            X[t, , drop = FALSE])
          out = private$.apply_treatment_and_noise(y_linear_model[t], w_t, rep_data)
          des_obj$add_one_subject_response(t, out$y, out$dead)
        }
      } else {
        # Fixed: all assignments known upfront — vectorize across all n subjects.
        des_obj$add_all_subjects_to_experiment(X)
        des_obj$assign_w_to_all_subjects()
        w   = des_obj$get_w()
        out = private$.apply_treatment_and_noise(y_linear_model, w, rep_data)
        des_obj$add_all_subject_responses(out$y, out$dead)
      }
      des_obj
    },
    .flush_pending_to_disk = function() {
      if (length(private$pending_file_rows) == 0L) return(invisible(NULL))
      combined = data.table::rbindlist(private$pending_file_rows, use.names = TRUE, fill = TRUE)
      private$pending_file_rows = list()
      private$.append_result_row_to_file(combined)
      invisible(NULL)
    },
    # Append multiple rows to raw_results and write to disk in one go.
    .record_batch = function(rows, skipped_count = 0L, keys = NULL) {
      is_dt = data.table::is.data.table(rows)
      n_rows = if (is_dt) nrow(rows) else length(rows)
      if (n_rows == 0L && skipped_count == 0L) return(invisible(NULL))
      if (n_rows > 0L) {
        # Deduplicate: workers run in fresh processes with an empty C++ key store,
        # so they re-execute combos that the main process already has on record.
        # Drop those rows here before writing to avoid corrupt data and inflated
        # progress_count (which would trigger a premature "already complete" exit).
        if (is_dt) {
          already_done = check_in_result_key_store_cpp(
            rows$response_type, rows$cond_exp_func_model, rows$n, rows$p, rows$betaT,
            rows$rep, rows$design, rows$inference, rows$inference_type
          )
          if (any(already_done)) {
            rows = rows[!already_done]
            n_rows = nrow(rows)
          }
        }
      }
      if (n_rows > 0L) {
        # 1. Update memory store
        # Optimization: Store data.tables directly in raw_results list
        private$results_idx = private$results_idx + 1L
        private$raw_results[[private$results_idx]] = rows
        # 2. Update seen keys
        if (is_dt) {
          add_to_result_key_store_cpp(
            rows$response_type, rows$cond_exp_func_model, rows$n, rows$p, rows$betaT,
            rows$rep, rows$design, rows$inference, rows$inference_type
          )
        } else {
          # rows is a list of lists - convert to vectors for C++
          add_to_result_key_store_cpp(
            vapply(rows, `[[`, "", "response_type"),
            vapply(rows, `[[`, "", "cond_exp_func_model"),
            vapply(rows, `[[`, 0L, "n"),
            vapply(rows, `[[`, 0L, "p"),
            vapply(rows, `[[`, 0, "betaT"),
            vapply(rows, `[[`, 0L, "rep"),
            vapply(rows, `[[`, "", "design"),
            vapply(rows, `[[`, "", "inference"),
            vapply(rows, `[[`, "", "inference_type")
          )
        }
        # 3. Buffer for periodic disk flush
        private$pending_file_rows[[length(private$pending_file_rows) + 1L]] = rows
      }
      # 4. Advance progress count
      # ONLY add n_rows (the new ones).
      # DO NOT add skipped_count because those tasks were already accounted for
      # in the initial progress_count calculation (the "scanning" phase).
      private$progress_count = private$progress_count + n_rows
      private$.draw_progress()
      invisible(NULL)
    },
    # ── Defaults ──────────────────────────────────────────────────────────────
    .default_design_classes = function() {
      list(
        # ── Fixed ──────────────────────────────────────────────────────────────
        DesignFixedBernoulli,
        DesignFixediBCRD,
        DesignFixedBinaryMatch,
        DesignFixedBlocking,                # strata_cols auto-injected if absent
        DesignFixedGreedy,
        DesignFixedMatchingGreedyPairSwitching,
        DesignFixedRerandomization,
        DesignFixedOptimalBlocks,
        DesignFixedCluster,                 # cluster_col auto-injected if absent
        DesignFixedBlockedCluster,          # strata_cols + cluster_col auto-injected if absent
        DesignFixedDOptimal,
        DesignFixedAOptimal,
        DesignFixedFactorial,               # factors auto-injected if absent
        # ── Sequential one-by-one ──────────────────────────────────────────────
        DesignSeqOneByOneBernoulli,
        DesignSeqOneByOneiBCRD,
        DesignSeqOneByOneEfron,
        DesignSeqOneByOneAtkinson,
        DesignSeqOneByOneUrn,
        DesignSeqOneByOneRandomBlockSize,   # strata_cols auto-injected if absent
        DesignSeqOneByOneSPBR,              # strata_cols auto-injected if absent
        DesignSeqOneByOnePocockSimon,       # strata_cols auto-injected if absent
        DesignSeqOneByOneKK21,
        DesignSeqOneByOneKK21stepwise,
        DesignSeqOneByOneKK14
      )
    },
    .default_inference_classes = function() {
      rt   = private$response_type_values[[1L]]
      univ = if (!(rt == "survival" && private$prob_censoring > 0)) list(InferenceAllSimpleMeanDiff) else list()
      type_specific = switch(rt,
        continuous = list(
          InferenceAllSimpleWilcox,
          InferenceContinOLS,
          InferenceContinLin,
          InferenceContinRobustRegr,
          InferenceContinRobustRegr,
          InferenceContinKKOLSIVWC,
          InferenceContinKKOLSOneLik
        ),
        incidence = list(
          InferenceIncidLogRegr,
          InferenceIncidLogRegr,
          InferenceIncidModifiedPoisson,
          InferenceIncidModifiedPoisson,
          InferenceIncidKKCondLogitIVWC,
          InferenceIncidKKCondLogitOneLik,
          InferenceIncidCMH,
          InferenceIncidExtendedRobins,
          InferenceIncidenceExactZhang,
          InferenceIncidExactFisher,
          InferenceIncidExactBinomial
        ),
        proportion = list(
          InferenceAllSimpleWilcox,
          InferencePropBetaRegr,
          InferencePropBetaRegr,
          InferencePropFractionalLogit,
          InferencePropFractionalLogit,
          InferencePropKKGEE,
          InferencePropKKQuantileRegrIVWC
        ),
        count = list(
          InferenceAllSimpleWilcox,
          InferenceCountPoisson,
          InferenceCountPoisson,
          InferenceCountRobustPoisson,
          InferenceCountRobustPoisson,
          InferenceCountKKGEE
        ),
        survival = list(
          InferenceSurvivalCoxPHRegr,
          InferenceSurvivalCoxPHRegr,
          InferenceSurvivalLogRank,
          InferenceSurvivalRestrictedMeanDiff,
          InferenceSurvivalKKStratCoxPHIVWC,
          InferenceSurvivalKKLWACoxPHIVWC
        ),
        ordinal = list(
          InferenceOrdinalPropOddsRegr,
          InferenceOrdinalOrderedProbitRegr,
          InferenceOrdinalCloglogRegr,
          InferenceOrdinalKKGEE
        ),
        stop("Unknown response_type: ", rt)
      )
      c(univ, type_specific)
    }
  )
)
