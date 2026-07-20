#' Jonckheere-Terpstra (JT) Test for Ordinal Responses
#'
#' Exact Jonckheere-Terpstra (JT) rank test for a two-arm ordered alternative with an
#' ordinal response. For treatment versus control, the test statistic is the
#' sum of Mann-Whitney U counts across groups. This class provides the exact
#' distribution-based p-value.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "ordinal",
#'   verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalJonckheereTerpstraTest$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalJonckheereTerpstraTest = R6::R6Class(
	"InferenceOrdinalJonckheereTerpstraTest",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description Initialize the JT test object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Returns the estimated treatment effect (JT superiority measure).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$compute_asymptotic_jt_components(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Returns the weighted JT superiority estimate for Bayesian-bootstrap re-estimation.
		#' @param subject_or_block_weights Bootstrap weights at the subject/block level.
		#' @param estimate_only If TRUE, skip exact p-value calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			private$cached_values$beta_hat_T = private$weighted_superiority(private$y, private$w, row_weights) - 0.5
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Returns the exact two-sided p-value.
		compute_exact_two_sided_pval_for_treatment_effect = function(){
			private$compute_exact_jt_components()
			private$cached_values$p_exact
		},
		#' @description Computes the asymptotic normal confidence interval using the same
		#' null-variance approximation as \code{clinfun::jonckheere.test()}.
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_asymptotic_jt_components()
			if (!is.finite(private$cached_values$s_beta_hat_T)) return(c(NA_real_, NA_real_))
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the asymptotic normal two-sided p-value using the same
		#' \eqn{Z}-approximation as \code{clinfun::jonckheere.test()}.
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$compute_asymptotic_jt_components()
			if (!is.finite(private$cached_values$s_beta_hat_T)) return(NA_real_)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# ordinal: no sharp-null shift supported
			if (delta != 0) return(NULL)
			# "smoothed" adds continuous Gaussian noise, which is not meaningful for integer
			# ordinal category codes; decline and let the R-level fallback (which truncates
			# via as.integer()) handle it, unchanged from before this kernel existed.
			if (length(rand_bootstrap_draws) > 0L && !is.null(rand_bootstrap_draws[[1L]][["smooth_noise"]])) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_jt_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), mats$i_mat, mats$w_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
		weighted_superiority = function(y_vals, w_vals, row_weights){
			y_vals = as.numeric(y_vals)
			w_vals = as.integer(w_vals)
			row_weights = as.numeric(row_weights)
			i_t = which(w_vals == 1L & is.finite(y_vals) & is.finite(row_weights) & row_weights > 0)
			i_c = which(w_vals == 0L & is.finite(y_vals) & is.finite(row_weights) & row_weights > 0)
			if (length(i_t) == 0L || length(i_c) == 0L) return(NA_real_)
			y_t = y_vals[i_t]
			y_c = y_vals[i_c]
			w_t = row_weights[i_t]
			w_c = row_weights[i_c]
			diffs = outer(y_t, y_c, "-")
			comp = (diffs > 0) + 0.5 * (diffs == 0)
			w_pair = outer(w_t, w_c, "*")
			den = sum(w_pair)
			if (!is.finite(den) || den <= 0) return(NA_real_)
			sum(comp * w_pair) / den
		},
		compute_asymptotic_jt_components = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			y = as.integer(private$y)
			w = as.integer(private$w)
			ok = is.finite(y) & is.finite(w)
			y = y[ok]
			w = w[ok]
			n_treat = sum(w == 1L)
			n_control = sum(w == 0L)
			if (n_treat == 0L || n_control == 0L) {
				private$cache_nonestimable_estimate("jt_empty_treatment_arm")
				return(invisible(NULL))
			}
			levs = sort(unique(y))
			treat_counts = tabulate(match(y[w == 1L], levs), nbins = length(levs))
			control_counts = tabulate(match(y[w == 0L], levs), nbins = length(levs))
			control_below = cumsum(c(0, head(control_counts, -1L)))
			u_stat = sum(treat_counts * (control_below + 0.5 * control_counts))
			superiority = u_stat / (n_treat * n_control)
			private$cached_values$superiority = superiority
			private$cached_values$beta_hat_T = superiority - 0.5
			private$cached_values$jt_u_stat = u_stat
			private$cached_values$jt_n_treat = n_treat
			private$cached_values$jt_n_control = n_control
			if (estimate_only) return(invisible(NULL))
			jtvar = n_treat * n_control * (n_treat + n_control + 1) / 12
			se = sqrt(jtvar) / (n_treat * n_control)
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$df = NA_real_
			invisible(NULL)
		},
		compute_exact_jt_components = function(){
			if (!is.null(private$cached_values$p_exact)) return(invisible(NULL))
			res = exact_jonckheere_terpstra_pval_cpp(as.integer(private$y), as.integer(private$w))
			if (is.null(private$cached_values$s_beta_hat_T)) {
				private$compute_asymptotic_jt_components()
			}
			private$cached_values$superiority = res$superiority
			private$cached_values$beta_hat_T = res$superiority - 0.5
			private$cached_values$p_exact = res$p_exact
			private$cached_values$p_lower = res$p_lower
			private$cached_values$p_upper = res$p_upper
			private$cached_values$jt_stat2 = res$stat2
			invisible(NULL)
		}
	)
)
