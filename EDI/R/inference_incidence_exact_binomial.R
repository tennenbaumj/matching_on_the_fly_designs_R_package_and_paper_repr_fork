#' Exact Binomial Incidence Inference for Matched-Pair Designs
#'
#' Performs exact matched-pair binomial inference for binary outcomes using only
#' discordant matched pairs. This class is available for
#' \code{DesignFixedBinaryMatch} and KK matching-on-the-fly designs. For KK
#' designs, only the matched-pair data are used and the reservoir is ignored.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidExactBinomial$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidExactBinomial = R6::R6Class("InferenceIncidExactBinomial",
	lock_objects = FALSE,
	inherit = InferenceExact,
	public = list(
		#' @description Initialize exact matched-pair binomial inference for incidence outcomes.
		#' @param des_obj A completed design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		#' @return A new \code{InferenceIncidExactBinomial} object.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (!private$design_supports_exact_binomial()) {
				stop("Exact binomial incidence inference requires DesignFixedBinaryMatch or KK matching designs.")
			}
			if (inherits(des_obj, "DesignFixedBinaryMatch")) {
				private$des_obj_priv_int$ensure_matching_structure_computed()
			}
		},
		#' @description Compute the matched-pair treatment estimate on the log-odds scale.
		#' @param estimate_only Ignored for this estimator.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			private$get_exact_binomial_log_or_estimate()
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		}
	),
	private = list(
		default_exact_type = "Binomial",
		resolve_exact_type = function(type){
			if (is.null(type)) type = private$default_exact_type
			if (should_run_asserts()) {
				assertChoice(type, c("Binomial"))
			}
			type
		},
		normalize_exact_inference_args = function(type, args_for_type = NULL){
			if (should_run_asserts()) {
				assertChoice(type, c("Binomial"))
				assertList(args_for_type, null.ok = TRUE)
			}
			utils::modifyList(setNames(list(list()), type), if (is.null(args_for_type)) list() else args_for_type)
		},
		assert_exact_inference_params = function(type, args_for_type){
			if (should_run_asserts()) {
				assertChoice(type, c("Binomial"))
				assertList(args_for_type)
				if (!(type %in% names(args_for_type))) stop("args_for_type must contain a list for ", type)
			}
			args = args_for_type[[type]]
			if (should_run_asserts()) {
				assertList(args)
				assertResponseType(private$des_obj$get_response_type(), "incidence")
				assertNoCensoring(private$any_censoring)
			}
			if (!private$design_supports_exact_binomial()) {
				stop("Exact binomial incidence inference requires DesignFixedBinaryMatch or KK matching designs.")
			}
			stats = private$get_exact_binomial_stats()
			if (should_run_asserts()) {
				if (stats$m <= 0L) {
					stop("Exact binomial incidence inference requires at least one matched pair.")
				}
				if (stats$d_plus + stats$d_minus <= 0L) {
					stop("Exact binomial incidence inference requires at least one discordant matched pair.")
				}
			}
			invisible(args)
		},
		compute_exact_confidence_interval_by_type = function(type, alpha, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Binomial = private$ci_exact_binomial(alpha)
			)
		},
		compute_exact_two_sided_pval_for_treatment_effect_by_type = function(type, delta, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Binomial = private$pval_exact_binomial(delta)
			)
		},
		design_supports_exact_binomial = function(){
			is(private$des_obj, "DesignFixedBinaryMatch") || is(private$des_obj, "DesignSeqOneByOneKK14")
		},
		pval_exact_binomial = function(delta_0){
			stats = private$get_exact_binomial_stats()
			if (stats$m <= 0L) {
				private$cache_nonestimable_estimate("exact_binomial_no_matched_pairs")
				return(NA_real_)
			}
			if (stats$d_plus + stats$d_minus <= 0L) {
				return(1)
			}
			zhang_exact_binom_pval_cpp(stats$d_plus, stats$d_minus, delta_0)
		},
		ci_exact_binomial = function(alpha){
			stats = private$get_exact_binomial_stats()
			d_total = stats$d_plus + stats$d_minus
			if (d_total <= 0L) {
				private$cache_nonestimable_estimate("exact_binomial_no_discordant_pairs")
				ci = c(NA_real_, NA_real_)
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}
			ci_prob = stats::binom.test(stats$d_plus, d_total, conf.level = 1 - alpha)$conf.int
			ci = stats::qlogis(ci_prob)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		get_exact_binomial_log_or_estimate = function(){
			stats = private$get_exact_binomial_stats()
			if (stats$m <= 0L) return(NA_real_)
			log((stats$d_plus + 0.5) / (stats$d_minus + 0.5))
		},
		get_exact_binomial_stats = function(){
			if (!is.null(private$cached_values$incidence_exact_binomial_stats)) {
				return(private$cached_values$incidence_exact_binomial_stats)
			}
			if (is(private$des_obj, "DesignFixedBinaryMatch")) {
				private$des_obj_priv_int$ensure_matching_structure_computed()
			}
			m_vec = private$des_obj_priv_int$m
			if (is.null(m_vec) || length(m_vec) == 0L) {
				stop("Matching structure is unavailable for exact binomial incidence inference.")
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			KKstats = compute_zhang_match_data_cpp(private$get_X(), private$y, private$w, m_vec)
			stats = list(
				m = as.integer(KKstats$m),
				d_plus = as.integer(KKstats$d_plus),
				d_minus = as.integer(KKstats$d_minus)
			)
			private$cached_values$incidence_exact_binomial_stats = stats
			stats
		}
	)
)
