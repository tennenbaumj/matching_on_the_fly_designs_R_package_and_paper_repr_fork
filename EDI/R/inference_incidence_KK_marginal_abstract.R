#' Abstract class for all-subject marginal incidence inference in KK designs
#'
#' @keywords internal
InferenceAbstractKKMarginalIncid = R6::R6Class("InferenceAbstractKKMarginalIncid",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = utils::modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},
		get_cluster_ids = function(){
			des_priv = private$des_obj_priv_int
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec_int = as.integer(m_vec)
			m_vec_int[is.na(m_vec_int)] = 0L
			# Normalize design's m_vec the same way
			des_m = des_priv$m
			if (is.null(des_m)) des_m = rep(NA_integer_, private$n)
			des_m_int = as.integer(des_m)
			des_m_int[is.na(des_m_int)] = 0L
			# Check design-level cache (only when m_vec matches design's m_vec)
			if (!is.null(des_priv$cluster_id) && identical(m_vec_int, des_m_int)){
				return(des_priv$cluster_id)
			}
			# Check inference-level cache (for bootstrap resamples)
			if (!is.null(private$cached_values$cluster_id) &&
				identical(m_vec_int, private$cached_values$cluster_id_m_vec)){
				return(private$cached_values$cluster_id)
			}
			cluster_id = des_priv$compute_matching_cluster_ids(m_vec_int)
			# Store at design level if this is the original m_vec
			if (identical(m_vec_int, des_m_int)){
				des_priv$cluster_id = cluster_id
				des_priv$cluster_id_m_vec = m_vec_int
			} else {
				private$cached_values$cluster_id = cluster_id
				private$cached_values$cluster_id_m_vec = m_vec_int
			}
			cluster_id
		}
	))
)
