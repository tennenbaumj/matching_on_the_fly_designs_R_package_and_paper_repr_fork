#' Internal Resampling Draw Contracts
#'
#' These contracts keep the draw-loader and distribution-cache vocabulary aligned
#' across randomization, nonparametric bootstrap, bootstrap randomization, and
#' Bayesian bootstrap execution paths.
#'
#' @keywords internal
#' @noRd
EDI_RESAMPLING_DRAW_CONTRACTS = list(
	rand = list(
		operation = "rand",
		draw_type = "assignment",
		loader = "load_randomization_draw_into_worker",
		estimator = "compute_randomization_worker_estimate",
		cache_name = "rand_distr_cache",
		cache_key_method = "build_randomization_distribution_cache_key"
	),
	non_param_boot = list(
		operation = "non_param_boot",
		draw_type = "row_sample",
		loader = "load_non_param_bootstrap_draw_into_worker",
		estimator = "compute_bootstrap_worker_estimate",
		cache_name = "boot_distr_cache",
		cache_key_method = NULL
	),
	rand_bootstrap = list(
		operation = "rand_bootstrap",
		draw_type = "row_sample_plus_assignment",
		loader = "load_rand_bootstrap_draw_into_worker",
		estimator = "compute_bootstrap_worker_estimate",
		cache_name = "rand_boot_distr_cache",
		cache_key_method = NULL
	),
	bayesian_boot = list(
		operation = "bayesian_boot",
		draw_type = "weights_plus_context",
		loader = "load_bayesian_bootstrap_draw_into_worker",
		estimator = "compute_bayesian_bootstrap_worker_estimate",
		cache_name = "bayes_boot_distr_cache",
		cache_key_method = "bayesian_bootstrap_cache_key"
	)
)

#' Returns a resampling draw-loader/cache contract.
#'
#' @keywords internal
#' @noRd
resampling_draw_contract = function(operation){
	if (!is.character(operation) || length(operation) != 1L || is.na(operation)) {
		stop("operation must be one resampling operation name.", call. = FALSE)
	}
	contract = EDI_RESAMPLING_DRAW_CONTRACTS[[operation]]
	if (is.null(contract)) {
		stop("Unknown resampling operation: ", operation, call. = FALSE)
	}
	contract
}

#' Reads a resampling distribution from the operation-specific cache.
#'
#' @keywords internal
#' @noRd
resampling_distribution_cache_get = function(cached_values, operation, cache_key){
	cache_name = resampling_draw_contract(operation)$cache_name
	cache = cached_values[[cache_name]]
	if (is.null(cache)) return(NULL)
	cache[[cache_key]]
}

#' Writes a resampling distribution to the operation-specific cache.
#'
#' @keywords internal
#' @noRd
resampling_distribution_cache_set = function(cached_values, operation, cache_key, value){
	cache_name = resampling_draw_contract(operation)$cache_name
	if (is.null(cached_values[[cache_name]])) cached_values[[cache_name]] = list()
	cached_values[[cache_name]][[cache_key]] = value
	cached_values
}

#' Ensures an operation-specific resampling distribution cache exists.
#'
#' @keywords internal
#' @noRd
resampling_distribution_cache_ensure = function(cached_values, operation){
	cache_name = resampling_draw_contract(operation)$cache_name
	if (is.null(cached_values[[cache_name]])) cached_values[[cache_name]] = list()
	cached_values
}
