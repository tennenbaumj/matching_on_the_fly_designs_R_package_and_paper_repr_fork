#' Internal Mixin Host Contracts
#'
#' Pattern-1 mixins are lists spliced into an R6 class's \code{public} and
#' \code{private} lists. \code{EDI_MIXIN_CONTRACTS} documents the host-private
#' methods and state which a mixin requires but does not define itself.
#' Empty vectors mean that the mixin is self-contained (or is an intentionally
#' empty future extension point). These contracts are deliberately narrow: a
#' method supplied by the mixin itself is not repeated as a host requirement.
#'
#' \code{EDI_MIXIN_COMPOSITIONS} lists every base class that combines two or
#' more mixins. \code{EDI_MIXIN_ALLOWED_COLLISIONS} records the sole deliberate
#' overwrite: the compound KK mixin replaces the pass-through implementation of
#' \code{compute_basic_match_data()}. Tests use these objects to guard against
#' silent method-name overwrites as new mixins are added.
#'
#' @keywords internal
#' @noRd
EDI_MIXIN_CONTRACTS = list(
	InferenceMixinBartlettApprox = list(
		file = "inference_mixin_bartlett_approx.R",
		private_methods = c(
			"supports_lik_ratio_param_bootstrap", "run_param_bootstrap_replicates",
			"param_bootstrap_lr_extreme"
		),
		private_state = "active_resampling_operation"
	),
	InferenceMixinBcaBootstrapCI = list(
		file = "inference_mixin_bca_bootstrap_ci.R",
		private_methods = character(),
		private_state = character()
	),
	InferenceMixinCIInversion = list(
		file = "inference_mixin_ci_inversion.R",
		private_methods = c(
			"cache_nonestimable_se", "compute_wald_confidence_interval_impl",
			"get_likelihood_test_spec", "get_memoized_likelihood_test_eval",
			"get_memoized_likelihood_test_pval", "get_standard_error",
			"normalize_testing_type"
		),
		private_state = c("likelihood_ci_max_abs", "testing_type")
	),
	InferenceMixinCordeiroFerrariApprox = list(
		file = "inference_mixin_cordeiro_ferrari_approx.R",
		private_methods = character(),
		private_state = character()
	),
	InferenceMixinCustomRandomizationStatistic = list(
		file = "inference_mixin_custom_randomization_statistic.R",
		private_methods = character(),
		private_state = "cached_values"
	),
	InferenceMixinInformationMatrix = list(
		file = "inference_mixin_information_matrix.R",
		private_methods = c("get_default_information_source", "get_likelihood_test_spec"),
		private_state = c("cached_mod", "information_preference", "information_source_used")
	),
	InferenceMixinKKGEEShared = list(
		file = "inference_mixin_kk_gee_shared.R",
		private_methods = c(
			"cache_nonestimable_estimate", "cache_nonestimable_se", "clear_nonestimable_state",
			"compute_z_or_t_ci_from_s_and_df", "compute_z_or_t_two_sided_pval_from_s_and_df",
			"create_design_matrix", "expand_subject_or_block_weights_to_row_weights",
			"fit_with_hardened_qr_column_dropping", "gee_family", "gee_response_type",
			"get_fit_warm_start_fisher", "get_fit_warm_start_for_length",
			"get_fit_warm_start_weights", "set_fit_warm_start", "shared_gee_dispatch"
		),
		private_state = c("any_censoring", "cached_values", "harden", "n", "w", "y")
	),
	InferenceMixinKKGLMMShared = list(
		file = "inference_mixin_kk_glmm_shared.R",
		private_methods = c(
			"cache_nonestimable_estimate", "compute_standard_error_from_information_matrix",
			"compute_z_or_t_ci_from_s_and_df", "compute_z_or_t_two_sided_pval_from_s_and_df",
			"create_design_matrix", "expand_subject_or_block_weights_to_row_weights",
			"fit_with_hardened_qr_column_dropping", "glmm_family", "glmm_response_type",
			"shared"
		),
		private_state = c("any_censoring", "cached_values", "harden", "n", "w", "y")
	),
	InferenceMixinKKPassThrough = list(
		file = "inference_mixin_kk_passthrough.R",
		private_methods = c(
			"assert_valid_bootstrap_type", "cache_nonestimable_estimate",
			"effective_parallel_cores",
			"expand_subject_or_block_weights_to_row_weights", "get_X", "has_private_method",
			"object_has_private_method", "par_lapply"
		),
		private_state = c("cached_values", "des_obj_priv_int", "has_match_structure", "n")
	),
	InferenceMixinKKPassThroughCompound = list(
		file = "inference_mixin_kk_passthrough_compound.R",
		private_methods = c("cache_nonestimable_estimate", "compute_basic_kk_match_data_impl"),
		private_state = c("cached_values", "has_match_structure")
	),
	InferenceMixinLemonteGradientApprox = list(
		file = "inference_mixin_lemonte_gradient_approx.R",
		private_methods = character(),
		private_state = character()
	),
	InferenceMixinLikelihoodTestMemoization = list(
		file = "inference_mixin_likelihood_test_memoization.R",
		private_methods = c(
			"cache_nonestimable_estimate", "cache_nonestimable_se",
			"get_bartlett_factor_approx", "get_bartlett_factor_exact",
			"get_likelihood_null_warm_state", "get_likelihood_test_eval_entry",
			"get_likelihood_test_spec", "get_score_test_information_matrix",
			"normalize_likelihood_test_delta", "set_likelihood_null_warm_state",
			"set_likelihood_test_eval_entry", "supports_bartlett_likelihood_ratio_approx",
			"supports_bartlett_likelihood_ratio_exact"
		),
		private_state = "null_fit_warm_start_enabled"
	),
	InferenceMixinQuantileRandCI = list(
		file = "inference_mixin_quantile_rand_ci.R",
		private_methods = c("compute_rand_pval_matched_pairs", "compute_rand_pval_reservoir"),
		private_state = "cached_values"
	),
	InferenceMixinSequentialMCPval = list(
		file = "inference_mixin_sequential_mc_pval.R",
		private_methods = "compute_two_sided_randomization_pval_band",
		private_state = character()
	)
)

EDI_MIXIN_COMPOSITIONS = list(
	InferenceAsympLik = c(
		"InferenceMixinCIInversion", "InferenceMixinInformationMatrix",
		"InferenceMixinLikelihoodTestMemoization"
	),
	InferenceKKPassThroughCompound = c(
		"InferenceMixinKKPassThrough", "InferenceMixinKKPassThroughCompound"
	),
	InferenceNonParamBootstrap = "InferenceMixinBcaBootstrapCI",
	InferenceParamBootstrap = "InferenceMixinBartlettApprox",
	InferenceRand = c(
		"InferenceMixinCustomRandomizationStatistic", "InferenceMixinSequentialMCPval"
	),
	InferenceRandQuantileCI = "InferenceMixinQuantileRandCI"
)

EDI_MIXIN_ALLOWED_COLLISIONS = list(
	InferenceKKPassThroughCompound = list(
		private = "compute_basic_match_data",
		public = character()
	)
)
