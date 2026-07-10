#' Inference for A Sequential Design
#'
#' An abstract R6 Class that estimates, tests and provides intervals for a
#' treatment effect in a completed design.
#' This class takes a completed \code{Design} object as an input where this object
#' contains data for a fully completed experiment (i.e. all treatment
#' assignments were allocated and all responses were collected).
#'
#' @keywords internal
Inference = R6::R6Class("Inference",
	lock_objects = FALSE,
	public = list(
		#' @description Initialize an estimation and test object after the design is completed.
		#' @param des_obj         A completed \code{Design} object whose entire n subjects are
		#'   assigned and response y is recorded within.
		#' @param verbose Whether to print progress messages.
		#' @param harden  Whether to apply robustness measures (default \code{TRUE}). When
		#'   \code{TRUE}, the inference methods employ defensive strategies including QR-based
		#'   rank reduction of the design matrix, progressive correlation-threshold dropping,
		#'   and fallback fits (e.g.\ robust survival regression, treatment-only models) to
		#'   avoid crashes on ill-conditioned data. When \code{FALSE}, the vanilla algorithm
		#'   runs on the full design matrix as supplied; any rank deficiency or convergence
		#'   failure will surface as an error rather than being silently worked around. Set to
		#'   \code{FALSE} when you want to verify that the raw model converges without
		#'   intervention.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param smart_cold_start_default Whether to use smart cold start values by default for
		#'   likelihood-based models. Explicit starts always override this object-level policy.
		#'   \code{NULL} (default) consults the global cold-start dispatch policy.
		#' @param seed Integer seed for reproducibility.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_cold_start_default = NULL, seed = NULL){
			if (is.null(smart_cold_start_default)) {
				smart_cold_start_default = edi_cold_start_dispatch_policy(class(self)[1])
			}
			if (should_run_asserts()) {
				assertClass(des_obj, "Design")
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(verbose)
				assertFlag(harden)
				assertFlag(smart_cold_start_default)
				des_obj$assert_all_responses_recorded()
			}
			private$harden = harden
			private$smart_cold_start_default = smart_cold_start_default
			private$cached_values = list()
			private$any_censoring = des_obj$any_censoring()
			private$des_obj = des_obj
			private$des_obj_priv_int = des_obj$.__enclos_env__$private
			private$y = private$des_obj_priv_int$y
			private$y_temp = private$y
			private$w = private$des_obj_priv_int$w
			private$dead = private$des_obj_priv_int$dead
			private$is_KK = inherits(des_obj, "DesignSeqOneByOneKK14")
			private$has_match_structure = private$is_KK || inherits(des_obj, "DesignFixedBinaryMatch")
			private$n = des_obj$get_n()
			private$prob_T = des_obj$get_prob_T()
			private$supports_design_resampling = isTRUE(des_obj$supports_resampling())
			# Handle model_formula and X matrix construction
			if (is.null(model_formula)) {
				private$model_formula = des_obj$get_model_formula()
				# Ensure the design matrix is built
				private$X = private$get_X()
			} else {
				# Ensure design is ready
				private$des_obj_priv_int$covariate_impute_if_necessary_and_then_create_model_matrix()
				
				if (should_run_asserts()) {
					all_features = names(des_obj$get_X_imp())
					required_vars = all.vars(model_formula)
					# '.' is a special symbol in formulas representing all variables
					required_vars = setdiff(required_vars, ".")
					if (length(required_vars) > 0 && !all(required_vars %in% all_features)) {
						stop("model_formula contains variables not present in the design's covariates: ", 
							 paste(setdiff(required_vars, all_features), collapse = ", "))
					}
				}
				private$model_formula = model_formula
				# Path #3: final numeric design matrix
				private$X = create_model_matrix_from_features(private$model_formula, des_obj$get_X_imp())
			}
			
			private$verbose = verbose
			private$seed = seed
			private$cached_values$rand_distr_cache = list()
			private$cached_values$m_cache = list()
			private$cached_values$likelihood_test_eval_cache = list()
			if (private$verbose){
				cat(paste0(
					"Initialized inference methods for a ",
					class(des_obj)[1],
					" design and response type ",
					des_obj$get_response_type(),
					".\n"
				))
			}
		},
		#' @description Computes an exact two-sided p-value. Subclasses that support exact inference override this.
		#' @param ... Other arguments passed to the method.
		compute_exact_two_sided_pval_for_treatment_effect = function(...){
			stop("Exact inference is only supported for exact inference classes.")
		},
		#' @description Computes an exact confidence interval. Subclasses that support exact inference override this.
		#' @param ... Other arguments passed to the method.
		compute_exact_confidence_interval = function(...){
			stop("Exact inference is only supported for exact inference classes.")
		},
		#' @description Computes an asymptotic two-sided p-value. Subclasses that support asymptotic inference override this.
		#' @param delta Null treatment effect.
		compute_asymp_two_sided_pval = function(delta = 0){
			stop("Asymptotic inference is not implemented for this inference class.")
		},
		#' @description Computes an asymptotic confidence interval. Subclasses that support asymptotic inference override this.
		#' @param alpha Significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			stop("Asymptotic inference is not implemented for this inference class.")
		},
		#' @description Computes the treatment estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		#' @return A numeric treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			stop("Must be implemented by concrete class.")
		},
		#' @description Returns whether the most recent inference attempt explicitly marked the
		#' result as non-estimable.
		#' @param type Which stage to query: \code{"any"}, \code{"estimate"}, or \code{"se"}.
		#' @return A logical scalar.
		is_nonestimable = function(type = c("any", "estimate", "se")){
			type = match.arg(type)
			if (!isTRUE(private$cached_values$nonestimable)) return(FALSE)
			stage = private$cached_values$nonestimable_stage
			if (identical(type, "any")) return(TRUE)
			if (identical(type, "estimate")) return(identical(stage, "estimate"))
			if (identical(type, "se")) return(identical(stage, "se"))
			FALSE
		},
		#' @description Returns the reason recorded for the most recent explicit non-estimability.
		#' @return A character scalar or \code{NULL}.
		get_nonestimable_reason = function(){
			private$cached_values$nonestimable_reason
		},
		#' @description Returns the stage recorded for the most recent explicit non-estimability.
		#' @return A character scalar or \code{NULL}.
		get_nonestimable_stage = function(){
			private$cached_values$nonestimable_stage
		},
		#' @description Duplicate this inference object
		#' @param verbose 	A flag indicating whether messages should be displayed.
		#' @param make_fork_cluster 	Whether the duplicate should be allowed to create a fork 
		#'   cluster. Default FALSE.
		#' @return 			A new \code{Inference} object with the same data
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			i = self$clone()
			i$.__enclos_env__$private$verbose = verbose
			i$.__enclos_env__$private$fork_cluster = NULL
			i$.__enclos_env__$private$cached_values = list()
			i$.__enclos_env__$private$cached_values$m_cache = private$cached_values$m_cache
			i$.__enclos_env__$private$cached_values$t0s_rand = private$cached_values$t0s_rand
			i$.__enclos_env__$private$cached_values$likelihood_test_eval_cache = list()
			if (!is.null(private$active_resampling_operation)) {
				ws_val = edi_warm_start_dispatch_policy(class(self)[1], private$active_resampling_operation, n = private$des_obj$get_t())
				i$.__enclos_env__$private$fit_warm_start_enabled = ws_val
				i$.__enclos_env__$private$null_fit_warm_start_enabled = ws_val
				i$.__enclos_env__$private$active_resampling_operation = private$active_resampling_operation
			}
			if (private$has_private_method("custom_randomization_statistic_function") &&
				!is.null(i$.__enclos_env__$private$custom_randomization_statistic_function)){
				clone_private = i$.__enclos_env__$private
				fn = clone_private$custom_randomization_statistic_function
				if (bindingIsLocked("custom_randomization_statistic_function", clone_private)) {
					# Use a trick to avoid CRAN check for unsafe unlockBinding
					get("unlockBinding", envir = asNamespace("base"))("custom_randomization_statistic_function", clone_private)
				}
				clone_private[["custom_randomization_statistic_function"]] = fn
			}
			if (private$has_private_method("compiled_cpp_stat_src") &&
				!is.null(private[["compiled_cpp_stat_src"]])){
				clone_private = i$.__enclos_env__$private
				src = private[["compiled_cpp_stat_src"]]
				if (bindingIsLocked("compiled_cpp_stat_src", clone_private)) {
					get("unlockBinding", envir = asNamespace("base"))("compiled_cpp_stat_src", clone_private)
				}
				clone_private[["compiled_cpp_stat_src"]] = src
			}
			i
		},
		#' @description Return the response vector used by inference extension classes.
		#'
		#' This accessor is part of the supported extension contract for user-defined
		#' R6 inference classes. Prefer this method over direct access to private
		#' fields.
		#'
		#' @return A numeric response vector.
		get_response = function(){
			private$y
		},
		#' @description Return the treatment-assignment vector used by inference extension classes.
		#'
		#' This accessor is part of the supported extension contract for user-defined
		#' R6 inference classes. Treatment is encoded as 0/1.
		#'
		#' @return A numeric or integer 0/1 treatment vector.
		get_treatment = function(){
			private$w
		},
		#' @description Return the processed covariate matrix used by inference extension classes.
		#'
		#' This accessor returns the design object's model-matrix covariates, after
		#' the package's missingness handling and encoding. It may be \code{NULL} if
		#' no covariates are available.
		#'
		#' @return A numeric matrix of covariates or \code{NULL}.
		get_covariates = function(){
			private$des_obj$get_X()
		},
		#' @description Return a data frame with response, treatment, censoring status, and covariates.
		#'
		#' This accessor is the preferred data interface for user-defined R6
		#' inference classes. It avoids reliance on private implementation fields.
		#' The returned data frame always contains \code{y}, \code{w}, and
		#' \code{dead}; covariate columns are appended when available.
		#'
		#' @return A data frame suitable for user-defined model fitting.
		get_analysis_data = function(){
			out = data.frame(
				y = private$y,
				w = private$w,
				dead = private$dead
			)
			X = self$get_covariates()
			if (!is.null(X) && length(X) > 0L) {
				X = as.data.frame(X)
				if (is.null(names(X)) || any(names(X) == "")) {
					names(X) = paste0("x", seq_len(ncol(X)))
				}
				out = cbind(out, X)
			}
			out
		},
		#' @description Return the completed design object backing this inference object.
		#'
		#' This accessor is part of the supported extension contract. Extension
		#' classes should use this method instead of \code{private$des_obj}.
		#'
		#' @return The completed \code{Design} object.
		get_design_object = function(){
			private$des_obj
		},
		#' @description Return the response type for the backing design.
		#'
		#' @return A character scalar such as \code{"continuous"},
		#'   \code{"incidence"}, \code{"proportion"}, \code{"count"},
		#'   \code{"survival"}, or \code{"ordinal"}.
		get_response_type = function(){
			private$des_obj$get_response_type()
		},
		#' @description Return the model formula used for covariate adjustment.
		#' @return A formula object or \code{NULL}.
		get_model_formula = function(){
			private$model_formula
		},
		#' @description Set the optimizer used by likelihood-based inference implementations.
		#' @param optimization_alg The optimizer name. Valid values are configured by
		#'   the concrete inference class.
		#' @param allow_irls 		Whether to allow IRLS (Iteratively Reweighted Least Squares)
		#'   as a fallback or primary optimization algorithm.
		#' @param default 			The default optimizer to use if none is specified.
		#' @return Invisibly returns \code{self}.
		set_optimization_alg = function(optimization_alg = NULL, allow_irls = private$optimization_alg_allow_irls, default = private$optimization_alg_default){
			private$optimization_alg_allow_irls = allow_irls
			if (missing(default) || is.null(default)) {
				default = edi_optimization_dispatch_policy(class(self)[1])
			}
			private$optimization_alg_default = default
			new_optimization_alg = .normalize_optimizer_algorithm(
				optimization_alg,
				allow_irls = allow_irls,
				default = default
			)
			if (!identical(private$optimization_alg, new_optimization_alg) && length(private$cached_values) > 0L){
				private$cached_values = list()
			}
			private$optimization_alg = new_optimization_alg
			invisible(self)
		},
		#' @description Return the optimizer used by likelihood-based inference implementations.
		get_optimization_alg = function(){
			private$optimization_alg
		},
		#' @description Set the seed for reproducibility.
		#' @param seed Integer seed for reproducibility.
		set_seed = function(seed) { private$seed = seed }
		),
	active = list(
		#' @field num_cores Current number of cores for this inference object.
		#'   Defaults to the global budget unless overridden on the object.
		num_cores = function(value) {
			if (missing(value)) {
				if (!is.null(private$num_cores_override)) return(private$num_cores_override)
				return(get_num_cores())
			}
			if (should_run_asserts()) {
				checkmate::assertCount(value, positive = TRUE)
			}
			private$num_cores_override = as.integer(value)
			invisible(self)
		}
	),
	private = list(
		finalize = function(){
			# We no longer own the cluster, it is global.
		},
		seed = NULL,
		harden = TRUE,
		des_obj = NULL,		des_obj_priv_int = NULL,
		m = NULL,
		is_KK = NULL,
		has_match_structure = NULL,
		supports_design_resampling = FALSE,
		any_censoring = NULL,
		warned_no_parallel = FALSE,
		fork_cluster = NULL,
		verbose = FALSE,
		n = NULL,
		p = NULL,
		prob_T = NULL,
		y = NULL,
		w = NULL,
		dead = NULL,
		y_temp = NULL,
		X = NULL,
		model_formula = NULL,
		smart_cold_start_default = NULL,
		optimization_alg = NULL,
		optimization_alg_allow_irls = FALSE,
		optimization_alg_default = "lbfgs",
		active_resampling_operation = NULL,
		fit_warm_start_enabled = TRUE,
		null_fit_warm_start_enabled = TRUE,
		reusable_bootstrap_worker_enabled = TRUE,
		fit_warm_start = NULL,
		fit_warm_start_type = NULL,
		fit_warm_start_fisher = NULL,
		fit_warm_start_weights = NULL,
		likelihood_null_warm_cache = NULL,
		reduced_design_keep_cache = NULL,
		fixed_covariate_keep_cache = NULL,
		cached_design_matrix = NULL,
		cached_w_for_design_matrix = NULL,
		cached_harden_for_design_matrix = NULL,
		cached_hardened_X_cov = NULL,
		cached_reduced_X = NULL,
		cached_X_full_for_reduced = NULL,
		cached_keep_for_reduced = NULL,
		cached_j_treat_for_reduced = NULL,
			cached_values = list(),
			num_cores_override = NULL,
		# Returns the number of C++ OpenMP threads to use for a parallel C++ function
		# with n_work_items items of work. Caps threads so that each thread handles
		# at least 10 items; for tiny r/B values this prevents thread-management
		# overhead from dominating over the actual computation.
		n_cpp_threads = function(n_work_items) {
			min(self$num_cores, max(1L, as.integer(n_work_items) %/% 10L))
		},
		parallel_dispatch_policy = function(operation) {
			edi_parallel_dispatch_policy(
				inference_class = class(self)[1],
				response_type = private$des_obj$get_response_type(),
				operation = operation
			)
		},
		effective_parallel_cores = function(operation, requested_cores = self$num_cores) {
			requested_cores = max(1L, as.integer(requested_cores))
			policy = private$parallel_dispatch_policy(operation)
			if (requested_cores > 1L && isTRUE(policy$force_serial)) {
				return(1L)
			}
			requested_cores
		},
		clear_nonestimable_state = function(){
			private$cached_values$nonestimable = FALSE
			private$cached_values$nonestimable_reason = NULL
			private$cached_values$nonestimable_stage = NULL
			invisible(NULL)
		},
		clear_fit_warm_start = function(){
			private$fit_warm_start = NULL
			private$fit_warm_start_type = NULL
			private$fit_warm_start_fisher = NULL
			private$fit_warm_start_weights = NULL
			invisible(NULL)
		},
		set_fit_warm_start = function(start, type = c("beta", "params"), fisher = NULL, weights = NULL, force_pd = TRUE){
			if (!isTRUE(private$fit_warm_start_enabled)) {
				private$clear_fit_warm_start()
				return(invisible(NULL))
			}
			
			# Resampling fits should not replace the primary MLE warm state.
			# Randomization is the intentional exception: it chains sequentially
			# across permutations by updating the worker warm start after each
			# successful null fit.
			if (!is.null(private$active_resampling_operation) &&
					!identical(private$active_resampling_operation, "rand")) {
				return(invisible(NULL))
			}
			
			if (isTRUE(force_pd)) {
				# Direct trusted path: skip all checks
				private$fit_warm_start = start
				private$fit_warm_start_type = if (is.character(type)) type[1L] else "beta"
				private$fit_warm_start_fisher = fisher
				private$fit_warm_start_weights = weights
				return(invisible(NULL))
			}
			
			type = match.arg(type)
			if (is.null(start)) {
				private$clear_fit_warm_start()
				return(invisible(NULL))
			}
			start = as.numeric(start)
			if (!length(start) || any(!is.finite(start))) {
				private$clear_fit_warm_start()
				return(invisible(NULL))
			}
			private$fit_warm_start = start
			private$fit_warm_start_type = type
			
			if (!is.null(fisher) && is.matrix(fisher)) {
				fisher = as.matrix(fisher)
				if (nrow(fisher) == length(start) && ncol(fisher) == length(start) && all(is.finite(fisher))) {
					fisher = (fisher + t(fisher)) / 2
					is_pd = isTRUE(tryCatch({
						chol(fisher)
						TRUE
					}, error = function(e) FALSE))
					private$fit_warm_start_fisher = if (is_pd) fisher else NULL
				} else {
					private$fit_warm_start_fisher = NULL
				}
			} else {
				private$fit_warm_start_fisher = NULL
			}
			
			if (!is.null(weights)) {
				weights = as.numeric(weights)
				if (length(weights) > 0 && all(is.finite(weights))) {
					private$fit_warm_start_weights = weights
				} else {
					private$fit_warm_start_weights = NULL
				}
			} else {
				private$fit_warm_start_weights = NULL
			}
			invisible(NULL)
		},
		get_fit_warm_start = function(type = c("beta", "params")){
			if (!isTRUE(private$fit_warm_start_enabled)) return(NULL)
			type = match.arg(type)
			if (!identical(private$fit_warm_start_type, type)) return(NULL)
			start = private$fit_warm_start
			if (is.null(start) || !length(start) || any(!is.finite(start))) return(NULL)
			start
		},
		get_fit_warm_start_for_length = function(type = c("beta", "params"), expected_length = NULL){
			start = private$get_fit_warm_start(match.arg(type))
			if (is.null(start) || is.null(expected_length)) return(start)
			expected_length = as.integer(expected_length)[1L]
			if (!is.finite(expected_length) || expected_length < 1L) return(NULL)
			if (length(start) != expected_length) return(NULL)
			start
		},
		get_fit_warm_start_fisher = function(expected_dim = NULL){
			if (!isTRUE(private$fit_warm_start_enabled)) return(NULL)
			fisher = private$fit_warm_start_fisher
			if (is.null(fisher)) return(NULL)
			if (!is.null(expected_dim)) {
				expected_dim = as.integer(expected_dim)[1L]
				if (nrow(fisher) != expected_dim || ncol(fisher) != expected_dim) return(NULL)
			}
			fisher
		},
		get_fit_warm_start_weights = function(expected_n = NULL){
			if (!isTRUE(private$fit_warm_start_enabled)) return(NULL)
			w = private$fit_warm_start_weights
			if (is.null(w)) return(NULL)
			if (!is.null(expected_n)) {
				expected_n = as.integer(expected_n)[1L]
				if (length(w) != expected_n) return(NULL)
			}
			w
		},
		clear_likelihood_null_warm_cache = function(){
			private$likelihood_null_warm_cache = list()
			invisible(NULL)
		},
		clear_likelihood_test_eval_cache = function(){
			private$cached_values$likelihood_test_eval_cache = list()
			invisible(NULL)
		},
		get_likelihood_test_eval_cache = function(){
			cache = private$cached_values$likelihood_test_eval_cache
			if (is.null(cache)) {
				cache = list()
				private$cached_values$likelihood_test_eval_cache = cache
			}
			cache
		},
		normalize_likelihood_test_delta = function(delta){
			delta = as.numeric(delta)[1L]
			if (!is.finite(delta)) return(delta)
			if (identical(delta, 0) || abs(delta) < .Machine$double.eps) return(0)
			delta
		},
		likelihood_test_delta_key = function(testing_type, delta){
			delta = private$normalize_likelihood_test_delta(delta)
			paste0(as.character(testing_type)[1L], "::", sprintf("%.17g", delta))
		},
		get_likelihood_test_eval_entry = function(testing_type, delta){
			cache = private$get_likelihood_test_eval_cache()
			cache[[private$likelihood_test_delta_key(testing_type, delta)]]
		},
		set_likelihood_test_eval_entry = function(testing_type, delta, entry){
			cache = private$get_likelihood_test_eval_cache()
			cache[[private$likelihood_test_delta_key(testing_type, delta)]] = entry
			private$cached_values$likelihood_test_eval_cache = cache
			invisible(entry)
		},
		get_likelihood_null_warm_state = function(key){
			if (!isTRUE(private$null_fit_warm_start_enabled)) return(NULL)
			cache = private$likelihood_null_warm_cache
			if (is.null(cache) || is.null(cache[[key]])) return(NULL)
			cache[[key]]
		},
		set_likelihood_null_warm_state = function(key, delta, start){
			if (!isTRUE(private$null_fit_warm_start_enabled)) return(invisible(NULL))
			if (is.null(private$likelihood_null_warm_cache)) private$likelihood_null_warm_cache = list()
			private$likelihood_null_warm_cache[[key]] = list(delta = as.numeric(delta)[1L], start = start)
			invisible(NULL)
		},
		use_reusable_bootstrap_worker = function(){
			isTRUE(private$reusable_bootstrap_worker_enabled) &&
				private$has_private_method("supports_reusable_bootstrap_worker") &&
				isTRUE(tryCatch(private$supports_reusable_bootstrap_worker(), error = function(e) FALSE))
		},
		cache_nonestimable_estimate = function(reason = "not_estimable"){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$nonestimable = TRUE
			private$cached_values$nonestimable_reason = as.character(reason)[1L]
			private$cached_values$nonestimable_stage = "estimate"
			if (is.null(private$cached_values$df)) private$cached_values$df = NA_real_
			invisible(NULL)
		},
		cache_nonestimable_se = function(reason = "standard_error_unavailable"){
			if (is.null(private$cached_values$beta_hat_T)) private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$nonestimable = TRUE
			private$cached_values$nonestimable_reason = as.character(reason)[1L]
			private$cached_values$nonestimable_stage = "se"
			if (is.null(private$cached_values$df)) private$cached_values$df = NA_real_
			invisible(NULL)
		},
			par_lapply = function(X, FUN, n_cores = self$num_cores, budget = 1L, show_progress = FALSE, export_list = NULL){
				if (length(X) == 0L) return(list())
				n_cores = max(1L, min(as.integer(n_cores), length(X)))
				budget = max(1L, as.integer(budget))
				chunk_count = min(length(X), max(1L, 4L * n_cores))
				chunk_size = max(1L, ceiling(length(X) / chunk_count))
				chunks = split(X, ceiling(seq_along(X) / chunk_size))
				# Run a whole chunk under the requested worker budget so we do not
				# pay scheduler/export overhead once per iteration.
				RUN_CHUNK = function(chunk) {
					ns = asNamespace("EDI")
					edi_env = ns$edi_env
					prev_override = edi_env$num_cores_override
					prev_threads = getOption(".edi_last_set_threads")
					if (is.null(prev_threads) || length(prev_threads) != 1L || !is.finite(prev_threads)) {
						prev_threads = 1L
					}
					assign("num_cores_override", budget, envir = edi_env)
					ns$set_package_threads(budget)
					on.exit({
						assign("num_cores_override", prev_override, envir = edi_env)
						ns$set_package_threads(prev_threads)
					}, add = TRUE)
					lapply(chunk, FUN)
				}
				flatten_chunk_results = function(results) {
					if (length(results) == 0L) return(list())
					unlist(results, recursive = FALSE, use.names = FALSE)
				}
				if (n_cores <= 1L) {
					if (!is.null(private$seed)) {
						had = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
						if (had) old = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
						on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
						set.seed(private$seed)
					}
					return(RUN_CHUNK(X))
				}
				global_cl = get_global_fork_cluster()
				global_mirai_cores = get_global_mirai_cores()
				if (!is.null(global_cl)){
					worker_cl = global_cl[seq_len(min(n_cores, length(global_cl)))]
					if (!is.null(private$seed)) parallel::clusterSetRNGStream(worker_cl, private$seed)
					tryCatch({
						if (!is.null(export_list) && length(export_list) > 0L) {
							export_env = list2env(export_list, parent = emptyenv())
							parallel::clusterExport(worker_cl, names(export_list), envir = export_env)
						}
						flatten_chunk_results(parallel::parLapply(worker_cl, chunks, RUN_CHUNK))
					}, error = function(e) {
					# If the persistent cluster has been killed externally (for example by a
					# timeout watchdog), replace it with a fresh one and retry once.
					# This prevents stale socket connections from poisoning subsequent calls.
					msg = conditionMessage(e)
						if (.Platform$OS.type == "unix" &&
							grepl("connection|serialize|unserialize|postNode|sendData|recvData", msg, ignore.case = TRUE)) {
							try(parallel::stopCluster(global_cl), silent = TRUE)
							fresh_cl = make_configured_fork_cluster(n_cores)
							edi_env$global_fork_cluster = fresh_cl
							fresh_worker_cl = fresh_cl[seq_len(min(n_cores, length(fresh_cl)))]
							if (!is.null(private$seed)) parallel::clusterSetRNGStream(fresh_worker_cl, private$seed)
							return(flatten_chunk_results(parallel::parLapply(fresh_worker_cl, chunks, RUN_CHUNK)))
						}
							stop(e)
					})
				} else if (!is.null(global_mirai_cores)){
					requested_mirai_cores = min(n_cores, global_mirai_cores)
					private$ensure_mirai_daemons(requested_mirai_cores)
					on.exit(private$ensure_mirai_daemons(global_mirai_cores), add = TRUE)
					seed_val = private$seed
					tasks = lapply(seq_along(chunks), function(i) {
						chunk = chunks[[i]]
						chunk_seed = if (!is.null(seed_val)) seed_val + i else NULL
						mirai::mirai({
							if (!is.null(chunk_seed)) set.seed(chunk_seed)
							RUN_CHUNK(chunk)
						}, RUN_CHUNK = RUN_CHUNK, chunk = chunk, chunk_seed = chunk_seed)
					})
					flatten_chunk_results(lapply(tasks, function(m) m[]))
				} else if (.Platform$OS.type != "unix"){
					if (!isTRUE(private$warned_no_parallel)){
						message("Parallelism (num_cores > 1) requires the 'mirai' package on non-Unix systems. Install it with install.packages('mirai'). Falling back to serial computation.")
						private$warned_no_parallel = TRUE
					}
					RUN_CHUNK(X)
				} else {
					# Unix with no pre-existing cluster: create one lazily and cache it for
					# subsequent calls, giving the same persistent-cluster performance path.
					lazy_cl = make_configured_fork_cluster(n_cores)
					edi_env$global_fork_cluster = lazy_cl
					lazy_worker_cl = lazy_cl[seq_len(min(n_cores, length(lazy_cl)))]
					if (!is.null(private$seed)) parallel::clusterSetRNGStream(lazy_worker_cl, private$seed)
					flatten_chunk_results(parallel::parLapply(lazy_worker_cl, chunks, RUN_CHUNK))
				}
			},
		ensure_mirai_daemons = function(n){
			s = tryCatch(mirai::status(), error = function(e) list(connections = 0L))
			n_running = if (is.numeric(s$connections) && length(s$connections) == 1L) as.integer(s$connections) else 0L
			if (n_running != n) mirai::daemons(n)
			invisible(NULL)
		},
		stable_signature = function(obj){
			raw_sig = serialize(obj, NULL, xdr = FALSE)
			ints = as.integer(raw_sig)
			if (length(ints) == 0L) return("0:0:0")
			modulus = 2147483647
			h1 = 0
			h2 = 0
			step = max(1L, floor(length(ints) / 64L))
			for (i in seq_along(ints)) {
				val = ints[i]
				h1 = (h1 * 131 + val) %% modulus
				if (i == 1L || i == length(ints) || (i %% step) == 0L) {
					h2 = (h2 * 65599 + val + i) %% modulus
				}
			}
			paste(length(ints), as.integer(h1), as.integer(h2), sep = ":")
		},
		extract_dollar_paths = function(expr){
			paths = list()
			if (is.call(expr)) {
				if (identical(expr[[1]], as.name("$")) && length(expr) == 3L) {
					path = private$resolve_dollar_path(expr)
					if (!is.null(path)) paths = c(paths, list(path))
				}
				for (i in seq_along(expr)[-1]) {
					paths = c(paths, private$extract_dollar_paths(expr[[i]]))
				}
			}
			paths
		},
		resolve_dollar_path = function(expr){
			if (is.symbol(expr)) return(as.character(expr))
			if (is.call(expr) && identical(expr[[1]], as.name("$")) && length(expr) == 3L) {
				parent_path = private$resolve_dollar_path(expr[[2]])
				child_name = if (is.symbol(expr[[3]])) as.character(expr[[3]]) else NULL
				if (is.null(parent_path) || is.null(child_name)) return(NULL)
				return(c(parent_path, child_name))
			}
			NULL
		},
		has_private_method = function(method_name){
			method_name %in% names(private)
		},
			object_has_private_method = function(obj, method_name){
				method_name %in% names(obj$.__enclos_env__$private)
			},
			get_or_create_fork_cluster = function(){
				cl = get_global_fork_cluster()
				if (should_run_asserts()) {
					if (is.null(cl)) {
						stop("No global fork cluster is initialized. Call set_num_cores() first.")
					}
				}
				cl
			},
			assert_design_supports_resampling = function(method_family){
				if (isTRUE(private$supports_design_resampling)) return(invisible(NULL))
				if (should_run_asserts()) {
					stop(method_family, " is not available for plain DesignFixed objects. Use asymptotic inference or a concrete design subclass.")
				}
			},
		create_design_matrix = function(){
			dm = private$cached_design_matrix
			if (!is.null(dm)) return(dm)
			
			Xc = private$cached_hardened_X_cov
			if (is.null(Xc)) {
				Xc_raw = private$get_X()
				if (is.null(Xc_raw) || !ncol(Xc_raw)) {
					Xc = NULL
				} else {
					Xc = as.matrix(Xc_raw)
					if (is.null(colnames(Xc))) {
						colnames(Xc) = paste0("x", seq_len(ncol(Xc)))
					}
					if (isTRUE(private$harden)){
						Xc = drop_highly_correlated_cols(Xc, threshold = 0.999)$M
					}
					private$cached_hardened_X_cov = Xc
				}
			}

			if (is.null(Xc)) {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w, Xc)
			}
			
			if (isTRUE(private$harden) && ncol(X_full) > 2L){
				X_full = drop_linearly_dependent_cols(X_full)$M
			}
			private$cached_design_matrix = X_full
			X_full
		},
		get_X = function(){
			if (!is.null(private$X)) return(private$X)  # transient bootstrap override
			if (is.null(private$des_obj_priv_int$X))
				private$des_obj_priv_int$covariate_impute_if_necessary_and_then_create_model_matrix()
			private$des_obj_priv_int$X
		},
		reduce_treatment_only_design_fast = function(X_full){
			if (is.null(dim(X_full)) || ncol(X_full) != 2L || nrow(X_full) == 0L) return(NULL)
			X_mat = as.matrix(X_full)
			intercept = as.numeric(X_mat[, 1L])
			treatment = as.numeric(X_mat[, 2L])
			if (any(!is.finite(intercept)) || any(!is.finite(treatment))) return(NULL)
			if (any(intercept != intercept[1L])) return(NULL)
			if (any(treatment != treatment[1L])) {
				private$reduced_design_keep_cache = c(1L, 2L)
				return(list(X = X_mat, keep = c(1L, 2L), j_treat = 2L))
			}
			list(X = NULL, keep = 1L, j_treat = NA_integer_)
		},
		try_cached_reduced_design_keep = function(X_full, keep = private$reduced_design_keep_cache){
			if (is.null(keep) || !length(keep)) return(NULL)
			
			# Fast track: if we have a fully cached reduced matrix and X_full is the same
			if (!is.null(private$cached_reduced_X) && 
				identical(X_full, private$cached_X_full_for_reduced) &&
				identical(keep, private$cached_keep_for_reduced)) {
				return(list(X = private$cached_reduced_X, keep = keep, j_treat = private$cached_j_treat_for_reduced))
			}

			keep = sort(unique(as.integer(keep)))
			if (any(!is.finite(keep)) || any(keep < 1L) || any(keep > ncol(X_full)) || !(2L %in% keep)) {
				return(NULL)
			}
			X_try = as.matrix(X_full[, keep, drop = FALSE])
			
			# Avoid QR for simple treatment-only designs
			if (ncol(X_try) == 2L) {
				fast = private$reduce_treatment_only_design_fast(X_try)
				if (!is.null(fast) && !is.null(fast$X)) {
					res = list(X = X_try, keep = keep, j_treat = match(2L, keep))
					private$cached_reduced_X = res$X
					private$cached_X_full_for_reduced = X_full
					private$cached_keep_for_reduced = keep
					private$cached_j_treat_for_reduced = res$j_treat
					return(res)
				}
				return(NULL)
			}
			
			if (qr(X_try)$rank != ncol(X_try)) return(NULL)
			res = list(X = X_try, keep = keep, j_treat = match(2L, keep))
			private$cached_reduced_X = res$X
			private$cached_X_full_for_reduced = X_full
			private$cached_keep_for_reduced = keep
			private$cached_j_treat_for_reduced = res$j_treat
			return(res)
		},
		reduce_design_matrix_preserving_treatment = function(X_full){
			if (!private$harden) {
				X_mat = as.matrix(X_full)
				return(list(X = X_mat, keep = seq_len(ncol(X_mat)), j_treat = 2L))
			}
			fast = private$reduce_treatment_only_design_fast(X_full)
			if (!is.null(fast)) return(fast)
			cached = private$try_cached_reduced_design_keep(X_full)
			if (!is.null(cached)) return(cached)
			reduced = qr_reduce_preserve_cols_cpp(as.matrix(X_full), c(1L, 2L))
			keep = as.integer(reduced$keep)
			if (!(2L %in% keep)) return(list(X = NULL, keep = keep, j_treat = NA_integer_))
			private$reduced_design_keep_cache = keep
			list(
				X = reduced$X_reduced,
				keep = keep,
				j_treat = match(2L, keep)
			)
		},
		reduce_design_matrix_preserving_treatment_fixed_covariates = function(X_full){
			if (!private$harden) {
				X_mat = as.matrix(X_full)
				return(list(X = X_mat, keep = seq_len(ncol(X_mat)), j_treat = 2L))
			}
			fast = private$reduce_treatment_only_design_fast(X_full)
			if (!is.null(fast)) return(fast)
			if (is.null(dim(X_full)) || ncol(X_full) <= 2L) return(private$reduce_design_matrix_preserving_treatment(X_full))
			cached = private$try_cached_reduced_design_keep(X_full)
			if (!is.null(cached)) return(cached)
			other_cols = c(1L, seq.int(3L, ncol(X_full)))
			keep_other = private$fixed_covariate_keep_cache
			if (is.null(keep_other) || !length(keep_other) || any(keep_other < 1L) || any(keep_other > ncol(X_full)) || !(1L %in% keep_other)) {
				X_other = as.matrix(X_full[, other_cols, drop = FALSE])
				reduced_other = qr_reduce_preserve_cols_cpp(X_other, 1L)
				keep_other = other_cols[as.integer(reduced_other$keep)]
				private$fixed_covariate_keep_cache = keep_other
			}
			X_other_reduced = as.matrix(X_full[, keep_other, drop = FALSE])
			treatment = as.numeric(X_full[, 2L])
			if (any(!is.finite(treatment)) || any(!is.finite(X_other_reduced))) {
				return(private$reduce_design_matrix_preserving_treatment(X_full))
			}
			X_trial = cbind(X_other_reduced, treatment)
			if (qr(X_trial)$rank <= ncol(X_other_reduced)) {
				return(private$reduce_design_matrix_preserving_treatment(X_full))
			}
			keep = sort(c(keep_other, 2L))
			private$reduced_design_keep_cache = keep
			list(
				X = as.matrix(X_full[, keep, drop = FALSE]),
				keep = keep,
				j_treat = match(2L, keep)
			)
		},
		reduce_design_matrix_preserving_treatment_matrix = function(X_full){
			private$reduce_design_matrix_preserving_treatment(X_full)$X
		},
		fit_with_hardened_qr_column_dropping = function(X_full, fit_fun, fit_ok, required_cols = 1L){
			X_mat = as.matrix(X_full)
			if (is.null(dim(X_mat))){
				X_mat = matrix(X_mat, ncol = 1L)
			}
			if (!ncol(X_mat)){
				return(list(
					fit = tryCatch(fit_fun(X_mat), error = function(e) NULL),
					X = X_mat,
					keep = integer()
				))
			}
			required_cols = sort(unique(as.integer(required_cols)))
			required_cols = required_cols[
				is.finite(required_cols) &
				required_cols >= 1L &
				required_cols <= ncol(X_mat)
			]
			fit_fun_formals = names(formals(fit_fun))
			fit_fun_accepts_keep = "keep" %in% fit_fun_formals || "..." %in% fit_fun_formals
			attempt_fit = function(keep){
				X_try = X_mat[, keep, drop = FALSE]
				colnames(X_try) = colnames(X_mat)[keep]
				list(
					fit = tryCatch({
						if (fit_fun_accepts_keep) fit_fun(X_try, keep) else fit_fun(X_try)
					}, error = function(e) NULL),
					X = X_try,
					keep = keep
				)
			}
			if (!private$harden || ncol(X_mat) <= length(required_cols)){
				return(attempt_fit(seq_len(ncol(X_mat))))
			}
			qr_X = qr(X_mat)
			keep = sort(unique(c(required_cols, qr_X$pivot[seq_len(qr_X$rank)])))
			if (!length(keep)) keep = seq_len(ncol(X_mat))
			removable = rev(setdiff(qr_X$pivot[qr_X$pivot %in% keep], required_cols))
			best_attempt = attempt_fit(keep)
			if (isTRUE(fit_ok(best_attempt$fit, best_attempt$X, best_attempt$keep))){
				return(best_attempt)
			}
			if (!length(removable)){
				return(best_attempt)
			}
			for (k in seq_along(removable)){
				keep_try = sort(setdiff(keep, removable[seq_len(k)]))
				if (!all(required_cols %in% keep_try)) next
				attempt = attempt_fit(keep_try)
				best_attempt = attempt
				if (isTRUE(fit_ok(attempt$fit, attempt$X, attempt$keep))){
					return(attempt)
				}
			}
			best_attempt
		}
	)
)
