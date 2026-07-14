#' An Abstract Experimental Design
#'
#' @name Design
#' @description Internal method.
#' An abstract R6 Class encapsulating the data and functionality for an experimental design.
#' This class takes care of data storage and response handling.
#'
#' @details
#' Throughout the package, treatment assignment vectors \eqn{w} use the
#' \eqn{\{-1, +1\}} encoding: \eqn{+1} indicates a treated subject and \eqn{-1}
#' a control subject.  All public methods that return or accept \eqn{w}
#' (e.g. \code{get_w()}, \code{draw_ws_according_to_design()}) use this
#' convention.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' seq_des = Design$new(n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
Design = R6::R6Class("Design",
	lock_objects = FALSE,
	public = list(
		#' @description Initialize an experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param prob_T    Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature    Flag for missingness indicators.
		#' @param n            The sample size (if fixed).
		#' @param verbose    Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates when building the
		#'   model matrix for inference. One of:
		#'   \describe{
		#'     \item{\code{"impute"} (default)}{Missing values are filled in using random-forest
		#'       imputation (\code{missRanger}, falling back to \code{missForest} on failure).
		#'       The response vector is included as an auxiliary predictor when available.
		#'       This preserves all covariates and all subjects but introduces imputed values
		#'       that influence inference.}
		#'     \item{\code{"drop_column"}}{Any covariate column that contains at least one
		#'       missing value is dropped entirely from the model matrix before inference.
		#'       No values are invented; the remaining complete columns are used as-is.
		#'       This is conservative but transparent.}
		#'     \item{\code{"error"}}{An error is thrown as soon as any missing value is
		#'       detected in the covariate matrix. Use this when you want to guarantee that
		#'       inference runs on exactly the data you supplied, with no silent modification.}
		#'   }
		#' @param model_formula A formula object used to create the design matrix from
		#'   covariates. Default is \code{~ .}.
		#' @param ordinal_levels If the response type is "ordinal", the labels for the levels.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `Design` object
		initialize = function(
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = FALSE,
				n = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				model_formula = ~ .,
				ordinal_levels = NULL,
				seed = NULL
			) {
			if (should_run_asserts()) {
				assertChoice(response_type, c("continuous", "incidence", "proportion", "count", "survival", "ordinal"))
				assertNumeric(prob_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
				assertFlag(include_is_missing_as_a_new_feature)
				assertFlag(verbose)
				assertCount(n, null.ok = TRUE)
				assertCount(seed, null.ok = TRUE)
				assertChoice(missingness_method, c("impute", "drop_column", "error"))
				assertFormula(model_formula)
				if (response_type == "ordinal" && !is.null(ordinal_levels)) {
					assertCharacter(ordinal_levels, min.len = 2, any.missing = FALSE)
				}
			}
			if (is.null(n)){
				private$fixed_sample = FALSE
			} else {
				n = as.integer(n)
				private$n = n
				private$fixed_sample = TRUE
			}
			private$prob_T = prob_T
			private$response_type = response_type
			private$response_type_original = response_type
			private$ordinal_levels = ordinal_levels
			private$original_ordinal_levels = ordinal_levels
			private$include_is_missing_as_a_new_feature = include_is_missing_as_a_new_feature
			private$missingness_method = missingness_method
			private$model_formula = model_formula
			# Ensure budget is respected among openmp and other packages
			private$verbose = verbose
			if (private$fixed_sample){
				private$y = 	     rep(NA_real_, n)
				private$y_original = rep(NA_real_, n)
				private$w = 	     rep(NA_real_, n)
				private$dead =       rep(NA_real_, n)
			}
			if (private$verbose){
				cat(paste0("Initialized a ",
				class(self)[1],
				" experiment with response type ",
				response_type,
				" and ",
				ifelse(private$fixed_sample, "fixed sample", "not fixed sample"),
				 ".\n"))
			}
			private$seed = seed
		},
		#' @description For CARA designs, add a single subject response.
		#'
		#' @param t The subject index.
		#' @param y The response value.
		#' @param dead If the response is censored (0 for survival).
		add_one_subject_response = function(t, y, dead = 1) {
			if (should_run_asserts()) {
				assertNumeric(t, len = 1)
				assertNumeric(y, len = 1)
				assertNumeric(dead, len = 1)
				assertChoice(dead, c(0, 1))
				assertCount(t, positive = TRUE)
				if (t > private$t){
					stop(paste("You cannot add response for subject", t, "when the most recent subjects' record added is", private$t))
				}
			}
			if (length(private$y) >= t & !is.na(private$y[t])){
				warning(paste("Overwriting previous response for t =", t, "y[t] =", private$y[t]))
			}
			if (private$response_type == "ordinal" && is.factor(y)){
				if (should_run_asserts()) {
					assertFactor(y, ordered = TRUE, any.missing = FALSE)
				}
				levs = levels(y)
				private$ordinal_levels = levs
				if (private$response_type_original == "ordinal" && is.null(private$original_ordinal_levels)){
					private$original_ordinal_levels = levs
				}
				y = as.integer(y)
			}
			
			private$assert_y(y, private$response_type)
			if (private$response_type == "survival" && y == 0){
				warning("0 survival responses not allowed --- recording .Machine$double.eps as its value instead")
				y = .Machine$double.eps
			}
			if (should_run_asserts()) {
				if (dead == 0 & private$response_type != "survival"){
					stop("censored observations are only available for survival response types")
				}
			}
			if (private$fixed_sample | t <= length(private$y)){
				private$y[t] = y
				private$y_original[t] = y
				private$dead[t] = dead
			} else if (t == length(private$y) + 1){
				private$y = c(private$y, y)
				private$y_original = c(private$y_original, y)
				private$dead = c(private$dead, dead)
			} else {
				if (should_run_asserts()) {
					stop("You cannot add a response for a subject that has not yet arrived when the sample size is not fixed in advance.")
				}
			}
			private$y_i_t_i[[t]] = private$t
		},
		#' @description For non-CARA designs, add all subject responses.
		#'
		#' @param ys The responses as a numeric vector.
		#' @param deads The binary vector indicating if dead/censored.
		add_all_subject_responses = function(ys, deads = NULL) {
			if (is.null(deads)){
				deads = rep(1, private$t)
			}
			if (should_run_asserts()) {
				if (private$response_type == "ordinal" && is.factor(ys)){
					assertFactor(ys, len = private$t, ordered = TRUE, any.missing = FALSE)
				} else {
					assertNumeric(ys, len = private$t)
				}
				private$assert_y(ys, private$response_type)
				assertNumeric(deads, len = private$t)
				
				if (private$response_type != "survival" && any(deads == 0)){
					stop("censored observations are only available for survival response types")
				}
			}
			
			if (private$response_type == "ordinal" && is.factor(ys)){
				levs = levels(ys)
				private$ordinal_levels = levs
				if (private$response_type_original == "ordinal" && is.null(private$original_ordinal_levels)){
					private$original_ordinal_levels = levs
				}
				ys = as.integer(ys)
			}
			private$y = as.numeric(ys)
			private$y_original = as.numeric(ys)
			private$dead = as.numeric(deads)
			private$y_i_t_i = as.list(seq_len(private$t))
		},
		#' @description For analysis on already-completed experimental data
		#'
		#' @param w A {-1,+1} vector of subject assignments (+1 = treated, -1 = control).
		overwrite_all_subject_assignments = function(w) {
			if (should_run_asserts()) {
				assertIntegerish(w, lower = -1, upper = 1, any.missing = FALSE, len = private$t)
				if (any(!(w %in% c(-1L, 1L)))) {
					stop("overwrite_all_subject_assignments: w must contain only -1 (control) or +1 (treated).")
				}
			}
			private$w = (as.numeric(w) + 1L) / 2L
		},
		#' @description Check if this design was initialized with a fixed sample size n
		#'
		#' @return TRUE if fixed.
		is_fixed_sample_size = function(){
			private$fixed_sample
		},
		#' @description Asserts if all subjects arrived.
		assert_all_subjects_arrived = function(){
			if (should_run_asserts()) {
				if (private$fixed_sample & private$t < private$n){
					stop("This experiment is incomplete as all n subjects haven't arrived yet.")
				}
			}
		},
		#' @description Asserts if all responses are recorded.
		assert_all_responses_recorded = function(){
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
				if (sum(!is.na(private$y)) != length(private$w)){
					stop("This experiment is incomplete as all responses aren't recorded yet.")
				}
			}
		},
		#' @description Checks if the experiment is completed.
		#'
		#' @return  \code{TRUE} if experiment is complete, \code{FALSE} otherwise.
		check_experiment_completed = function(){
			if (private$fixed_sample & private$t < private$n){
				FALSE
			} else if (sum(!is.na(private$y)) != length(private$w)){
				FALSE
			} else {
				TRUE
			}
		},
		#' @description Checks if the experiment has a 50-50 allocation.
		assert_even_allocation = function(){
			if (should_run_asserts()) {
				if (private$prob_T != 0.5){
					stop("This type of design currently only works with even treatment allocation, i.e. you must set prob_T = 0.5 upon initialization")
				}
			}
		},
		#' @description Checks if the experiment has a fixed sample size.
		assert_fixed_sample = function(){
			if (should_run_asserts()) {
				if (!private$fixed_sample){
					stop("This type of design currently only works with fixed sample, i.e., you must specify n upon initialization")
				}
			}
		},
		#' @description Checks if the experiment has any censored responses
		#'
		#' @return  \code{TRUE} if any censored.
		any_censoring = function(){
			sum(private$dead) < length(private$dead)
		},
		#' @description Get t
		#'
		#' @return 			The current number of subjects.
		get_t = function(){
			private$t
		},
		#' @description Get raw X information
		#'
		#' @return 			A data frame of subject data.
		get_X_raw = function(){
			private$Xraw
		},
		#' @description Get imputed X information
		#'
		#' @return 		Same as \code{Xraw} except with imputations.
		get_X_imp = function(){
			private$Ximp
		},
		#' @description Get X matrix
		#'
		#' @return 			A numeric matrix of subject data.
		get_X = function(){
			private$X
		},
		#' @description Get y
		#'
		#' @return 			A numeric vector of subject responses.
		get_y = function(){
			private$y
		},
		#' @description Get y_original
		#'
		#' @return 			A numeric vector of the original subject responses.
		get_y_original = function(){
			private$y_original
		},
		#' @description Get w
		#'
		#' @return 			A {-1,+1} vector of subject assignments (+1 = treated, -1 = control).
		get_w = function(){
			2L * private$w - 1L
		},
		#' @description Draw treatment assignment vectors according to the design.
		#'
		#' @param r Number of vectors to draw. Default is 1.
		#' @return A matrix of size n x r with {-1,+1} entries (+1 = treated, -1 = control).
		draw_ws_according_to_design = function(r = 1L){
			result = private$draw_ws_raw(r)
			2L * result - 1L
		},
		#' @description Get n, the sample size
		#'
		#' @return 			The number of subjects.
		get_n = function(){
			ifelse(private$fixed_sample, private$n, private$t)
		},
		#' @description Get dead
		#'
		#' @return 			A binary vector of whether the subject is dead.
		get_dead = function(){
			private$dead
		},
		#' @description Get probability of treatment
		#'
		#' @return 			The specified probability.
		get_prob_T = function(){
			private$prob_T
		},
		#' @description Checks if the design contains any covariates.
		#'
		#' @return \code{TRUE} if there are covariates, \code{FALSE} otherwise.
		has_covariates = function(){
			!is.null(private$Xraw) && ncol(private$Xraw) > 0L
		},
		#' @description Get response type
		#'
		#' @return 			The specified response type.
		get_response_type = function(){
			private$response_type
		},
		#' @description Get the original response type
		#'
		#' @return 			The original specified response type.
		get_response_type_original = function(){
			private$response_type_original
		},
		#' @description Get ordinal levels
		#'
		#' @return 			The levels of the ordinal response.
		get_ordinal_levels = function(){
			private$ordinal_levels
		},
		#' @description Get original ordinal levels
		#'
		#' @return 			The labels for the levels of the original ordinal response.
		get_original_ordinal_levels = function(){
			private$original_ordinal_levels
		},
		#' @description Get the missingness method
		#'
		#' @return 			The missingness handling method: \code{"impute"}, \code{"drop_column"},
		#'   or \code{"error"}.
		get_missingness_method = function(){
			private$missingness_method
		},
		#' @description Transform the response vector y
		#'
		#' @param transform_fun A function that takes y_original and returns a new y.
		#' @param transformed_response_type The response type of the transformed y.
		#' @param ordinal_levels If the transformed response type is "ordinal", the labels for the levels.
		transform_y = function(transform_fun, transformed_response_type, ordinal_levels = NULL) {
			if (should_run_asserts()) {
				assertFunction(transform_fun)
				if (!identical(names(formals(transform_fun))[1], "y_original")) {
					stop("transform_fun must have its first argument named 'y_original'")
				}
				assertChoice(transformed_response_type, c("continuous", "incidence", "proportion", "count", "survival", "ordinal"))
				if (transformed_response_type == "ordinal" && !is.null(ordinal_levels)) {
					assertCharacter(ordinal_levels, min.len = 2, any.missing = FALSE)
				}
			}
			y_temp = transform_fun(y_original = private$y_original)
			if (should_run_asserts()) {
				assertNumeric(y_temp, len = length(private$y_original))
				private$assert_y(y_temp, transformed_response_type)
			}
			private$y = as.numeric(y_temp)
			private$response_type = transformed_response_type
			private$ordinal_levels = ordinal_levels
			invisible(private$y)
		},
		#' @description Get the model formula
		#'
		#' @return 			The model formula.
		get_model_formula = function(){
			private$model_formula
		},
		#' @description Duplicate this design object
		#'
		#' @param verbose 	A flag for verbosity.
		#' @return 			A new `Design` object with the same data
		duplicate = function(verbose = FALSE){
			if (should_run_asserts()) {
				self$assert_all_responses_recorded() #can't duplicate without the experiment being done
			}
			# Use the built-in R6 clone method (shallow by default) to bypass $new() logic.
			d = self$clone()
			d$.__enclos_env__$private$seed = NULL
			d$.__enclos_env__$private$verbose = verbose
			d
		}
	),
	active = list(
		#' @field num_cores Current number of cores in the global budget.
		num_cores = function() get_num_cores()
	),
	private = list(
		seed = NULL,
		all_subject_data_cache = list(),
		t = 0L,
		n = NULL,
		Xraw = data.table(),
		Ximp = data.table(),
		X = NULL,
		p_raw_t = NULL,
		w = numeric(),
		y = numeric(),
		y_original = numeric(),
		dead = numeric(),
		permutations_cache   = list(),
		lin_centered_covariates = NULL,
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			list(i_b = sample_int_replace_cpp(private$n, private$n), m_vec_b = NULL)
		},
		prob_T = NULL,
		response_type = NULL,
		response_type_original = NULL,
		ordinal_levels = NULL,
		original_ordinal_levels = NULL,
		fixed_sample = NULL,
		include_is_missing_as_a_new_feature = NULL,
		missingness_method = "impute",
		model_formula = NULL,
		verbose = NULL,
		y_i_t_i = list(),	 #at what point during the experiment are the subjects recorded?
		uses_covariates = FALSE, #does this design use the covariates to make assignments? The default is FALSE
		resample_assignment = function(){
			n = private$n
			i_b = sample_int_replace_cpp(n, n)
			private$w    = private$w[i_b]
			private$y    = private$y[i_b]
			private$y_original = private$y_original[i_b]
			private$dead = private$dead[i_b]
			invisible(self)
		},
		assert_y = function(y, response_type) {
			if (should_run_asserts()) {
				if (response_type == "incidence") {
					assertIntegerish(y, lower = 0, upper = 1, any.missing = FALSE)
				} else if (response_type == "proportion") {
					assertNumeric(y, lower = 0, upper = 1, any.missing = FALSE)
				} else if (response_type == "count") {
					assertIntegerish(y, lower = 0, any.missing = FALSE)
				} else if (response_type == "survival") {
					assertNumeric(y, lower = 0, any.missing = FALSE)
				} else if (response_type == "ordinal") {
					assertIntegerish(y, lower = 1, any.missing = FALSE)
				}
			}
		},
		covariate_impute_if_necessary_and_then_create_model_matrix = function(){
			#make a copy... sometimes the raw will be the same as the imputed if there are no imputations
			private$Ximp = copy(private$Xraw)
			column_has_missingness = columns_have_missingness_cpp(private$Xraw)
			if (any(column_has_missingness)){
				if (private$missingness_method == "error"){
					if (should_run_asserts()) {
						missing_names = names(private$Xraw)[column_has_missingness]
						stop("Missing values detected in covariate(s): ",
							paste(missing_names, collapse = ", "),
							". Set missingness_method = \"impute\" or \"drop_column\" to handle missing data automatically.")
					}
				} else if (private$missingness_method == "drop_column"){
					cols_to_keep = which(!column_has_missingness)
					private$Ximp = private$Ximp[, ..cols_to_keep]
				} else {
					# "impute": random-forest imputation (missRanger with missForest fallback)
					#deal with include_is_missing_as_a_new_feature here
					if (private$include_is_missing_as_a_new_feature){
						missing_cols_idx = which(column_has_missingness)
						if (length(missing_cols_idx) > 0){
							# Use C++ function to create missingness indicators efficiently
							missingness_indicators = create_missingness_indicators_cpp(private$Ximp, missing_cols_idx)
							# Add the new columns to Ximp
							for (col_name in names(missingness_indicators)) {
								private$Ximp[[col_name]] = missingness_indicators[[col_name]]
							}
						}
					}
					#we need to convert characters into factor for the imputation to work
					col_types = get_column_types_cpp(private$Ximp)
					idx_cols_to_convert_to_factor = which(col_types == "character")
					private$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]
					#now do the imputation here by using missRanger (fast but fragile) and if that fails, use missForest (slow but more robust)
					private$Ximp = tryCatch({
											if (any(!is.na(private$y))){
												suppressWarnings(missRanger(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), verbose = FALSE, num.threads = self$num_cores)[, 1 : ncol(private$Ximp)])
											} else {
												suppressWarnings(missRanger(private$Ximp, verbose = FALSE, num.threads = self$num_cores))
											}
										}, error = function(e){
											if (any(!is.na(private$y))){
												suppressWarnings(missForest(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), num.threads = self$num_cores)$ximp[, 1 : ncol(private$Ximp)])
											} else {
												suppressWarnings(missForest(private$Ximp, num.threads = self$num_cores)$ximp)
											}
										}
									)
				}
			}
			analysis_col_names = names(private$Ximp)[!startsWith(names(private$Ximp), ".assignment_only_")]
			Ximp_for_model = if (length(analysis_col_names)) {
				private$Ximp[, ..analysis_col_names]
			} else {
				private$Ximp[, .SD, .SDcols = integer(0)]
			}
			#now let's drop any columns that don't have any variation
			num_unique_values_per_column = count_unique_values_cpp(Ximp_for_model)
			Ximp_for_model = Ximp_for_model[, .SD, .SDcols = which(num_unique_values_per_column > 1)]
			# now we need to update the numeric model matrix which may have expanded due to new factors, new missingness cols, etc
			private$X = create_model_matrix_from_features(private$model_formula, Ximp_for_model)
			# Ensure it is a numeric matrix (not character)
			if (should_run_asserts()) {
				if (ncol(private$X) > 0 && is.character(private$X)){
					stop("model.matrix returned a character matrix - this should not happen.")
				}
			}
			if (ncol(private$X) > 0){
				if (should_run_asserts()) {
					if (nrow(private$X) != nrow(private$Xraw) | nrow(private$X) != nrow(private$Ximp) | nrow(private$Ximp) != nrow(private$Xraw)){
						stop("improper sizing for the internal X representation")
					}
				}
			}
		},
		compute_all_subject_data = function(){
			i_present_y = which(!is.na(private$y))
			i_all = 1 : private$t
			i_all_y_present = intersect(i_all, i_present_y)
			
			# Cache lookup
			# Since covariates are fixed and NA positions in y are fixed during randomization,
			# the set of subjects with responses up to t is constant for a given t.
			cache_key = as.character(private$t)
			if (!is.null(private$all_subject_data_cache[[cache_key]])) {
				cpp_result = private$all_subject_data_cache[[cache_key]]
			} else {
				# Call consolidated C++ function for all matrix computations
				cpp_result = compute_all_subject_data_cpp(
					as.matrix(private$X[1:private$t, , drop = FALSE]),
					private$t,
					as.integer(i_all_y_present)
				)
				# Restore column names
				X_names = colnames(private$X)
				if (length(cpp_result$cols_prev) > 0) {
					nms = X_names[cpp_result$cols_prev]
					colnames(cpp_result$X_prev) = nms
					names(cpp_result$xt_prev) = nms
				}
				if (length(cpp_result$cols_all) > 0) {
					colnames(cpp_result$X_all) = X_names[cpp_result$cols_all]
				}
				if (length(cpp_result$cols_all_scaled) > 0) {
					nms = X_names[cpp_result$cols_all_scaled]
					colnames(cpp_result$X_all_scaled) = nms
					names(cpp_result$xt_all_scaled) = nms
				}
				if (length(cpp_result$cols_all_with_y_scaled) > 0) {
					colnames(cpp_result$X_all_with_y_scaled) = X_names[cpp_result$cols_all_with_y_scaled]
				}
				
				if (is.null(private$all_subject_data_cache)) private$all_subject_data_cache = list()
				private$all_subject_data_cache[[cache_key]] = cpp_result
			}
			# Add the simple array slices that don't need C++ optimization
			# These MUST NOT be cached because w and y change during randomization!
			cpp_result$w_all_with_y_scaled = private$w[i_all_y_present]
			cpp_result$y_all = private$y[i_all_y_present]
			cpp_result$dead_all = private$dead[i_all_y_present]
			cpp_result
		},
		assign_wt_Bernoulli = function(){
			rbinom(1, 1, private$prob_T)
		},
		has_private_method = function(method_name) {
			exists(method_name, envir = self$.__enclos_env__$private, inherits = FALSE)
		},
		maybe_set_seed = function() { if (!is.null(private$seed)) set.seed(private$seed) }
	)
)
