#' A Fixed Design
#'
#' An abstract R6 Class encapsulating the data and functionality for a fixed experimental design.
#' This class takes care of whole-experiment randomization.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' des = DesignFixed$new(n = 10, response_type = 'continuous')
#' des$add_all_subjects_to_experiment(data.frame(x1 = rnorm(10)))
#' des$assign_w_to_all_subjects()
#' }
DesignFixed = R6::R6Class("DesignFixed",
	lock_objects = FALSE,
	inherit = DesignMatching,
	public = list(
		#' @description Initialize a fixed experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param  prob_T  Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param  n  		The sample size.
		#' @param verbose A flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param model_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#' @param ... Extra arguments passed to the \code{Design} superclass.
		#'
		#' @return  A new `DesignFixed` object
		initialize = function(
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				model_formula = ~ .,
				seed = NULL,
				...
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, model_formula, seed = seed, ...)
		},
		#' @description Assign treatment to all subjects in the fixed experiment.
		#' @param w_precomputed Optional {-1,+1} numeric vector of length n. If supplied the
		#'   allocation is used directly (converted to internal {0,1} storage) and
		#'   \code{draw_ws_according_to_design} is not called (avoids e.g. the Java
		#'   round-trip for \code{DesignFixedGreedy}).
		assign_w_to_all_subjects = function(w_precomputed = NULL){
			if (!is.null(w_precomputed)) {
				private$w[1:self$get_n()] = (as.numeric(w_precomputed) + 1L) / 2L
			} else {
				private$w[1:self$get_n()] = private$draw_ws_raw(1)[, 1]
			}
		},
		#' @description Add all subjects' covariates to a fixed design at once.
		#'
		#' @param X_all A data frame containing the full covariate matrix.
		#' @return Invisibly returns the design object.
		add_all_subjects_to_experiment = function(X_all){
			n_all = nrow(X_all)
			n_expected = self$get_n()
			if (should_run_asserts()) {
				assertClass(X_all, "data.frame")
				self$assert_fixed_sample()
				if (n_all != n_expected){
					stop("X_all must have exactly ", n_expected, " rows for this fixed design.")
				}
				if (private$t > 0L || nrow(private$Xraw) > 0L){
					stop("Subjects have already been added to this design.")
				}
			}
			private$Xraw = copy(as.data.table(X_all))
			private$p_raw_t = ncol(private$Xraw)
			private$t = n_all
			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			invisible(self)
		},
		#' @description Add all subject responses for a fixed design.
		#'
		#' @param ys The responses as a numeric vector.
		#' @param deads The binary vector indicating if dead/censored.
		add_all_subject_responses = function(ys, deads = NULL){
			if (is.null(deads)){
				deads = rep(1, private$t)
			}
			if (should_run_asserts()) {
				self$assert_fixed_sample()
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
		#' @description Overwrite all subject assignments for a fixed design.
		#'
		#' @param w A {-1,+1} vector of subject assignments (+1 = treated, -1 = control).
		overwrite_all_subject_assignments = function(w){
			if (should_run_asserts()) {
				assertIntegerish(w, lower = -1, upper = 1, any.missing = FALSE, len = private$t)
				if (any(!(w %in% c(-1L, 1L)))) {
					stop("overwrite_all_subject_assignments: w must contain only -1 (control) or +1 (treated).")
				}
			}
			private$w = (as.numeric(w) + 1L) / 2L
		},
		#' @description Check if the design supports resampling.
		#'
		#' @return 	TRUE if supported.
		supports_resampling = function(){
			class(self)[1] != "DesignFixed"
		}
	),
	private = list()
)
