#' Sequential One-by-One Experimental Design
#'
#' Abstract R6 class encapsulating data and functionality for a sequential one-by-
#' one experimental design.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' seq_des = DesignSeqOneByOne$new(n = 6, response_type = 'continuous')
#' seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
DesignSeqOneByOne = R6::R6Class("DesignSeqOneByOne",
	lock_objects = FALSE,
	inherit = DesignMatching,
	public = list(
		#' @description Initialize a sequential one-by-one design.
		#' @param response_type The data type of response values.
		#' @param prob_T The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature If missing data is present, include a dummy
		#'   variable for it.
		#' @param n The sample size (if fixed).
		#' @param verbose Whether to print progress messages.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#' @param ... Extra arguments passed to the \code{Design} superclass.
		initialize = function(
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL,
				...
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed, ...)
			private$maybe_set_seed()
			private$equal_block_sizes = FALSE
		},
		#' @description Checks if the design supports resampling.
		#' @return Always TRUE for sequential designs.
		supports_resampling = function(){
			TRUE
		},
		#' @description Add subject-specific measurements for the next subject entrant.
		#'
		#' @param x_new A data frame with one row representing the new subject's covariates.
		#' @param allow_new_cols Allow new features in the new subject's covariates.
		add_one_subject = function(x_new, allow_new_cols = TRUE){
			if (should_run_asserts()) {
				assertClass(x_new, "data.frame")
			}
			x_new = as.data.table(x_new)
			if (should_run_asserts()) {
				if (nrow(x_new) != 1){
					stop("You can only add one subject at a time.")
				}
			
				if (private$t == 0 && !is.null(private$strata_cols)) {
					for (col in private$strata_cols) {
						if (is.numeric(x_new[[col]])) {
							stop("Error: Continuous covariates are not allowed for stratification in sequential designs because stable binning cannot be determined on-the-fly. Please pre-discretize the numeric column(s) into factors/categories (e.g., using fixed clinical thresholds) before adding subjects to the experiment.")
						}
					}
				}
			}
			j_with_NAs = is.na(unlist(x_new))
			if (any(j_with_NAs) & private$t == 0){
				x_new = x_new[which(!j_with_NAs)]
				if (!allow_new_cols){
					warning("There is missing data in the first subject's covariate value(s). Setting the flag allow_new_cols = FALSE will disallow additional subjects")
				}
			}
			xnew_data_types = get_column_types_cpp(x_new)
			if (should_run_asserts()) {
				if ("ordered" %in% xnew_data_types){
					stop("Ordered factor data type is not supported; please convert to either an unordered factor or numeric.")
				}
				if ("Date" %in% xnew_data_types){
					stop("Date data type is not supported; please convert to numeric.")
				}
			}
			if (private$t > 0){
				Xraw_data_types = get_column_types_cpp(private$Xraw)
				colnames_Xraw = names(private$Xraw)
				colnames_xnew = names(x_new)
				if (setequal(colnames_Xraw, colnames_xnew)){
					if (should_run_asserts()) {
						idx_data_types_that_changed = which(xnew_data_types != Xraw_data_types)
						if (length(idx_data_types_that_changed) > 0){
							for (e in idx_data_types_that_changed){
								warning("You entered data type ", xnew_data_types[e], " for attribute named ", colnames_Xraw[e], " that was previously entered with data type ", Xraw_data_types[e])
							}
						}
					}
				} else {
					if (allow_new_cols){ #make NA's in appropriate places
						new_Xraw_cols = setdiff(colnames_xnew, colnames_Xraw)
						if (length(new_Xraw_cols) > 0){
							new_Xraw_col_types = xnew_data_types[new_Xraw_cols]
							for (j in 1 : length(new_Xraw_cols)){
								private$Xraw[, (new_Xraw_cols[j]) := switch(new_Xraw_col_types[j],
									character = NA_character_,
									factor =    NA_character_, #I think this is correct
									numeric =   NA_real_,
									logical =   NA_real_, #just let it be zero or one
									integer =   NA_real_  #I don't want to take the risk on a decimal popping up somewhere
								)]
							}
						}
						new_xnew_cols = setdiff(colnames_Xraw, colnames_xnew)
						if (length(new_xnew_cols) > 0){
							new_xnew_cols_types = Xraw_data_types[new_xnew_cols]
							for (j in 1 : length(new_xnew_cols)){
								x_new[, (new_xnew_cols[j]) := switch(new_xnew_cols_types[j],
									character = NA_character_,
									factor =    NA_character_, #I think this is correct
									numeric =   NA_real_,
									logical =   NA_real_, #just let it be zero or one
									integer =   NA_real_  #I don't want to take the risk on a decimal popping up somewhere
								)]
							}
						}
					} else {
						stop(paste(
							"The new subject vector has columns:\n  ",
							paste(colnames_xnew, collapse = ", "),
							"\nwhich are not the same as the current dataset's columns:\n  ",
							paste(colnames_Xraw, collapse = ", "),
							"\nIf you want to allow new columns on-the-fly, run this function again with the option\n  'allow_new_cols = TRUE'"
						))
					}
				}
			}
			#add new subject's measurements to the raw data frame (there should be the same exact columns even if there are new ones introduced)
			private$Xraw = rbindlist(list(private$Xraw, x_new))
			private$p_raw_t = ncol(private$Xraw)
			#iterate t
			private$t = private$t + 1L #t must be an integer for data.table's fast "set" function below to work
			#we only bother with imputation and model matrices if we have enough data otherwise it's a huge mess
			#thus, designs cannot utilize imputations nor model matrices until this condition is met
			#luckily, those are the designs implemented herein so we have complete control (if you are extending this package, you'll have to deal with this issue here)
			if (private$t > (ncol(private$Xraw) + 2) & private$uses_covariates){ #we only need to impute if we need the X's to make the allocation decisions
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
			}
		},
		#' @description Adds a subject and assigns treatment.
		#' @param x_new A data frame with one row representing the new subject's covariates.
		#' @return The treatment assignment as {-1,+1} (+1 = treated, -1 = control).
		add_one_subject_to_experiment_and_assign = function(x_new){
			self$add_one_subject(x_new)
			w_t = self$assign_wt()
			if (private$fixed_sample){
				private$w[private$t] = w_t
			} else {
				if (length(private$w) < private$t) {
					private$w = c(private$w, w_t)
				} else {
					private$w[private$t] = w_t
				}
			}
			2L * private$w[private$t] - 1L
		},
		#' @description Assigns treatment to the current subject.
		#' @return The treatment assignment (0 or 1).
		assign_wt = function(){
			stop("Must be implemented by subclass.")
		},
		#' @description Prints the current subject's assignment.
		print_current_subject_assignment = function(){
			cat("Subject number", private$t, "is assigned to", ifelse(private$w[private$t] == 1, "TREATMENT", "CONTROL"), "via design", class(self)[1], "\n")
		}
	),
	private = list(
		draw_ws_raw = function(r = 100){
			w_mat = matrix(NA_real_, nrow = private$t, ncol = r)
			for (j in 1 : r){
				if (private$fixed_sample) {
					private$w = rep(NA_real_, private$n)
				} else {
					private$w = rep(NA_real_, private$t)
				}
				for (t in 1 : private$t){
					private$w[t] = self$assign_wt()
				}
				w_mat[, j] = private$w[1:private$t]
			}
			w_mat
		}
	)
)
