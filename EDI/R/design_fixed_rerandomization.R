#' A Rerandomization Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed rerandomization
#' experimental design.
#' This design generates random allocations and only accepts those that meet a
#' covariate balance criterion.
#' For balanced designs (prob_T = 0.5, even n) uses a native C++ parallel rejection sampler.
#'
#' @examples
#' \dontrun{
#' des = DesignFixedRerandomization$new(n = 10, response_type = 'continuous')
#' }
#' @export
DesignFixedRerandomization = R6::R6Class("DesignFixedRerandomization",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a rerandomization fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T  The probability of the treatment assignment.
		#' @param obj_val_cutoff 	The maximum allowable objective value. Cannot be specified together with prop_acceptable.
		#' @param prop_acceptable	The proportion of randomizations to accept (draws r/prop_acceptable total, returns r lowest). Cannot be specified together with obj_val_cutoff.
		#' @param objective 	The objective function to use. Default is "mahal_dist".
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedRerandomization` object
		#'
		initialize = function(
				response_type,
				prob_T = 0.5,
				obj_val_cutoff = NULL,
				prop_acceptable = NULL,
				objective = "mahal_dist",
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
			if (!is.null(obj_val_cutoff) && !is.null(prop_acceptable)) {
				stop("Cannot specify both obj_val_cutoff and prop_acceptable.")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$obj_val_cutoff = obj_val_cutoff
			private$prop_acceptable = prop_acceptable
			private$objective = objective
			private$uses_covariates = TRUE
			#note: we are not setting private$m as this is not a blocking design
		}
	),
	private = list(
		obj_val_cutoff = NULL,
		prop_acceptable = NULL,
		objective = NULL,
		S_inv = NULL,
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				n_T = round(n * private$prob_T)
				return(replicate(r, sample(c(rep(1, n_T), rep(0, n - n_T)))))
			}
			# prop_acceptable path: draw r/prop_acceptable randomizations, return r with lowest obj vals
			if (!is.null(private$prop_acceptable)) {
				if (private$objective == "mahal_dist" && is.null(private$S_inv)){
					X = private$X[1:n, , drop = FALSE]
					S = var(X)
					if (abs(det(S)) < 1e-10){
						S = S + diag(1e-6, ncol(X))
					}
					private$S_inv = solve(S)
				}
				X = private$X[1:n, , drop = FALSE]
				n_draw = round(r / private$prop_acceptable)
				ext_seed = if (!is.null(private$seed)) private$seed else sample.int(.Machine$integer.max, 1L)
				if (private$prob_T == 0.5) {
					indicTs = complete_randomization_forced_balanced_cpp(n, n_draw, ext_seed)
				} else {
					n_T = round(n * private$prob_T)
					indicTs = complete_randomization_imbalanced_cpp(n, n_T, n_draw, ext_seed)
				}
				obj_vals = compute_objective_vals_cpp(X, indicTs, private$objective, private$S_inv)
				ord = order(obj_vals)
				# indicTs is n_draw x n; subset rows for best r, then transpose to n x r
				w_mat = t(indicTs[ord[seq_len(r)], , drop = FALSE])
				storage.mode(w_mat) = "numeric"
				return(w_mat)
			}
			# C++ parallel rejection sampler for balanced even-n case
			if (private$prob_T == 0.5 && n %% 2 == 0) {
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
				X = private$X[1:n, , drop = FALSE]
				cutoff_user = if (is.null(private$obj_val_cutoff)) Inf else private$obj_val_cutoff
				# The C++ function uses the M-matrix objective which is scaled relative to the
				# user-facing objective: mahal_dist → f_cpp = user_mahal/4;
				# abs_sum_diff → f_cpp = GED-standardised/2.
				cutoff_scale = if (private$objective == "mahal_dist") 4.0 else 2.0
				cutoff_cpp   = if (is.infinite(cutoff_user)) Inf else cutoff_user / cutoff_scale
				max_draws = max(r * 1000L, 100000L)
				w_mat = rerandomization_search_cpp(
					X_raw     = X,
					r         = as.integer(r),
					objective = private$objective,
					cutoff    = as.double(cutoff_cpp),
					max_draws = as.integer(max_draws)
				)
				w_mat = private$validate_allocation_matrix(w_mat, n = n, r = ncol(w_mat), require_balanced = TRUE)
				# Recycle if fewer were found than requested (tight cutoff)
				if (ncol(w_mat) < r) {
					w_mat = w_mat[, rep(seq_len(ncol(w_mat)), length.out = r), drop = FALSE]
				}
				storage.mode(w_mat) = "numeric"
				return(w_mat[, seq_len(r), drop = FALSE])
			}
			# Fallback pure-R for unbalanced or odd-n cases
			if (private$objective == "mahal_dist" && is.null(private$S_inv)){
				X = private$X[1:n, , drop = FALSE]
				S = var(X)
				if (abs(det(S)) < 1e-10){
					S = S + diag(1e-6, ncol(X))
				}
				private$S_inv = solve(S)
			}
			w_mat = matrix(NA_real_, nrow = n, ncol = r)
			for (j in seq_len(r)){
				w_mat[, j] = private$generate_one_rerandomized_w()
			}
			w_mat
		},
		validate_allocation_matrix = function(w_mat, n, r, require_balanced = FALSE){
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (should_run_asserts()) {
				if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < 1L) {
					stop("DesignFixedRerandomization returned an unexpected allocation matrix shape.")
				}
			}
			storage.mode(w_mat) = "numeric"
			if (should_run_asserts()) {
				if (any(!is.finite(w_mat)) || any(is.na(w_mat))) {
					stop("DesignFixedRerandomization returned non-finite treatment assignments.")
				}
				if (any(!(w_mat %in% c(0, 1)))) {
					stop("DesignFixedRerandomization returned an invalid treatment assignment matrix.")
				}
				if (isTRUE(require_balanced)) {
					treated_counts = colSums(w_mat)
					if (any(treated_counts != n / 2)) {
						stop("DesignFixedRerandomization returned an unbalanced allocation.")
					}
				}
			}
			w_mat[, seq_len(min(r, ncol(w_mat))), drop = FALSE]
		},
		generate_one_rerandomized_w = function(){
			n = self$get_n()
			X = private$X[1:n, , drop = FALSE]
			n_T = round(n * private$prob_T)
			repeat {
				if (private$prob_T == 0.5){
					w_cand = sample(c(rep(1, n_T), rep(0, n - n_T)))
				} else {
					w_cand = rbinom(n, 1, private$prob_T)
				}
				obj_val = private$compute_obj(X, w_cand)
				if (is.null(private$obj_val_cutoff) || obj_val <= private$obj_val_cutoff){
					return(w_cand)
				}
			}
		},
		compute_obj = function(X, w){
			if (private$objective == "mahal_dist"){
				diff = colMeans(X[w == 1, , drop = FALSE]) - colMeans(X[w == 0, , drop = FALSE])
				return(as.numeric(diff %*% private$S_inv %*% diff))
			} else if (private$objective == "abs_sum_diff"){
				diff = colMeans(X[w == 1, , drop = FALSE]) - colMeans(X[w == 0, , drop = FALSE])
				return(sum(abs(diff)))
			}
			stop("Unsupported objective")
		}
	)
)
