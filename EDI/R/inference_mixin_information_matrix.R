#' Mixin for Information-Matrix Inference
#'
#' A Pattern-1 mixin (plain list with code{$public} and code{$private} slots)
#' providing information-source selection and treatment-coefficient standard
#' errors for likelihood-backed inference. Consumers must provide
#' code{get_likelihood_test_spec()} and code{get_default_information_source()}
#' private methods, as well as code{information_preference} and
#' code{information_source_used} private fields.
#'
#' Splice into a class with
#' code{private = c(InferenceMixinInformationMatrix$private, list(...))}.
#'
#' @keywords internal
#' @noRd
InferenceMixinInformationMatrix = list(
	public = list(),
	private = list(
		get_information_matrix = function(spec = NULL, fit = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec)) return(NULL)
			if (is.null(fit)) {
				fit = spec$full_fit %||% private$cached_mod
			}
			if (is.null(fit)) return(NULL)

			extract_fisher = function(){
				tryCatch({
					if (!is.null(spec$fisher_information)) {
						spec$fisher_information(fit)
					} else if (!is.null(fit$fisher_information)) {
						fit$fisher_information
					} else if (identical(fit$information_type %||% "", "fisher") && !is.null(fit$information)) {
						fit$information
					} else {
						NULL
					}
				}, error = function(e) NULL)
			}

			extract_observed = function(){
				tryCatch({
					if (!is.null(spec$observed_information)) {
						spec$observed_information(fit)
					} else if (!is.null(fit$observed_information)) {
						fit$observed_information
					} else if (identical(fit$information_type %||% "", "observed") && !is.null(fit$information)) {
						fit$information
					} else {
						NULL
					}
				}, error = function(e) NULL)
			}

			extract_legacy = function(){
				tryCatch({
					if (!is.null(spec$information)) {
						spec$information(fit)
					} else {
						fit$information
					}
				}, error = function(e) NULL)
			}

			preference = private$information_preference
			if (identical(preference, "auto")) {
				preference = private$get_default_information_source()
			}
			if (identical(preference, "fisher")) {
				information = extract_fisher()
				if (is.null(information)) {
					stop(class(self)[1], " does not expose Fisher information for information-backed inference.", call. = FALSE)
				}
				private$information_source_used = "fisher"
				return(information)
			}

			if (identical(preference, "observed")) {
				information = extract_observed()
				if (is.null(information)) {
					stop(class(self)[1], " does not expose observed information for information-backed inference.", call. = FALSE)
				}
				private$information_source_used = "observed"
				return(information)
			}

			information = extract_fisher()
			if (!is.null(information)) {
				private$information_source_used = "fisher"
				return(information)
			}
			information = extract_observed()
			if (!is.null(information)) {
				private$information_source_used = "observed"
				return(information)
			}
			information = extract_legacy()
			if (!is.null(information)) {
				private$information_source_used = "legacy"
			}
			information
		},

		compute_variance_from_information_matrix = function(information, j){
			information = as.matrix(information)
			if (!is.matrix(information) || nrow(information) != ncol(information) || length(j) != 1L ||
				!is.finite(j) || j < 1L || j > nrow(information)) {
				return(NA_real_)
			}
			if (nrow(information) == 1L) {
				val = as.numeric(information[1L, 1L])
				return(if (is.finite(val) && val > 0) 1 / val else NA_real_)
			}
			res = tryCatch(eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(information, as.integer(j)), error = function(e) NA_real_)
			if (is.finite(res) && res > 0) return(res)

			vcov = tryCatch(solve(information), error = function(e) NULL)
			if (is.null(vcov) || any(!is.finite(vcov))) {
				vcov = tryCatch(qr.solve(information, diag(nrow(information))), error = function(e) NULL)
			}
			if (is.null(vcov) || any(!is.finite(vcov))) return(NA_real_)
			vcov = (vcov + t(vcov)) / 2
			as.numeric(vcov[j, j])
		},

		compute_standard_error_from_information_matrix = function(spec = NULL, fit = NULL, j = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec)) return(NA_real_)
			if (is.null(j)) {
				j = as.integer(spec$j)
			}
			information = tryCatch(private$get_information_matrix(spec = spec, fit = fit), error = function(e) NULL)
			if (is.null(information)) return(NA_real_)
			variance = private$compute_variance_from_information_matrix(information, j)
			if (is.finite(variance) && variance >= 0) sqrt(variance) else NA_real_
		},

		get_score_test_information_matrix = function(spec, fit){
			private$get_information_matrix(spec = spec, fit = fit)
		}
	)
)
