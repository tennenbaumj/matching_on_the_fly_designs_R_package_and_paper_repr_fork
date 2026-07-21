#' Pattern-1 mixin: KK compound pass-through behavior
#'
#' Provides the five private methods shared by all KK compound estimators:
#' \code{supports_likelihood_tests}, \code{compute_basic_match_data},
#' \code{compute_estimate_from_matched_and_reservoir}, \code{only_matches},
#' \code{only_reservoir}, and \code{compute_reservoir_and_match_statistics}.
#'
#' Splice into a daughter class via
#' \code{private = c(InferenceMixinKKPassThroughCompound$private, list(...))}.
#' The capability flag \code{private$kk_passthrough_compound} is set to \code{TRUE}.
#'
#' @keywords internal
#' @noRd
InferenceMixinKKPassThroughCompound = list(
	public = list(),
	private = list(
		kk_passthrough_compound = TRUE,

		supports_likelihood_tests = function() FALSE,

		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),

		compute_estimate_from_matched_and_reservoir = function(run_matched, run_reservoir){
			if (!isTRUE(private$has_match_structure)) {
				private$cache_nonestimable_estimate("kk_design_required")
				return(invisible(NULL))
			}
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			if (is.null(private$cached_values$KKstats)) {
				return(invisible(NULL))
			}
			if (private$only_matches()){
				run_matched()
			} else if (private$only_reservoir()){
				run_reservoir()
			} else {
				run_matched()
				run_reservoir()
			}
		},

		only_matches = function(){
			if (is.null(private$cached_values$KKstats)) return(FALSE)
			nRT = private$cached_values$KKstats$nRT
			nRC = private$cached_values$KKstats$nRC
			if (!is.finite(nRT) || !is.finite(nRC)) return(FALSE)
			nRT <= 1 || nRC <= 1
		},

		only_reservoir = function(){
			if (is.null(private$cached_values$KKstats)) return(FALSE)
			m = private$cached_values$KKstats$m
			is.finite(m) && m <= 1
		},

		compute_reservoir_and_match_statistics = function(){
			nRC = private$cached_values$KKstats$nRC
			nRT = private$cached_values$KKstats$nRT
			nR = nRT + nRC
			m = private$cached_values$KKstats$m

			y_reservoir_T = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 1]
			y_reservoir_C = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 0]

			ssqD_bar = if (is.finite(m) && m > 1){
							var_cpp(private$cached_values$KKstats$y_matched_diffs) / m
						} else {
							NA_real_
						}
			ssqR = if (is.finite(nRT) && is.finite(nRC) && nRT > 1 && nRC > 1 && nR > 2){
						(var_cpp(y_reservoir_T) * (nRT - 1) + var_cpp(y_reservoir_C) * (nRC - 1)) /
							(nR - 2) * (1 / nRT + 1 / nRC)
					} else {
						NA_real_
					}

			private$cached_values$KKstats$d_bar = if (is.finite(m) && m > 0) mean_cpp(private$cached_values$KKstats$y_matched_diffs) else {
				NA_real_
			}
			private$cached_values$KKstats$ssqD_bar = ssqD_bar
			private$cached_values$KKstats$r_bar = if (is.finite(nRT) && is.finite(nRC) && nRT > 0 && nRC > 0) mean_cpp(y_reservoir_T) - mean_cpp(y_reservoir_C) else {
				NA_real_
			}
			private$cached_values$KKstats$ssqR = ssqR
			private$cached_values$KKstats$w_star = if (is.finite(ssqR) && is.finite(ssqD_bar) && (ssqR + ssqD_bar) > 0) {
				ssqR / (ssqR + ssqD_bar)
			} else {
				NA_real_
			}
		}
	)
)
