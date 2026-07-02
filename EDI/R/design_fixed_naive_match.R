#' A Naive Match Fixed Design
#'
#' An R6 Class encapsulating a dummy fixed matched-pair experimental design.
#' Extends \code{DesignFixedBinaryMatch} but skips the (expensive)
#' \pkg{nbpMatching} non-bipartite optimal matching solve entirely: subjects
#' are paired consecutively as they arrive, i.e. \code{<1,2>, <3,4>, ...,
#' <n-1,n>}, with no regard to covariate distance. Useful in simulations
#' where only a valid matched-pair data structure is needed (e.g. to exercise
#' the KK-prefixed paired inference classes) and the quality of the matching
#' itself is not under study, since it avoids both the KK21 sequential
#' weight-fitting cost and the \pkg{nbpMatching} solve.
#'
#' @examples
#' des = DesignFixedNaiveMatch$new(n = 10, response_type = 'continuous')
#' des$add_all_subjects_to_experiment(data.frame(x1 = rnorm(10)))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedNaiveMatch = R6::R6Class("DesignFixedNaiveMatch",
	inherit = DesignFixedBinaryMatch,
	private = list(
		ensure_matching_structure_computed = function(){
			n = self$get_n()
			if (is.null(private$bms)) {
				if (should_run_asserts()) {
					if (n %% 2L != 0L) {
						stop("Design matrix must have an even number of rows for naive matching.")
					}
				}
				pairs = matrix(seq_len(n), ncol = 2L, byrow = TRUE)
				private$bms = list(
					indicies_pairs = pairs,
					indices_pairs  = pairs,
					n = n,
					p = if (is.null(private$X)) 0L else ncol(private$X)
				)
				m_vec = integer(n)
				for (i in seq_len(nrow(pairs))){
					m_vec[pairs[i, 1]] = i
					m_vec[pairs[i, 2]] = i
				}
				private$m = m_vec
				private$reset_matching_caches()
			}
			invisible(NULL)
		}
	)
)
