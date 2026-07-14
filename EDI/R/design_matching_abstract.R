#' An Abstract Matching Experimental Design
#'
#' @name DesignMatching
#' @description Internal method.
#' An abstract R6 class encapsulating shared matching-specific caches and
#' utilities for designs that may expose matched-pair / reservoir structure.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' des = DesignMatching$new(n = 6, response_type = "continuous")
#' }
DesignMatching = R6::R6Class("DesignMatching",
	lock_objects = FALSE,
	inherit = DesignBlocking,
	public = list(
		#' @description Check whether this design currently has matching structure.
		#'
		#' @return \code{TRUE} if the design advertises matching support.
		is_matching_design = function(){
			isTRUE(private$matching_capable)
		},
		#' @description Assert that this design supports matching-specific operations.
		assert_matching_design = function(){
			if (should_run_asserts()) {
				if (!self$is_matching_design()) {
					stop("This design requires a matching design.")
				}
			}
		},
		#' @description Return cluster IDs implied by the current matching structure.
		#'
		#' @param m_vec Optional integer match vector. Defaults to this design's match vector.
		#'
		#' @return Integer cluster IDs for matched pairs plus singleton reservoir units.
		get_matching_cluster_ids = function(m_vec = private$m){
			private$compute_matching_cluster_ids(m_vec)
		}
	),
	private = list(
		draw_ws_raw = function(r = 100){
			stop("draw_ws_raw must be implemented by a concrete design subclass.")
		},
		xm_structural     = NULL,
		xm_m_vec          = NULL,
		lin_xm_structural = NULL,
		lin_xm_m_vec      = NULL,
		cluster_id        = NULL,
		cluster_id_m_vec  = NULL,
		boot_pair_rows    = NULL,
		boot_i_reservoir  = NULL,
		boot_n_reservoir  = NULL,
		matching_capable  = FALSE,
		ensure_matching_structure_computed = function(){
			invisible(NULL)
		},
		reset_matching_caches = function(){
			private$xm_structural     = NULL
			private$xm_m_vec          = NULL
			private$lin_xm_structural = NULL
			private$lin_xm_m_vec      = NULL
			private$cluster_id        = NULL
			private$cluster_id_m_vec  = NULL
			private$boot_pair_rows    = NULL
			private$boot_i_reservoir  = NULL
			private$boot_n_reservoir  = NULL
			invisible(NULL)
		},
		init_matching_bootstrap_structure = function(){
			if (!is.null(private$boot_pair_rows)) return(invisible(NULL))
			m_vec = private$m
			n = private$n
			if (is.null(m_vec)){
				private$boot_i_reservoir = seq_len(n)
				private$boot_n_reservoir = n
				private$boot_pair_rows   = matrix(integer(0), nrow = 0L, ncol = 2L)
				return(invisible(NULL))
			}
			m_vec_int = as.integer(m_vec)
			m_vec_int[is.na(m_vec_int)] = 0L
			i_reservoir = which(m_vec_int == 0L)
			m_max = max(m_vec_int)
			pair_rows = if (m_max > 0L) {
				pr = matrix(integer(0), nrow = m_max, ncol = 2L)
				for (pid in seq_len(m_max)) pr[pid, ] = which(m_vec_int == pid)
				pr
			} else {
				matrix(integer(0), nrow = 0L, ncol = 2L)
			}
			private$boot_i_reservoir = i_reservoir
			private$boot_n_reservoir = length(i_reservoir)
			private$boot_pair_rows   = pair_rows
			invisible(NULL)
		},
		draw_matching_bootstrap_indices = function(){
			private$init_matching_bootstrap_structure()
			draw_matching_bootstrap_sample_cpp(
				i_reservoir = private$boot_i_reservoir,
				pair_rows   = private$boot_pair_rows,
				n_reservoir = private$boot_n_reservoir
			)
		},
		compute_matching_cluster_ids = function(m_vec = private$m){
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec_int = as.integer(m_vec)
			m_vec_int[is.na(m_vec_int)] = 0L
			des_m = private$m
			if (is.null(des_m)) des_m = rep(NA_integer_, private$n)
			des_m_int = as.integer(des_m)
			des_m_int[is.na(des_m_int)] = 0L
			if (!is.null(private$cluster_id) && identical(m_vec_int, des_m_int)){
				return(private$cluster_id)
			}
			cluster_id = compute_cluster_ids_cpp(m_vec_int)
			if (identical(m_vec_int, des_m_int)){
				private$cluster_id = cluster_id
				private$cluster_id_m_vec = m_vec_int
			}
			cluster_id
		},
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			private$ensure_matching_structure_computed()
			if (!isTRUE(private$matching_capable) || is.null(private$m)){
				n = self$get_n()
				return(list(i_b = sample_int_replace_cpp(n, n), m_vec_b = NULL))
			}
			private$draw_matching_bootstrap_indices()
		}
	)
)
