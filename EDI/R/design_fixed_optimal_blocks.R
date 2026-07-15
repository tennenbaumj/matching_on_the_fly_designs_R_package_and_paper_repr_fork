#' An Optimal-Blocks Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed experimental
#' design that first partitions subjects into \code{B} covariate-homogeneous blocks
#' by solving a balanced clustering problem with \pkg{ompr} and \pkg{glpk}, then randomizes
#' treatment within those blocks. When \code{B} is omitted and \code{n} is known
#' at initialization, the default is \code{floor(sqrt(n))}, truncated below at 1.
#'
#' @examples
#' des = DesignFixedOptimalBlocks$new(n = 9, response_type = 'continuous')
#' des$add_all_subjects_to_experiment(data.frame(x = rnorm(9)))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedOptimalBlocks = R6::R6Class("DesignFixedOptimalBlocks",
	inherit = DesignFixed,
	public = list(
			#' @description Returns TRUE so the framework pre-generates all w vectors
			#'   once per cell (paying the clustering cost once, reusing across reps).
			supports_batch_w_pregeneration = function() TRUE,
			#' @description Initialize a fixed optimal-blocks design.
			#' @param B Number of blocks to form. If omitted and \code{n} is supplied,
			#'   defaults to \code{floor(sqrt(n))}, with a minimum of 1.
			#' @param method Algorithm used to partition subjects into blocks.
			#'   \describe{
			#'     \item{\code{"K-way"} (default)}{Balanced k-means anticlustering via
			#'       \code{anticlust::balanced_clustering}.  Requires the
			#'       \pkg{anticlust} package.  Produces well-spread blocks and is
			#'       significantly faster than \code{"greedy"} (e.g., ~10x faster
			#'       for \eqn{n=200, p=10, B=10}) while achieving a better
			#'       within-block distance objective (e.g., ~4\% lower).}
			#'     \item{\code{"greedy"}}{Greedy nearest-neighbour matching
			#'       via \code{blockTools::block}.  Requires the \pkg{blockTools}
			#'       package.  Fast even for large \eqn{n}.}
			#'     \item{\code{"ompr"}}{Exact mixed-integer programme solved with
			#'       GLPK via \pkg{ompr}.  Globally optimal but scales as
			#'       \eqn{O(n^2 B)} in variables and is only practical for small
			#'       \eqn{n}.}
			#'   }
			#' @param dist Distance specification used only when \code{method = "ompr"}.
			#'   Either a function or one of \code{"euclidean"}, \code{"sum_abs_diff"},
			#'   or \code{"mahal"}. Default is \code{"mahal"}.
			#' @param response_type The response type for the design.
			#' @param prob_T Treatment assignment probability within each block.
			#' @param include_is_missing_as_a_new_feature Whether to include missingness indicators.
			#' @param n Planned sample size.
			#' @param verbose Whether to print progress messages.
			#' @param missingness_method How to handle missing values in covariates.
			#' @param design_formula A formula object.
			#' @param seed Integer seed for reproducibility.
			#' @return A new \code{DesignFixedOptimalBlocks} object.
			initialize = function(
				B = NULL,
				method = "K-way",
				dist = "mahal",
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
				if (should_run_asserts()) {
					assertChoice(method, c("ompr", "greedy", "K-way"))
				}
				if (should_run_asserts()) {
					if (method == "ompr") {
						if (!(is.function(dist) || (is.character(dist) && length(dist) == 1L)))
							stop("dist must be a function or one of 'euclidean', 'sum_abs_diff', or 'mahal'.")
						if (is.character(dist))
							assertChoice(dist, c("euclidean", "sum_abs_diff", "mahal"))
						assert_optimal_blocks_libraries_installed("DesignFixedOptimalBlocks with method='ompr'")
					} else if (method == "greedy") {
						assert_blocktools_installed("DesignFixedOptimalBlocks with method='greedy'")
					} else {
						assert_anticlust_installed("DesignFixedOptimalBlocks with method='K-way'")
					}
				}
				if (is.null(B)) {
					if (is.null(n)) {
						if (should_run_asserts()) {
							stop("DesignFixedOptimalBlocks requires B when n is not supplied.")
						}
					}
					B = max(1L, floor(sqrt(as.integer(n))))
				}
				if (should_run_asserts()) {
					assertCount(B, positive = TRUE)
				}
				super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
				private$blocking_capable = TRUE
				private$B      = as.integer(B)
				private$method = method
				private$dist_spec = dist
				private$uses_covariates = TRUE
				if (!is.null(n))
					if (should_run_asserts()) {
						private$assert_feasible_block_sizes(as.integer(n))
					}
			}
	),
	private = list(
		B = NULL,
		method = NULL,
		dist_spec = NULL,
		block_ids = NULL,
		distance_matrix = NULL,
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}
			block_ids = private$get_or_compute_block_ids()
			w_mat = replicate(r, as.numeric(as.character(randomizr::block_ra(blocks = block_ids, prob = private$prob_T))))
			storage.mode(w_mat) = "numeric"
			w_mat
		},
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			block_ids = private$get_or_compute_block_ids()
			if (is.null(bootstrap_type) || bootstrap_type == "within_blocks") {
				group_id = match(block_ids, unique(block_ids))
				list(i_b = stratified_bootstrap_indices_cpp(as.integer(group_id)), m_vec_b = NULL)
			} else {
				group_id = match(block_ids, unique(block_ids))
				i_b = resample_group_rows_cpp(as.integer(group_id), length(unique(group_id)))
				list(i_b = as.integer(i_b), m_vec_b = NULL)
			}
		},
		assert_feasible_block_sizes = function(n){
			private$assert_min_block_size(n, private$B)
			invisible(NULL)
		},
		get_or_compute_block_ids = function(){
			if (!is.null(private$block_ids)) {
				return(private$block_ids)
			}
			n = self$get_n()
			if (should_run_asserts()) {
				private$assert_feasible_block_sizes(n)
			}
			if (is.null(private$X)) {
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
			}
			X = private$X[seq_len(n), , drop = FALSE]
			if (ncol(X) == 0L) {
				private$block_ids = factor(rep(seq_len(private$B), length.out = n))
				return(private$block_ids)
			}
			private$block_ids = switch(private$method,
				greedy = private$solve_greedy_blocks(X),
				`K-way` = private$solve_kway_blocks(X),
				ompr = {
					D = private$get_or_compute_distance_matrix(X)
					private$solve_optimal_blocks(D)
				}
			)
			private$block_ids
		},
			get_or_compute_distance_matrix = function(X){
				if (!is.null(private$distance_matrix)) {
					return(private$distance_matrix)
				}
				if (is.function(private$dist_spec)) {
					D = optimal_blocks_distance_matrix_cpp(
						X = X,
						dist_code = 0L,
						dist_fn = private$dist_spec
					)
				} else {
					D = switch(private$dist_spec,
						euclidean = optimal_blocks_distance_matrix_cpp(X, dist_code = 1L),
						sum_abs_diff = optimal_blocks_distance_matrix_cpp(X, dist_code = 2L),
						mahal = optimal_blocks_distance_matrix_cpp(X, dist_code = 3L)
					)
				}
			storage.mode(D) = "double"
			private$distance_matrix = D
			D
		},
		solve_greedy_blocks = function(X) {
			n    = nrow(X)
			B    = private$B
			n_tr = as.integer(floor(n / B))
			df   = data.frame(`.__id__` = seq_len(n), as.data.frame(X), check.names = FALSE)
			block_out = blockTools::block(
				data       = df,
				n.tr       = n_tr,
				id.vars    = ".__id__",
				block.vars = colnames(X)
			)
			# $assg[[1]] is a data.frame: rows=blocks, cols=Treatment 1..n_tr + Max Distance
			assg_df   = blockTools::assignment(block_out)$assg[[1L]]
			id_cols   = grep("^Treatment", colnames(assg_df), value = TRUE)
			block_ids = integer(n)
			for (b in seq_len(nrow(assg_df))) {
				ids = as.integer(unlist(assg_df[b, id_cols]))
				ids = ids[!is.na(ids)]
				block_ids[ids] = b
			}
			# Assign any leftover subjects (n not divisible by n_tr) to nearest assigned neighbour
			unassigned = which(block_ids == 0L)
			if (length(unassigned) > 0L) {
				assigned = which(block_ids > 0L)
				for (i in unassigned) {
					dists = rowSums(sweep(X[assigned, , drop = FALSE], 2, X[i, ])^2)
					block_ids[i] = block_ids[assigned[which.min(dists)]]
				}
			}
			factor(block_ids, levels = seq_len(B))
		},
		solve_kway_blocks = function(X) {
			assignments = anticlust::balanced_clustering(X, K = private$B)
			factor(assignments, levels = seq_len(private$B))
		},
		solve_optimal_blocks = function(D){
			n = nrow(D)
			B = private$B
			lower_size = floor(n / B)
			upper_size = ceiling(n / B)
			model = ompr::MIPModel()
			model = ompr::add_variable(model, x[i, k], i = 1:n, k = 1:B, type = "binary")
			model = ompr::add_variable(model, z[i, j, k], i = 1:n, j = 1:n, k = 1:B, type = "binary", i < j)
			model = ompr::set_objective(
				model,
				ompr::sum_expr(D[i, j] * z[i, j, k], i = 1:n, j = 1:n, k = 1:B, i < j),
				"min"
			)
			model = ompr::add_constraint(model, ompr::sum_expr(x[i, k], k = 1:B) == 1, i = 1:n)
			model = ompr::add_constraint(model, ompr::sum_expr(x[i, k], i = 1:n) >= lower_size, k = 1:B)
			model = ompr::add_constraint(model, ompr::sum_expr(x[i, k], i = 1:n) <= upper_size, k = 1:B)
			model = ompr::add_constraint(model, x[k, k] == 1, k = 1:B)
			model = ompr::add_constraint(model, z[i, j, k] <= x[i, k], i = 1:n, j = 1:n, k = 1:B, i < j)
			model = ompr::add_constraint(model, z[i, j, k] <= x[j, k], i = 1:n, j = 1:n, k = 1:B, i < j)
			model = ompr::add_constraint(model, z[i, j, k] >= x[i, k] + x[j, k] - 1, i = 1:n, j = 1:n, k = 1:B, i < j)
			result = ompr::solve_model(model, ompr.roi::with_ROI(solver = "glpk", verbose = FALSE))
			solution = ompr::get_solution(result, x[i, k])
			if (should_run_asserts()) {
				if (nrow(solution) == 0L) {
					stop("ompr failed to produce a block assignment solution.")
				}
			}
			solution = solution[solution$value > 0.5, , drop = FALSE]
			if (should_run_asserts()) {
				if (nrow(solution) != n) {
					stop("ompr returned an incomplete block assignment.")
				}
			}
			labels = integer(n)
			labels[solution$i] = solution$k
			factor(labels, levels = seq_len(B))
		}
	)
)
