#' An Abstract Blocking Experimental Design
#'
#' @name DesignBlocking
#' @description Internal method.
#' An abstract R6 Class encapsulating shared blocking state and functionality for
#' experimental designs that may define block identifiers.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' des = DesignBlocking$new(n = 6, response_type = "continuous")
#' }
DesignBlocking = R6::R6Class("DesignBlocking",
	lock_objects = FALSE,
	inherit = Design,
	public = list(
		#' @description Checks whether this design currently has blocking structure.
		#'
		#' @return \code{TRUE} if the design supports blocking or manual block IDs
		#'   have been recorded, \code{FALSE} otherwise.
		is_blocking_design = function(){
			isTRUE(private$blocking_capable) || !is.null(private$m)
		},
		#' @description Inject a precomputed matrix of assignment vectors for CMH SE estimation.
		#' @param w_mat Integer matrix (n x se_est_num_vectors).
		inject_cmh_se_w_mat = function(w_mat){
			private$cmh_se_w_mat = w_mat
			invisible(self)
		},
		#' @description Retrieve the precomputed CMH SE w_mat, or NULL if not set.
		get_cmh_se_w_mat = function() private$cmh_se_w_mat,
		#' @description Checks whether this design is a blocking design.
		assert_blocking_design = function(){
			if (should_run_asserts()) {
				if (!self$is_blocking_design()){
					stop("This design requires a blocking design.")
				}
			}
		},
		#' @description Check whether the current blocking structure is complete.
		#'
		#' @return \code{TRUE} if blocking is defined and all block IDs are positive.
		is_complete_blocking_design = function(){
			if (!self$is_blocking_design()) {
				return(FALSE)
			}
			block_ids = self$get_block_ids()
			length(block_ids) > 0L && all(is.finite(block_ids)) && all(block_ids > 0L)
		},
		#' @description Checks whether all blocks have the same number of subjects.
		#'
		#' @return `TRUE` invisibly if the block sizes are equal.
		assert_equal_block_sizes = function(){
			self$assert_blocking_design()
			block_ids = self$get_block_ids()
			block_sizes = as.integer(table(block_ids))
			if (should_run_asserts()) {
				if (length(block_sizes) > 1L && any(block_sizes != block_sizes[1L])) {
					stop("All blocks must have the same number of subjects.")
				}
			}
			invisible(TRUE)
		},
		#' @description Record externally-supplied matched-pair or block identifiers.
		#'
		#' This is primarily for post-hoc analysis of completed experiments whose
		#' block IDs were defined outside the built-in design classes.
		#'
		#' @param m A positive-integer vector of block IDs with one entry per subject.
		#'
		#' @return Invisibly returns \code{self}.
		add_all_subject_matched_pair_ids = function(m){
			if (should_run_asserts()) {
				assertIntegerish(m, lower = 1, any.missing = FALSE, len = private$t)
			}
			private$m = as.integer(m)
			invisible(self)
		},
		#' @description Set the strata (block identifiers) for this design.
		#'
		#' Can only be called when \code{private$m} is still \code{NULL} (i.e. before strata
		#' have been derived or assigned). The design must advertise blocking support.
		#'
		#' @param m A non-negative-integer vector of block IDs, one entry per subject already
		#'   recorded. Zero is reserved as a reservoir flag used by sequential KK designs.
		#'
		#' @return Invisibly returns \code{self} for method chaining.
		set_m = function(m){
			if (should_run_asserts()) {
				if (!isTRUE(private$blocking_capable)){
					stop("set_m() can only be used with a blocking-capable design.")
				}
				if (!is.null(private$m)){
					stop("set_m() can only be called when the strata have not yet been set (private$m is not NULL).")
				}
				assertIntegerish(m, lower = 0, any.missing = FALSE, min.len = 1)
			}
			private$m = as.integer(m)
			invisible(self)
		},
		#' @description If the design is a block design, get block identifiers (otherwise halts)
		#'
		#' @return An integer vector of block identifiers.
		get_block_ids = function(){
			if (!is.null(private$m) && length(private$m) == length(private$y)) {
				return(private$m)
			}
			block_ids = private$m
			Xraw = private$Xraw
			strata_cols = private$strata_cols
			if (is.null(block_ids) && is(self, "DesignFixedOptimalBlocks")) {
				block_ids = private$get_or_compute_block_ids()
			}
			if (is.null(block_ids) &&
					(is(self, "DesignFixedBlocking") ||
					 is(self, "DesignFixedBlockedCluster") ||
					 (!is.null(strata_cols) && length(strata_cols) > 0L)) &&
					nrow(Xraw) == length(private$y)) {
				strata_keys = private$get_strata_keys()
				if (length(strata_keys) == length(private$y)) {
					block_ids = match(strata_keys, unique(strata_keys))
				}
			}
			if (is.null(block_ids)) {
				stop("Block identifiers are undefined for this design.")
			}
			block_ids = as.integer(block_ids)
			if (length(block_ids) != length(private$y)) {
				stop("Block identifiers are improperly sized for this design.")
			}
			private$m = block_ids
			block_ids
		},
		#' @description Summarize the covariates within blocks.
		#'
		#' @param block_ids A vector of block identifiers to summarize. If NULL, defaults to all blocks.
		#'
		#' @return A list of data.tables, one for each block.
		summarize_blocks = function(block_ids = NULL) {
			if (should_run_asserts()) {
				self$assert_blocking_design()
			}
			all_block_ids = self$get_block_ids()
			if (is.null(block_ids)) {
				block_ids = sort(unique(all_block_ids))
			}
			Xraw = self$get_X_raw()
			res = list()
			for (bid in block_ids) {
				idx = which(all_block_ids == bid)
				if (length(idx) == 0) {
					warning(paste("Block ID", bid, "not found in the design."))
					next
				}
				X_block = Xraw[idx, ]
				col_summaries = list()
				for (col_name in names(X_block)) {
					col_data = X_block[[col_name]]
					if (is.numeric(col_data)) {
						col_summaries[[col_name]] = data.table::data.table(
							variable = col_name,
							level = NA_character_,
							mean_or_pct = mean(col_data, na.rm = TRUE),
							sd = stats::sd(col_data, na.rm = TRUE)
						)
					} else {
						counts = table(col_data, useNA = "no")
						pcts = prop.table(counts) * 100
						if (length(pcts) == 0) {
							col_summaries[[col_name]] = data.table::data.table(
								variable = col_name,
								level = "N/A",
								mean_or_pct = NA_real_,
								sd = NA_real_
							)
						} else {
							col_summaries[[col_name]] = data.table::data.table(
								variable = col_name,
								level = names(pcts),
								mean_or_pct = as.numeric(pcts),
								sd = NA_real_
							)
						}
					}
				}
				res[[as.character(bid)]] = data.table::rbindlist(col_summaries)
			}
			res
		}
	),
	private = list(
		assert_min_block_size = function(n, B) {
			if (should_run_asserts() && floor(n / B) < 2L) {
				stop(
					"Cannot use B = ", B, " with n = ", n, ": ",
					"floor(n / B) = ", floor(n / B), " < 2. Minimum block size is 2."
				)
			}
			invisible(NULL)
		},
		m = NULL,
		strata_cols = NULL,
		preferred_num_bins_for_continuous_covariate = NULL,
		B_target = NULL,
		exact_num_blocks = FALSE,
		equal_block_sizes = TRUE,
		blocking_capable = FALSE,
		cmh_se_w_mat = NULL,
		get_strata_keys = function(){
			n = private$t
			if (n == 0) return(character(0))
			strata_cols = if (is.null(private$strata_cols)) names(private$Xraw) else private$strata_cols
			target = if (!is.null(private$B_target)) {
				if (is.na(private$B_target)) floor(sqrt(n)) else as.integer(private$B_target)
			} else {
				NULL
			}
			col_to_str = function(col, num_bins = private$preferred_num_bins_for_continuous_covariate) {
				vec = private$Xraw[[col]]
				if (is.numeric(vec)) {
					probs = seq(0, 1, length.out = num_bins + 1)
					breaks = unique(stats::quantile(vec, probs = probs, na.rm = TRUE))
					s = if (length(breaks) > 1) as.character(cut(vec, breaks = breaks, include.lowest = TRUE)) else as.character(vec)
				} else {
					s = as.character(vec)
				}
				s[is.na(s)] = "NA"
				s
			}
			append_key = function(keys, col_str) {
				if (all(nchar(keys) == 0L)) col_str else paste(keys, col_str, sep = "|")
			}
			has_equal_block_sizes = function(keys) {
				block_counts = as.integer(table(keys))
				length(block_counts) <= 1L || all(block_counts == block_counts[1L])
			}
			choose_column_keys = function(keys, col) {
				vec = private$Xraw[[col]]
				if (!isTRUE(private$exact_num_blocks) || is.null(target) || !is.numeric(vec)) {
					return(append_key(keys, col_to_str(col)))
				}
				max_bins = max(1L, min(n, target))
				candidate_bins = unique(c(
					as.integer(private$preferred_num_bins_for_continuous_covariate),
					seq_len(max_bins)
				))
				best_keys = keys
				best_num_blocks = length(unique(keys))
				for (num_bins in candidate_bins) {
					candidate_keys = append_key(keys, col_to_str(col, num_bins = num_bins))
					num_blocks = length(unique(candidate_keys))
					if (num_blocks > target) next
					if (isTRUE(private$equal_block_sizes) && !has_equal_block_sizes(candidate_keys)) next
					if (num_blocks > best_num_blocks) {
						best_keys = candidate_keys
						best_num_blocks = num_blocks
					}
					if (num_blocks == target) break
				}
				best_keys
			}
			keys = rep("", n)
			if (!is.null(target)) {
				for (col in strata_cols) {
					new_keys = choose_column_keys(keys, col)
					if (length(unique(new_keys)) <= target) {
						keys = new_keys
					}
				}
			} else {
				for (col in strata_cols) {
					keys = append_key(keys, col_to_str(col))
				}
			}
			num_blocks = length(unique(keys))
			if (isTRUE(private$equal_block_sizes)) {
				block_counts = as.integer(table(keys))
				if (length(block_counts) > 1L && any(block_counts != block_counts[1L])) {
					stop("equal_block_sizes = TRUE but the strata produced unequal block sizes (",
						paste(sort(unique(block_counts)), collapse = ", "),
						"). Set equal_block_sizes = FALSE or adjust your covariate binning.")
				}
			}
			if (isTRUE(private$exact_num_blocks)) {
				if (is.null(private$B_target)) {
					stop("exact_num_blocks requires B_target.")
				}
				target = if (is.na(private$B_target)) floor(sqrt(n)) else as.integer(private$B_target)
				if (num_blocks != target) {
					stop("exact_num_blocks = TRUE but the greedy blocking key construction produced ", num_blocks,
						" blocks instead of the requested ", target, ".")
				}
			}
			if (should_run_asserts()) {
				if (num_blocks > n) {
					stop("Number of blocks (", num_blocks, ") exceeds sample size (", n, "). Reduce the number of strata columns or use fewer bins.")
				}
			}
			keys
		}
	)
)
