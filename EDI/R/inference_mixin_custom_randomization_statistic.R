#' Mixin for User-Supplied Custom Randomization Statistics
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private}
#' slots) bundling the machinery that lets a user override the default
#' randomization-inference statistic with either a plain R function or a
#' compiled/compilable C++ function. Splice into a daughter class (in
#' practice, just \code{InferenceRand} itself, once) via
#' \code{public = c(InferenceMixinCustomRandomizationStatistic$public, list(...))}
#' and
#' \code{private = c(InferenceMixinCustomRandomizationStatistic$private, list(...))}.
#'
#' Self-contained: the only host-class state it touches beyond its own three
#' private fields is \code{private$cached_values}, the shared invalidation
#' cache already present on every \code{Inference} subclass.
#'
#' @keywords internal
#' @noRd
InferenceMixinCustomRandomizationStatistic = list(
	public = list(
		#' @description Set Custom Randomization Statistic Computation
		#' @param custom_randomization_statistic_function  A function that returns a scalar value.
		set_custom_randomization_statistic_function = function(custom_randomization_statistic_function){
			if (!is.null(custom_randomization_statistic_function) && !is.null(private[["compiled_cpp_stat_fn"]])) {
				stop("Cannot specify both custom_randomization_statistic_function and custom_randomization_statistic_cpp.")
			}
			if (should_run_asserts()) {
				assertFunction(custom_randomization_statistic_function, null.ok = TRUE)
			}
			private[["custom_randomization_statistic_function"]] = custom_randomization_statistic_function
			private$cached_values$t0s_rand = NULL
			private$cached_values$rand_distr_cache = list()
			private$cached_values$custom_stat_analysis = NULL
		},
		#' @description Set Custom Randomization Statistic as a C++ Function
		#' @param fn  Either a C++ source code string or a pre-compiled Rcpp function returning a
		#'   scalar \code{double}. Must accept either \code{(NumericVector y, IntegerVector w)} or
		#'   \code{(NumericVector y, IntegerVector w, IntegerVector dead)} as arguments.
		#'   Passing a source string (recommended) enables safe use with parallel workers; passing
		#'   a pre-compiled function works only in the main process (XPtrs become NULL when
		#'   serialized to worker nodes). Pass \code{NULL} to clear.
		#'   Cannot be combined with \code{set_custom_randomization_statistic_function}.
		set_custom_randomization_statistic_cpp = function(fn){
			if (!is.null(fn) && !is.null(private[["custom_randomization_statistic_function"]])) {
				stop("Cannot specify both custom_randomization_statistic_function and custom_randomization_statistic_cpp.")
			}
			if (!is.null(fn)) {
				if (is.character(fn) && length(fn) == 1L) {
					compiled = Rcpp::cppFunction(fn)
					arity = length(formals(compiled))
					if (!arity %in% c(2L, 3L)) stop("custom_randomization_statistic_cpp source must define a function with 2 arguments (y, w) or 3 arguments (y, w, dead); got ", arity, ".")
					private[["compiled_cpp_stat_src"]] = fn
					private[["compiled_cpp_stat_fn"]] = compiled
				} else {
					if (!is.function(fn)) stop("custom_randomization_statistic_cpp must be a C++ source string or compiled Rcpp function, not a ", class(fn)[1], ".")
					arity = length(formals(fn))
					if (!arity %in% c(2L, 3L)) stop("custom_randomization_statistic_cpp must accept 2 arguments (y, w) or 3 arguments (y, w, dead); got ", arity, ".")
					private[["compiled_cpp_stat_src"]] = NULL
					private[["compiled_cpp_stat_fn"]] = fn
				}
			} else {
				private[["compiled_cpp_stat_src"]] = NULL
				private[["compiled_cpp_stat_fn"]] = NULL
			}
			private$cached_values$t0s_rand = NULL
			private$cached_values$rand_distr_cache = list()
			private$cached_values$custom_stat_analysis = NULL
		}
	),
	private = list(
		custom_randomization_statistic_function = NULL,
		compiled_cpp_stat_fn = NULL,
		compiled_cpp_stat_src = NULL,
		get_compiled_cpp_stat = function() private[["compiled_cpp_stat_fn"]],
		analyze_custom_randomization_statistic = function(){
			if (!is.null(private$cached_values$custom_stat_analysis)) return(private$cached_values$custom_stat_analysis)
			if (is.null(private$custom_randomization_statistic_function) && is.null(private[["compiled_cpp_stat_fn"]])) {
				analysis = list(can_use_lightweight_yw_only = FALSE, needs_match_data = TRUE)
				private$cached_values$custom_stat_analysis = analysis; return(analysis)
			}
			# C++ stat is always lightweight: it only ever receives (y, w) or (y, w, dead).
			if (!is.null(private[["compiled_cpp_stat_fn"]])) {
				cpp_src = private[["compiled_cpp_stat_src"]]
				# Build a self-compiling closure that is safe to serialize to parallel workers.
				# It captures only the source string (not the XPtr), so workers compile their
				# own local copy on first use rather than receiving an invalid null pointer.
				get_cpp_fn = if (!is.null(cpp_src)) {
					local({
						.src = cpp_src
						.fn  = NULL
						function() {
							if (is.null(.fn)) .fn <<- Rcpp::cppFunction(.src)
							.fn
						}
					})
				} else {
					NULL
				}
				analysis = list(can_use_lightweight_yw_only = TRUE, needs_match_data = FALSE, get_cpp_fn = get_cpp_fn)
				private$cached_values$custom_stat_analysis = analysis; return(analysis)
			}
			# Basic analysis: does it only use y and w?
			body_str = paste(deparse(body(private$custom_randomization_statistic_function)), collapse = " ")
			# Look for access to other members of private$des_obj_priv_int
			can_use_lightweight = !grepl("private\\$des_obj_priv_int\\$(?!y|w|dead)", body_str, perl = TRUE)
			analysis = list(can_use_lightweight_yw_only = can_use_lightweight, needs_match_data = FALSE)
			private$cached_values$custom_stat_analysis = analysis
			analysis
		},
		evaluate_lightweight_custom_randomization_statistic = function(lightweight_custom_context, y, w, dead, cpp_fn_override = NULL){
			# Fast path: compiled C++ function — no R interpreter overhead per permutation.
			# cpp_fn_override is used by parallel workers (compiled fresh from source on each
			# worker) so they never dereference the XPtr that was serialized from the main process.
			cpp_fn = if (!is.null(cpp_fn_override)) cpp_fn_override else private$get_compiled_cpp_stat()
			if (!is.null(cpp_fn)) {
				arity = length(formals(cpp_fn))
				return(as.numeric(
					if (arity >= 3L) cpp_fn(y, as.integer(w), as.integer(dead))
					else cpp_fn(y, as.integer(w))
				)[1L])
			}
			# We simulate the environment for the custom statistic
			fn = private$custom_randomization_statistic_function
			old_env = environment(fn)
			on.exit(environment(fn) <- old_env, add = TRUE)
			eval_env = new.env(parent = environment(fn))

			private_proxy = new.env(parent = emptyenv())
			seq_priv_proxy = new.env(parent = emptyenv())
			seq_priv_proxy$y = y; seq_priv_proxy$w = w; seq_priv_proxy$dead = dead
			private_proxy$des_obj_priv_int = seq_priv_proxy

			eval_env$private = private_proxy
			eval_env$inf_priv = private_proxy
			eval_env$des_priv = seq_priv_proxy
			eval_env$des_obj_priv_int = seq_priv_proxy

			environment(fn) = eval_env
			eval_env$.custom_randomization_statistic_function = fn
			eval(quote(.custom_randomization_statistic_function()), envir = eval_env)
		}
	)
)
