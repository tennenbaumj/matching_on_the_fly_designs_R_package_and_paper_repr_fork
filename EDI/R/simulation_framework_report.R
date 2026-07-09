#' Reporting class for SimulationFramework results
#'
#' @description An R6 class for accessing and summarizing the results of a
#' \code{SimulationFramework} run.  It can be constructed either from a completed
#' \code{SimulationFramework} object (via \code{SimulationFrameworkReport$new(sim)})
#' or by loading results from a previously saved CSV / CSV.BZ2
#' file (via \code{SimulationFrameworkReport$new("path/to/results.csv")}).
#'
#' @details When constructed from a \code{SimulationFramework} object all
#' design/inference parameter metadata is preserved, so \code{$summarize()} can
#' annotate each row with human-readable parameter strings.  When constructed from
#' a file only the raw results are available; parameter annotation columns will be
#' empty strings.
#'
#' @examples
#' \donttest{
#' sim <- SimulationFramework$new(
#'   response_type = "continuous",
#'   design_classes_and_params = list(DesignFixedBernoulli),
#'   inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
#'   n = 20L, Nrep = 5L, betaT = 1,
#'   results_filename = tempfile(fileext = ".csv"),
#'   verbose = FALSE, continue_from_last_result_row = FALSE
#' )
#' sim$run()
#' report <- SimulationFrameworkReport$new(sim)
#' report$get_results()
#' report$summarize()
#' }
#' @export
SimulationFrameworkReport = R6::R6Class("SimulationFrameworkReport",
  lock_objects = FALSE,
  public = list(
    #' @description Create a new \code{SimulationFrameworkReport}.
    #'
    #' @param sim_or_filename Either a completed \code{SimulationFramework} object
    #'   or a character string giving the path to a \code{.csv} or \code{.csv.bz2}
    #'   results file written by \code{SimulationFramework}.
    #' @param alpha Numeric in \eqn{(0,1)}.  Significance level for coverage and
    #'   power calculations.  When \code{sim_or_filename} is a
    #'   \code{SimulationFramework} object and \code{alpha} is \code{NULL}
    #'   (default), the framework's own alpha is used.  When loading from a file,
    #'   defaults to \code{0.05}.
    initialize = function(sim_or_filename, alpha = NULL) {
      if (is.character(sim_or_filename)) {
        private$.init_from_file(sim_or_filename, alpha %||% 0.05)
      } else if (inherits(sim_or_filename, "SimulationFramework")) {
        private$.init_from_sim(sim_or_filename, alpha)
      } else {
        stop("'sim_or_filename' must be a SimulationFramework object or a character filename")
      }
    },
    #' @description Get the raw per-replication results.
    #'
    #' @return A \code{data.table} with one row per
    #'   (replication, design, inference class, inference type).
    get_results = function() {
      private$results_dt
    },
    #' @description Return all errors captured during the simulation run.
    #'
    #' @return A list of named lists, one per captured error.  Each element
    #'   includes the simulation cell metadata, replication number,
    #'   design / inference path, user-supplied parameters, error stage, and
    #'   error message.  Empty when constructed from a file.
    get_errors = function() {
      private$errors
    },
    #' @description Aggregate and summarize simulation results.
    #'
    #' @return A \code{data.table} with one row per unique
    #'   (response_type, cond_exp_func_model, n, p, betaT, design, inference,
    #'   inference_type) combination.  Columns include \code{MSE},
    #'   \code{coverage}, \code{ci_length} (when CI types were run),
    #'   \code{power} (when betaT != 0 and p-value types were run),
    #'   \code{size} and \code{size_pval} (when betaT == 0 and p-value types
    #'   were run; \code{size_pval} is the exact two-sided binomial test
    #'   p-value of H0: true size = alpha, suitable for multiplicity-corrected
    #'   calibration checks across settings), and parameter annotation strings.
    summarize = function() {
      if (length(private$valid_combos) == 0L) {
        message("No results."); return(invisible(NULL))
      }
      # ── Reference grid ───────────────────────────────────────────────────────
      ref_grid = data.table::rbindlist(lapply(private$valid_combos, as.list), use.names = TRUE, fill = TRUE)
      ref_grid = unique(ref_grid[, setdiff(names(ref_grid), "rep"), with = FALSE])
      # ── Per-class params strings ──────────────────────────────────────────────
      if (!is.null(private$design_labels)) {
        design_params_map = stats::setNames(
          lapply(seq_along(private$design_labels), function(di)
            private$.params_to_str(if (!is.null(private$design_params)) private$design_params[[di]] else NULL)),
          private$design_labels)
      } else {
        unique_designs = unique(private$results_dt$design)
        design_params_map = stats::setNames(rep("", length(unique_designs)), unique_designs)
      }
      if (!is.null(private$inference_labels)) {
        inf_params_map = stats::setNames(
          lapply(seq_along(private$inference_labels), function(ii)
            private$.params_to_str(private$inference_constructor_params[[ii]])),
          private$inference_labels)
      } else {
        unique_infs = unique(private$results_dt$inference)
        inf_params_map = stats::setNames(rep("", length(unique_infs)), unique_infs)
      }
      ref_grid[, design_params    := unlist(design_params_map[design])]
      ref_grid[, inference_params := unlist(inf_params_map[inference])]
      ref_grid[, inference_type_params := vapply(inference_type, private$.params_for_inference_type_to_str, "")]
      # ── Aggregate raw results ─────────────────────────────────────────────────
      dt         = private$results_dt
      alpha      = private$alpha
      report_cov = any(grepl("_ci$",   private$inf_types))
      report_pow = any(grepl("_pval$", private$inf_types))
      has_sim_mode = "simulation_mode" %in% names(dt)
      by_cols = c("response_type", "cond_exp_func_model", "n", "p", "betaT", "design", "inference", "inference_type",
                  if (has_sim_mode) "simulation_mode")
      if (nrow(dt) > 0L) {
        agg = dt[, {
          est_fin = is.finite(estimate)
          m_row = list(
            MSE   = if (any(est_fin)) mean((estimate[est_fin] - true_estimand[est_fin])^2) else NA_real_,
            n_est = sum(est_fin)
          )
          if (report_cov) {
            ci_fin = is.finite(ci_lo) & is.finite(ci_hi)
            m_row$coverage  = if (any(ci_fin)) mean(ci_lo[ci_fin] <= true_estimand[ci_fin] & true_estimand[ci_fin] <= ci_hi[ci_fin]) else NA_real_
            m_row$n_cov     = sum(ci_fin)
            m_row$ci_length = if (any(ci_fin)) mean(ci_hi[ci_fin] - ci_lo[ci_fin]) else NA_real_
          }
          if (report_pow) {
            pv_fin = is.finite(pval)
            # For custom DGPs, use mean(true_estimand) to classify null vs alternative;
            # fall back to betaT for the standard DGP path.
            uses_custom_estimand = has_sim_mode &&
              length(simulation_mode) > 0L && simulation_mode[1L] != "standard"
            is_zero = if (uses_custom_estimand) {
              abs(mean(true_estimand[is.finite(true_estimand)])) < 1e-12
            } else {
              abs(betaT) < 1e-12
            }
            pv_fin_zero    = pv_fin & is_zero
            pv_fin_nonzero = pv_fin & !is_zero
            m_row$power  = if (any(pv_fin_nonzero)) mean(pval[pv_fin_nonzero] < alpha) else NA_real_
            m_row$n_pow  = sum(pv_fin_nonzero)
            m_row$size   = if (any(pv_fin_zero)) mean(pval[pv_fin_zero] < alpha) else NA_real_
            m_row$n_size = sum(pv_fin_zero)
            # Exact two-sided binomial test of H0: true size = alpha, so users
            # can flag miscalibrated tests across settings (e.g. with Bonferroni)
            m_row$size_pval = if (any(pv_fin_zero)) {
              stats::binom.test(sum(pval[pv_fin_zero] < alpha), sum(pv_fin_zero), p = alpha)$p.value
            } else NA_real_
          }
          m_row
        }, by = by_cols]
      } else {
        agg = data.table::data.table(
          response_type = character(), cond_exp_func_model = character(),
          n = integer(), p = integer(), betaT = numeric(),
          design = character(), inference = character(),
          inference_type = character(), MSE = numeric(),
          n_est = integer(), power = numeric(), n_pow = integer()
        )
      }
      # ── Right-join: every valid combo appears, NA for those with no data ──────
      data.table::setkeyv(agg,      by_cols)
      data.table::setkeyv(ref_grid, by_cols)
      result = agg[ref_grid]
      if (!"n_est"  %in% names(result)) result[, n_est  := 0L]
      if (!"n_pow"  %in% names(result)) result[, n_pow  := 0L]
      if (!"n_size" %in% names(result)) result[, n_size := 0L]
      result[is.na(n_est),  n_est  := 0L]
      result[is.na(n_pow),  n_pow  := 0L]
      result[is.na(n_size), n_size := 0L]
      result[order(cond_exp_func_model, n, p, betaT, design, inference, inference_type)]
    },
    #' @description Print a concise summary of the report.
    print = function() {
      cat("SimulationFrameworkReport\n")
      sm = self$summarize()
      if (!is.null(sm) && nrow(sm) > 0L) {
        cat(sprintf("Summary (alpha = %g):\n", private$alpha))
        print(sm)
      } else {
        cat("  No results.\n")
      }
      invisible(self)
    }
  ),
  # ── private ──────────────────────────────────────────────────────────────────
  private = list(
    results_dt                   = NULL,
    errors                       = NULL,
    valid_combos                 = NULL,
    alpha                        = NULL,
    inf_types                    = NULL,
    design_labels                = NULL,
    design_params                = NULL,
    inference_labels             = NULL,
    inference_constructor_params = NULL,
    inference_type_params        = NULL,

    .init_from_sim = function(sim, alpha_override) {
      p = sim$.__enclos_env__$private
      if (!p$has_run) stop("SimulationFramework must be run first (call $run()).")
      private$results_dt = if (p$results_idx == 0L) {
        data.table::data.table(
          rep = integer(), response_type = character(),
          cond_exp_func_model = character(), n = integer(),
          p = integer(), betaT = numeric(), design = character(),
          inference = character(), inference_type = character(),
          estimate = numeric(), ci_lo = numeric(), ci_hi = numeric(),
          pval = numeric(), true_estimand = numeric()
        )
      } else {
        data.table::rbindlist(p$raw_results[seq_len(p$results_idx)], use.names = TRUE, fill = TRUE)
      }
      private$errors                       = p$error_log
      private$valid_combos                 = p$valid_combos
      private$alpha                        = if (!is.null(alpha_override)) alpha_override else p$alpha
      private$inf_types                    = p$inf_types
      private$design_labels                = p$design_labels
      private$design_params                = p$design_params
      private$inference_labels             = p$inference_labels
      private$inference_constructor_params = p$inference_constructor_params
      private$inference_type_params        = p$inference_type_params
    },

    .init_from_file = function(filename, alpha) {
      if (!file.exists(filename)) stop("File not found: ", filename)
      fmt = if (grepl("\\.csv\\.bz2$", filename, ignore.case = TRUE)) "csv.bz2"
            else if (grepl("\\.csv$", filename, ignore.case = TRUE)) "csv"
            else stop("Unsupported file format; must end in .csv or .csv.bz2")
      if (identical(fmt, "csv")) {
        dt = data.table::fread(filename, showProgress = FALSE)
      } else {
        extracted = tempfile("sfr_", fileext = ".csv")
        input_con  = bzfile(filename, open = "rb")
        output_con = file(extracted, open = "wb")
        tryCatch({
          repeat {
            bytes = readBin(input_con, what = "raw", n = 1024L * 1024L)
            if (length(bytes) == 0L) break
            writeBin(bytes, output_con)
          }
        }, finally = {
          try(close(input_con),  silent = TRUE)
          try(close(output_con), silent = TRUE)
        })
        dt = data.table::fread(extracted, showProgress = FALSE)
        unlink(extracted, force = TRUE)
      }
      for (col in c("response_type", "cond_exp_func_model", "design", "inference", "inference_type"))
        if (col %in% names(dt)) dt[, (col) := as.character(get(col))]
      for (col in c("rep", "n", "p"))
        if (col %in% names(dt)) dt[, (col) := as.integer(get(col))]
      for (col in c("betaT", "estimate", "ci_lo", "ci_hi", "pval", "true_estimand"))
        if (col %in% names(dt)) dt[, (col) := as.numeric(get(col))]
      private$results_dt = dt
      private$errors     = list()
      private$alpha      = alpha
      if (nrow(dt) > 0L) {
        combo_cols = intersect(
          c("response_type", "cond_exp_func_model", "n", "p", "betaT",
            "design", "inference", "inference_type", "simulation_mode"),
          names(dt)
        )
        combos_dt = unique(dt[, combo_cols, with = FALSE])
        private$valid_combos = lapply(seq_len(nrow(combos_dt)), function(i) as.list(combos_dt[i]))
        private$inf_types = if ("inference_type" %in% names(dt)) unique(dt$inference_type) else character(0L)
      } else {
        private$valid_combos = list()
        private$inf_types    = character(0L)
      }
      private$design_labels                = NULL
      private$design_params                = NULL
      private$inference_labels             = NULL
      private$inference_constructor_params = NULL
      private$inference_type_params        = NULL
    },

    .params_to_str = function(p) {
      if (is.null(p) || length(p) == 0L) return("")
      kv = mapply(function(k, v) paste0(k, "=", paste(deparse(v), collapse = "")),
                  names(p), p, SIMPLIFY = TRUE)
      paste(kv, collapse = ", ")
    },

    .params_for_inference_type_to_str = function(inference_type) {
      if (is.null(private$inference_type_params)) return("")
      ps = private$.params_to_str(private$inference_type_params[[inference_type]])
      if (nchar(ps) > 0L) paste0(inference_type, "(", ps, ")") else ""
    }
  )
)
