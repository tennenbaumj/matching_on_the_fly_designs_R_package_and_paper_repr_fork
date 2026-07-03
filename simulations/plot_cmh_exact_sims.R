rm(list = ls())
pacman::p_load(data.table, R.utils, ggplot2, scales, gridExtra, xtable)

Nrep = 10002
raw_results_dt = fread(sprintf("cmh_exact_sims_plus_greedy_results_high_signal_Nrep_%d.csv.bz2", Nrep))
raw_results_dt[, reject := pval < 0.05]
raw_results_dt[, covers := ci_lo <= true_estimand & true_estimand <= ci_hi]
raw_results_dt[, ci_length := ci_hi - ci_lo]
table(raw_results_dt$design)
results_dt = raw_results_dt[,
  .(
    pow_avg = mean(reject, na.rm = TRUE), 
    n_pow = sum(!is.na(reject)),
    cov_avg = mean(covers, na.rm = TRUE), 
    n_cov = sum(!is.na(covers)),
    len_avg = mean(ci_length, na.rm = TRUE), 
    len_sd = sd(ci_length, na.rm = TRUE), 
    n_len = sum(!is.na(ci_length)),
    true_estimand = first(true_estimand)
  ),
  by = c("n", "p", "betaT", "design", "inference", "cond_exp_func_model")                         
]
estimands = results_dt[,
  .(true_estimand = first(true_estimand)),
  by = c("n", "p", "betaT", "cond_exp_func_model")
]
#table(results_dt$design, results_dt$inference)
#table(results_dt$n, results_dt$p, results_dt$betaT)

# Wald CIs for power and coverage
z = qnorm(0.975)
results_dt[, `:=`(
  pow_a = pow_avg - z * sqrt(pow_avg * (1 - pow_avg) / n_pow),
  pow_b = pow_avg + z * sqrt(pow_avg * (1 - pow_avg) / n_pow),
  cov_a = cov_avg - z * sqrt(cov_avg * (1 - cov_avg) / n_cov),
  cov_b = cov_avg + z * sqrt(cov_avg * (1 - cov_avg) / n_cov),
  len_a = len_avg - z * len_sd / sqrt(n_len),
  len_b = len_avg + z * len_sd / sqrt(n_len)
)]

# Shorten labels
results_dt[, design_short    := gsub("DesignFixed", "", design)]
table(results_dt$design)
table(results_dt$design_short)
results_dt[, inference_short := inference]
results_dt[, inference_short := gsub("InferenceIncid", "", inference_short)]

# Row 1: iBCRD, OptimalBlocks B=4/8/16/32; Row 2: BinaryMatch, Blocking B_target=4/8/16/32
design_levels = c(
  "iBCRD",
  "OptimalBlocks (B=2)",
  "OptimalBlocks (B=4)",
  "OptimalBlocks (B=8)",
  "OptimalBlocks (B=16)",
  "OptimalBlocks (B=32)",
  "OptimalBlocks (B=64)",
  "OptimalBlocks (B=128)",
  "Blocking (B_target=2, exact_num_blocks=TRUE)",
  "Blocking (B_target=4, exact_num_blocks=TRUE)",
  "Blocking (B_target=8, exact_num_blocks=TRUE)",
  "Blocking (B_target=16, exact_num_blocks=TRUE)",
  "Blocking (B_target=32, exact_num_blocks=TRUE)",
  "Blocking (B_target=64, exact_num_blocks=TRUE)",
  "Blocking (B_target=128, exact_num_blocks=TRUE)",
  "BinaryMatch",
  "Rerandomization (prop_acceptable=0.01)",
  'Greedy (objective=""abs_sum_diff"")',
  'MatchingGreedyPairSwitching (objective=""abs_sum_diff"")',
  'Greedy (objective=""mahal_dist"")',
  'MatchingGreedyPairSwitching (objective=""mahal_dist"")'
)
design_labels = c(
  "BCRD",
  "Optimal B=2",
  "Optimal B=4",
  "Optimal B=8",
  "Optimal B=16",
  "Optimal B=32",
  "Optimal B=64",
  "Optimal B=128",
  "Naive B=2",
  "Naive B=4",
  "Naive B=8",
  "Naive B=16",
  "Naive B=32",
  "Naive B=64",
  "Naive B=128",
  "BinaryMatch",
  "Rerandomization",
  "GreedyAbs",
  "BinaryMatchThenGreedyAbs",
  "GreedyMD",
  "BinaryMatchThenGreedyMD"
)

results_dt[, design_short := factor(design_short, levels = design_levels, labels = design_labels)]
results_dt = results_dt[!(design_short %in% c("GreedyAbs", "BinaryMatchThenGreedyAbs"))]
table(results_dt$design_short)
#results_dt = results_dt[design_short %in% c("Rerandomization", "Greedy", "BinaryMatchThenGreedy")]

extract_legend = function(plot_obj) {
  plot_grob = ggplotGrob(plot_obj)
  guide_idx = which(vapply(plot_grob$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(guide_idx) == 0L) {
    return(NULL)
  }
  plot_grob$grobs[[guide_idx[1L]]]
}

build_metric_plot = function(dat, metric, xlab, hline = NULL, show_y_axis = TRUE, show_legend = FALSE, hline_linetype = "dashed", show_n_strip = TRUE, free_metric_axis = FALSE, x_label_digits = NULL, x_breaks = NULL) {
  dat = copy(dat)
  dat[, y := get(metric)]

  has_ci_cols = paste0(substr(metric, 1, 3), "_a") %in% names(dat)
  if (has_ci_cols) {
    dat[, y_lo := get(paste0(substr(metric, 1, 3), "_a"))]
    dat[, y_hi := get(paste0(substr(metric, 1, 3), "_b"))]
  }

  dodge = position_dodge(width = 0.4)
  dat[, n_facet := factor(sprintf("n = %s", n), levels = sprintf("n = %s", sort(unique(n))))]
  plot_obj = ggplot(dat, aes(x = design_short, y = y,
                      color = inference_short, group = inference_short)) +
    geom_point(size = 1.0, position = dodge) +
    facet_grid(rows = vars(n_facet), scales = if (free_metric_axis) "free_y" else "fixed") +
    scale_color_brewer(palette = "Set1") +
    coord_flip() +
    labs(x = NULL, y = xlab, color = "Inference") +
    theme_bw(base_size = 9) +
    theme(
      strip.text.y     = if (show_n_strip) element_text(size = 10) else element_blank(),
      strip.background = if (show_n_strip) element_rect() else element_blank(),
      legend.position  = if (show_legend) "bottom" else "none",
      legend.text      = element_text(size = 12),
      legend.title     = element_text(size = 12),
      axis.title       = element_text(size = 13),
      axis.text        = element_text(size = 11),
      axis.title.x     = element_text(size = 13),
      axis.title.y     = if (show_y_axis) element_blank() else element_blank(),
      axis.text.y      = if (show_y_axis) element_text(size = 8) else element_blank(),
      axis.ticks.y     = if (show_y_axis) element_line() else element_blank(),
      panel.grid.minor = element_blank()
    )
  if (!is.null(hline)) {
    plot_obj = plot_obj + geom_hline(yintercept = hline, linetype = hline_linetype, color = "gray60", linewidth = 0.4)
  }
  if (has_ci_cols) {
    plot_obj = plot_obj + geom_errorbar(aes(ymin = y_lo, ymax = y_hi), width = 0.2, linewidth = 0.4, position = dodge)
  }
  if (!is.null(x_label_digits) || !is.null(x_breaks)) {
    acc = if (!is.null(x_label_digits)) 10^(-x_label_digits) else 0.01
    lbl_fn = function(x) {
      lbl = scales::label_number(accuracy = acc)(x)
      ifelse(abs(x - 1) < 1e-9, "", lbl)
    }
    plot_obj = plot_obj + scale_y_continuous(
      breaks = if (!is.null(x_breaks)) x_breaks else waiver(),
      labels = lbl_fn
    )
  }
  plot_obj
}

make_two_column_plot = function(
  dat_left, metric_left, xlab_left, hline_left,
  dat_right, metric_right, xlab_right, hline_right,
  title_str, filename_stem, plot = TRUE, save_PDF = FALSE, save_PNG = FALSE,
  hline_linetype_left = "dashed", hline_linetype_right = "dashed",
  free_metric_axis_left = FALSE, free_metric_axis_right = FALSE,
  x_label_digits_left = NULL, x_label_digits_right = NULL,
  x_breaks_left = NULL, x_breaks_right = NULL
) {
  left_plot_for_legend = build_metric_plot(
    dat_left, metric_left, xlab_left, hline = hline_left,
    show_y_axis = TRUE, show_legend = TRUE, hline_linetype = hline_linetype_left,
    show_n_strip = FALSE, free_metric_axis = free_metric_axis_left,
    x_label_digits = x_label_digits_left, x_breaks = x_breaks_left
  )
  legend_grob = extract_legend(left_plot_for_legend)

  left_plot = build_metric_plot(
    dat_left, metric_left, xlab_left, hline = hline_left,
    show_y_axis = TRUE, show_legend = FALSE, hline_linetype = hline_linetype_left,
    show_n_strip = FALSE, free_metric_axis = free_metric_axis_left,
    x_label_digits = x_label_digits_left, x_breaks = x_breaks_left
  )
  right_plot = build_metric_plot(
    dat_right, metric_right, xlab_right, hline = hline_right,
    show_y_axis = FALSE, show_legend = FALSE, hline_linetype = hline_linetype_right,
    show_n_strip = TRUE, free_metric_axis = free_metric_axis_right,
    x_label_digits = x_label_digits_right, x_breaks = x_breaks_right
  )

  combined_row = gridExtra::arrangeGrob(left_plot, right_plot, ncol = 2)
  if (is.null(legend_grob)) {
    combined_plot = gridExtra::arrangeGrob(combined_row, top = title_str)
  } else {
    combined_plot = gridExtra::arrangeGrob(
      combined_row,
      legend_grob,
      ncol = 1,
      # top = title_str,
      heights = c(10, 1.2)
    )
  }

  if (plot) {
    grid::grid.newpage()
    grid::grid.draw(combined_plot)
  }
  if (save_PDF) {
    ggsave(sprintf("%s.pdf", filename_stem), combined_plot, width = 6.5, height = 8)
  }
  if (save_PNG) {
    ggsave(sprintf("%s.png", filename_stem), combined_plot, width = 6.5, height = 8, dpi = 150)
  }
  invisible(combined_plot)
}

make_original_plot = function(dat, metric, ylab, hline, title_str, filename_stem, plot = TRUE, save_PDF = FALSE, save_PNG = FALSE, free_y = FALSE, hline_linetype = "dashed") {
  dat = copy(dat)
  dat[, y := get(metric)]

  has_ci_cols = paste0(substr(metric, 1, 3), "_a") %in% names(dat)
  if (has_ci_cols) {
    dat[, y_lo := get(paste0(substr(metric, 1, 3), "_a"))]
    dat[, y_hi := get(paste0(substr(metric, 1, 3), "_b"))]
  }

  dodge = position_dodge(width = 0.4)
  plot_obj = ggplot(dat, aes(x = design_short, y = y,
                      color = inference_short, group = inference_short)) +
    geom_point(size = 0.5, position = dodge) +
    facet_wrap(~ factor(n), nrow = 2, ncol = 5, scales = if (free_y) "free_y" else "fixed") +
    scale_color_brewer(palette = "Set1") +
    labs(x = "design", y = ylab, color = "Inference", title = title_str) +
    theme_bw(base_size = 9) +
    theme(
      strip.text       = element_text(size = 10),
      legend.position  = "bottom",
      legend.text      = element_text(size = 12),
      legend.title     = element_text(size = 12),
      axis.title       = element_text(size = 13),
      axis.text        = element_text(size = 11),
      axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_blank()
    )
  if (!is.null(hline)) {
    plot_obj = plot_obj + geom_hline(yintercept = hline, linetype = hline_linetype, color = "gray60", linewidth = 0.4)
  }
  if (has_ci_cols) {
    plot_obj = plot_obj + geom_errorbar(aes(ymin = y_lo, ymax = y_hi), width = 0.2, linewidth = 0.4, position = dodge)
  }
  if (plot) {
    plot(plot_obj)
  }
  if (save_PDF) {
    ggsave(sprintf("%s.pdf", filename_stem), plot_obj, width = 12, height = 6)
  }
  if (save_PNG) {
    ggsave(sprintf("%s.png", filename_stem), plot_obj, width = 12, height = 6, dpi = 150)
  }
  invisible(plot_obj)
}

Nrep = max(results_dt$n_pow, na.rm = TRUE)

for (p_ in unique(results_dt$p)) {
  for (dt_val in unique(results_dt$cond_exp_func_model)) {
    sub = results_dt[p == p_ & cond_exp_func_model == dt_val]
    sub_beta1 = sub[betaT != 0]
    sub_beta0 = sub[betaT == 0]

    # if (nrow(sub_beta1) > 0L) {
    #   stem_pow = sprintf("plot_power_p%s_%s_Nrep_%d",     p_, dt_val, Nrep)
    #   stem_cov = sprintf("plot_coverage_p%s_%s_Nrep_%d",  p_, dt_val, Nrep)
    #   stem_len = sprintf("plot_ci_length_p%s_%s_Nrep_%d", p_, dt_val, Nrep)
    #   make_original_plot(sub_beta1, "pow_avg", "Power", 0.05,
    #     sprintf("Power (betaT=1) | p=%s, log_odds_model=%s", p_, dt_val), stem_pow,
    #     save_PDF = TRUE, plot = FALSE
    #   )
    #   make_original_plot(sub_beta1, "cov_avg", "Coverage", 0.95,
    #     sprintf("Coverage (betaT=1) | p=%s, log_odds_model=%s", p_, dt_val), stem_cov,
    #     save_PDF = TRUE, plot = FALSE
    #   )
    #   make_original_plot(sub_beta1, "len_avg", "Confidence Interval Length", NULL,
    #     sprintf("Confidence Interval Length (betaT=1) | p=%s, log_odds_model=%s", p_, dt_val), stem_len,
    #     save_PDF = TRUE, plot = FALSE, free_y = TRUE
    #   )
    # }
    # 
    # if (nrow(sub_beta0) > 0L) {
    #   stem_size = sprintf("plot_size_p%s_%s_Nrep_%d", p_, dt_val, Nrep)
    #   make_original_plot(sub_beta0, "pow_avg", "Size", 0.05,
    #     sprintf("Size (betaT=0) | p=%s, log_odds_model=%s", p_, dt_val), stem_size,
    #     save_PDF = TRUE, plot = FALSE, hline_linetype = "dotted"
    #   )
    # }

    if (nrow(sub_beta1) > 0L) {
      stem_cov_len = sprintf("plot_coverage_ci_length_p%s_%s_Nrep_%d", p_, dt_val, Nrep)
      make_two_column_plot(
        dat_left = sub_beta1,  metric_left = "cov_avg", xlab_left = "Coverage",  hline_left = 0.95,
        dat_right = sub_beta1, metric_right = "len_avg", xlab_right = "Confidence Interval Length", hline_right = NULL,
        title_str = "", #sprintf("Coverage and Confidence Interval Length (betaT=1) | p=%s, log_odds_model=%s", p_, dt_val),
        filename_stem = stem_cov_len,
        plot = FALSE,
        save_PDF = TRUE,
        free_metric_axis_left = FALSE,
        free_metric_axis_right = TRUE
      )
    }

    if (nrow(sub_beta1) > 0L && nrow(sub_beta0) > 0L) {
      stem_pow_size = sprintf("plot_power_size_p%s_%s_Nrep_%d", p_, dt_val, Nrep)
      make_two_column_plot(
        dat_left = sub_beta1,  metric_left = "pow_avg", xlab_left = "Power", hline_left = NULL,
        dat_right = sub_beta0, metric_right = "pow_avg", xlab_right = "Size", hline_right = 0.05,
        title_str = "", #sprintf("Power and Size | p=%s, log_odds_model=%s", p_, dt_val),
        filename_stem = stem_pow_size,
        plot = FALSE,
        save_PDF = TRUE,
        hline_linetype_left = "dashed",
        hline_linetype_right = "dotted",
        free_metric_axis_left = TRUE,
        free_metric_axis_right = FALSE,
        x_label_digits_left = 2,
        x_label_digits_right = 2,
        x_breaks_right = c(0.05, 0.10, 0.15)
      )
    }
  }
}

n_ = 64
p_ = 1
cond_exp_func_model_ = "linear"
pow =    results_dt[betaT != 0 & n == n_ & p == p_ & cond_exp_func_model == cond_exp_func_model_, 
                 .(power = pow_avg, design = design_short, inference = inference_short)]
pow
size =   results_dt[betaT == 0 & n == n_ & p == p_ & cond_exp_func_model == cond_exp_func_model_, 
                  .(size = pow_avg, design = design_short, inference = inference_short)]
size
cov =    results_dt[betaT != 0 & n == n_ & p == p_ & cond_exp_func_model == cond_exp_func_model_, 
                    .(coverage = cov_avg, design = design_short, inference = inference_short)]
cov
ci_len = results_dt[betaT != 0 & n == n_ & p == p_ & cond_exp_func_model == cond_exp_func_model_, 
                  .(cilength = len_avg, design = design_short, inference = inference_short)]
ci_len

tab_results = merge(
                merge(
                  merge(pow, size, by = c("design", "inference"), all.x = TRUE), 
                    cov, by = c("design", "inference"), all.x = TRUE), 
                      ci_len, by = c("design", "inference"), all.x = TRUE)

design_table_order = c(setdiff(design_labels, "BCRD"), "BCRD")
inference_table_order = c("CMH", "Wald", "ExtendedRobins")

cmh_vs_wald = dcast(
  tab_results[inference %in% c("CMH", "Wald")],
  design ~ inference,
  value.var = c("power", "cilength")
)
cmh_vs_wald[, percent_power_gain := (power_CMH - power_Wald) / power_Wald * 100]
cmh_vs_wald[, percent_length_reduction := (cilength_Wald - cilength_CMH) / cilength_Wald * 100]

tab_results[, design := factor(as.character(design), levels = design_table_order)]
tab_results[, inference := factor(as.character(inference), levels = inference_table_order)]
tab_results = merge(
  tab_results,
  cmh_vs_wald[, .(design, percent_power_gain, percent_length_reduction)],
  by = "design",
  all.x = TRUE
)
tab_results[inference != "CMH", `:=`(
  percent_power_gain = NA_real_,
  percent_length_reduction = NA_real_
)]
setorder(tab_results, design, inference)
tab_results[, c("power", "size", "coverage") := lapply(.SD, function(x) 100 * x), .SDcols = c("power", "size", "coverage")]
numeric_cols = names(tab_results)[vapply(tab_results, is.numeric, logical(1))]
tab_results[, (numeric_cols) := lapply(.SD, signif, digits = 4), .SDcols = numeric_cols]
tab_results

tab_results_latex = copy(tab_results)
tab_results_latex[, `:=`(
  design    = as.character(design),
  inference = as.character(inference)
)]
tab_results_latex[, inference := gsub("ExtendedRobins", "Robins", inference)]
tab_results_latex = tab_results_latex[inference != "Robins"]
tab_results_latex[, design    := gsub("BinaryMatchThenGreedyMD", "BMGMD", design)]
tab_results_latex[, c("size", "coverage") := NULL]
tab_results_latex = tab_results_latex[, .(design, inference, power, cilength, percent_power_gain, percent_length_reduction)]

design_block_ends = which(tab_results_latex[, c(design[-1L] != design[-.N], TRUE)])

latex_table = xtable::xtable(
  as.data.frame(tab_results_latex),
  align = c("l", "l", "l|", "r", "r|", "r", "r")
)

two_row_header = paste0(
  "\\hline\n",
  "Design & Inference & Power & CI & Power & CI Length \\\\\n",
  " & & & Length & Gain (\\%) & Reduction (\\%) \\\\\n",
  "\\hline\n"
)
all_pos  = as.list(c(-1L, design_block_ends, nrow(tab_results_latex)))
all_cmds = c(two_row_header, rep("\\hline\n", length(design_block_ends) + 1L))

print(
  latex_table,
  include.rownames = FALSE,
  include.colnames = FALSE,
  add.to.row = list(pos = all_pos, command = all_cmds),
  hline.after = NULL,
  sanitize.text.function = identity,
  file = sprintf("tab_results_n_%s_p_%s_%s.tex", n_, p_, cond_exp_func_model_)
)

# CMH power vs B for OptimalBlocks at n=64, coloured by p
opt_cmh_n64 = raw_results_dt[
  (grepl("OptimalBlocks", design) | grepl("BCRD", design)) &
  inference == "InferenceIncidCMH" &
  inference_type == "asymp_pval" &
  betaT != 0 &
  n == 64 &
  p %in% c(1, 5, 10) &
  !is.na(pval),
  .(power = mean(pval < 0.05), N = .N),
  by = .(p, design)
]
opt_cmh_n64[design == "DesignFixediBCRD", B := 1]
opt_cmh_n64[design != "DesignFixediBCRD", B := as.integer(sub(".*\\bB=(\\d+).*", "\\1", design))]
opt_cmh_n64[, p_label := factor(paste0("p = ", p), levels = paste0("p = ", c(1, 5, 10)))]
opt_cmh_n64[, pow_lo := power - z * sqrt(power * (1 - power) / N)]
opt_cmh_n64[, pow_hi := power + z * sqrt(power * (1 - power) / N)]

p_cmh_opt_n64 = ggplot(opt_cmh_n64, aes(x = B, y = power, color = p_label, fill = p_label, group = p_label)) +
  geom_ribbon(aes(ymin = pow_lo, ymax = pow_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(opt_cmh_n64$B)), trans = "log2") +
  scale_color_brewer(palette = "Set1", name = "Covariates") +
  scale_fill_brewer(palette = "Set1", name = "Covariates") +
  labs(
    x = "B",
    y = "Power"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("plot_cmh_opt_power_vs_B_n64.pdf", p_cmh_opt_n64, width = 5, height = 4)
p_cmh_opt_n64

