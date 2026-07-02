#!/bin/bash
# EDI Benchmark Profiling Runner
# Generates perf stat + perf record + perf annotate for every EDI C++ path.
#
# Outputs per kernel (KERNEL = kernel name):
#   /tmp/perf_stat_<KERNEL>.txt      -- perf stat hardware counters
#   /tmp/perf_<KERNEL>.data          -- perf record raw samples
#   /tmp/perf_annotate_<KERNEL>.txt  -- perf annotate --stdio output
#
# All R scripts use seed=42; see profile/edi_kernel_profiler.R for full code.

set -uo pipefail

RSCRIPT="Rscript --no-save --no-restore"
PROFILER="$(dirname "$0")/edi_kernel_profiler.R"
PERF_FREQ=199
ANNO_TIMEOUT=1800   # 30 minutes for perf annotate

# All kernel names to profile
KERNELS=(
    # --- estimate-only paths ---
    logistic_est
    ols_est
    poisson_est
    negbin_est
    beta_est
    robust_est
    logbin_est
    probit_est
    identbin_est
    hurdle_p_est
    zip_est
    zinb_est
    hurdle_nb_est
    prop_odds_est
    adj_cat_est
    cont_ratio_est
    ord_probit_est
    ord_cloglog_est
    ord_cauchit_est
    coxph_est
    strat_coxph_est
    weibull_est
    logrank_est
    km_diff_est
    rmean_diff_est
    wilcox_est
    gcomp_logistic_rd_est
    gcomp_logistic_rr_est
    gcomp_frac_logit_est
    gcomp_ordinal_est
    # --- variance / full-inference paths ---
    ols_var
    ols_lin_var
    robust_var
    logistic_var
    logbin_var
    probit_var
    newcombe_var
    mn_var
    poisson_var
    quasi_var
    negbin_var
    hurdle_nb_var
    hurdle_p_var
    zip_var
    zinb_var
    poisson_robust_var
    beta_var
    prop_odds_var
    adj_cat_var
    cont_ratio_var
    ridit_var
    jt_var
    coxph_var
    strat_coxph_var
    weibull_var
    logrank_var
    wilcox_var
    gcomp_logistic_post_fit_var
    gcomp_ordinal_var
    # --- weighted estimate-only paths ---
    logistic_weighted_est
    poisson_weighted_est
    probit_weighted_est
    beta_weighted_est
    logbin_weighted_est
    identbin_weighted_est
    prop_odds_weighted_est
    # --- missing full-inference paths ---
    identbin_var_full
    ord_cauchit_var
    ord_cloglog_var
    ord_probit_var
    stereotype_est
    stereotype_var
    # --- new model types ---
    trunc_negbin_est
    zero_one_inflated_beta_est
    weibull_frailty_est
    weibull_frailty_var
    # --- GLMM paths ---
    logistic_glmm_est
    logistic_glmm_var
    poisson_glmm_est
    poisson_glmm_var
    gaussian_lmm_est
    gaussian_lmm_var
    ordinal_glmm_est
    ordinal_glmm_var
    ordinal_clmm_est
    hurdle_p_glmm_est
    hurdle_p_glmm_var
    # --- GEE paths ---
    gee_pairs_singletons_logistic
    gee_pairs_singletons_weighted_logistic
    # --- KK21 weight functions ---
    kk21_continuous_wts
    kk21_logistic_wts
    kk21_beta_wts
    kk21_negbin_wts
    kk21_ordinal_wts
    kk21_survival_wts
    kk21_stepwise_continuous_wts
    kk21_stepwise_logistic_wts
    # --- stats helpers ---
    newcombe_paired
    mn_ci
    zhang_binom_pval
    zhang_fisher_pval
    # --- post-fit variance helpers ---
    gcomp_frac_logit_post_fit_var
    glm_sandwich_post_fit_var
    ordinal_gcomp_post_fit_var
    # --- bootstrap index generators ---
    bootstrap_indices
    stratified_bootstrap_indices
    bootstrap_m_indices
    # --- randomization primitives ---
    complete_randomization_balanced
    complete_randomization_imbalanced
    shuffle_w
    efron_redraw
    spbr_redraw
    random_block_size_redraw
    redraw_w_kk14
    atkinson_redraw
    pocock_simon_assign
    pocock_simon_assign_and_update
    pocock_simon_redraw_w
    # --- permutation generators ---
    generate_permutations_bernoulli
    generate_permutations_efron
    generate_permutations_ibcrd
    generate_permutations_blocking
    generate_permutations_matching
    generate_permutations_atkinson
    generate_permutations_pocock_simon
    generate_permutations_cluster
    generate_permutations_spbr
    # --- bootstrap/randomization loops ---
    draw_binary_match_assignments
    draw_matching_bootstrap_sample
    randomization_loop
    base_bootstrap_loop
)

TOTAL=${#KERNELS[@]}
IDX=0

for KERNEL in "${KERNELS[@]}"; do
    IDX=$((IDX + 1))
    STAT_OUT="/tmp/perf_stat_${KERNEL}.txt"
    DATA_OUT="/tmp/perf_${KERNEL}.data"
    ANNO_OUT="/tmp/perf_annotate_${KERNEL}.txt"

    echo "[$IDX/$TOTAL] === $KERNEL ==="

    # ---- perf stat ----
    echo "  perf stat -> $STAT_OUT"
    perf stat \
        --output="$STAT_OUT" \
        -- $RSCRIPT "$PROFILER" "$KERNEL" \
        > /tmp/perf_stat_stdout_${KERNEL}.txt 2>&1 \
        || { echo "  WARN: perf stat failed for $KERNEL"; }

    # ---- perf record ----
    echo "  perf record -> $DATA_OUT"
    perf record \
        -F "$PERF_FREQ" \
        --output="$DATA_OUT" \
        -q \
        -- $RSCRIPT "$PROFILER" "$KERNEL" \
        > /tmp/perf_record_stdout_${KERNEL}.txt 2>&1 \
        || { echo "  WARN: perf record failed for $KERNEL"; continue; }

    # ---- perf annotate ----
    echo "  perf annotate -> $ANNO_OUT"
    timeout "$ANNO_TIMEOUT" perf annotate \
        --stdio \
        -i "$DATA_OUT" \
        --no-source \
        2>/dev/null \
        > "$ANNO_OUT" \
        || { echo "  WARN: perf annotate timed out or failed for $KERNEL"; }

    echo "  Done $KERNEL"
done

echo ""
echo "=== All $TOTAL kernels done ==="
echo "Outputs:"
echo "  perf stat:     /tmp/perf_stat_<kernel>.txt"
echo "  perf record:   /tmp/perf_<kernel>.data"
echo "  perf annotate: /tmp/perf_annotate_<kernel>.txt"
