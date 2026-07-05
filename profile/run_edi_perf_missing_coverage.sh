#!/bin/bash
# Profiles the C++ paths identified in TODO-98 that had zero coverage in the original
# 129-kernel sweep.  Covers: ClogitPlusGLMM, dep-cens-transform survival,
# d-optimal search, KK compound distribution, BAI parallel, rerandomization helpers,
# and CMH block SE.

set -uo pipefail

RSCRIPT="Rscript --no-save --no-restore"
PROFILER="$(dirname "$0")/edi_kernel_profiler.R"
PERF_FREQ=199
ANNO_TIMEOUT=1800   # 30 minutes

KERNELS=(
    clogit_glmm_est
    clogit_glmm_var
    dep_cens_transform_est
    dep_cens_transform_var
    d_optimal_search
    kk_compound_distr
    bai_distr
    rerandomization_search
    rerandomization_obj_vals
    cmh_block_se
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
echo "=== All $TOTAL missing-coverage kernels done ==="
echo "Outputs:"
echo "  perf stat:     /tmp/perf_stat_<kernel>.txt"
echo "  perf record:   /tmp/perf_<kernel>.data"
echo "  perf annotate: /tmp/perf_annotate_<kernel>.txt"
