#!/bin/bash
# Reruns perf annotate only (perf stat + perf record already done) with 30min timeout.
set -uo pipefail

TIMEOUT=1800

mapfile -t KERNELS < /tmp/timed_out_kernels.txt
TOTAL=${#KERNELS[@]}
IDX=0

for KERNEL in "${KERNELS[@]}"; do
    IDX=$((IDX + 1))
    DATA_OUT="/tmp/perf_${KERNEL}.data"
    ANNO_OUT="/tmp/perf_annotate_${KERNEL}.txt"

    echo "[$IDX/$TOTAL] === $KERNEL ==="

    if [[ ! -f "$DATA_OUT" ]]; then
        echo "  SKIP: no .data file"
        continue
    fi

    echo "  perf annotate -> $ANNO_OUT"
    timeout "$TIMEOUT" perf annotate \
        --stdio \
        -i "$DATA_OUT" \
        --no-source \
        2>/dev/null \
        > "$ANNO_OUT" \
        && echo "  Done $KERNEL" \
        || echo "  WARN: perf annotate timed out or failed for $KERNEL"
done

echo ""
echo "=== Retry annotate complete for $TOTAL kernels ==="
