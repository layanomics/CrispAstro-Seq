#!/bin/bash

# ---------- CONFIG ----------
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../../config/config.sh"

# ---------- DYNAMIC THREAD ADJUSTMENT ----------
AVAILABLE_THREADS=$(nproc)
if [ "$THREADS" -gt "$AVAILABLE_THREADS" ]; then
    echo "⚠️ Configured THREADS=$THREADS exceeds available cores ($AVAILABLE_THREADS). Using $AVAILABLE_THREADS threads."
    THREADS=$AVAILABLE_THREADS
fi

# ---------- LOG SETUP ----------
LOG_DIR="$SCRIPT_DIR/../../logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/quantification_$(date +%Y%m%d_%H%M%S).log"

# ---------- Run featureCounts ----------
mkdir -p "$COUNTS_DIR"

featureCounts \
    -T "$THREADS" \
    -a "$GENOME_GTF" \
    -o "$COUNTS_DIR/gene_counts.txt" \
    -g gene_id \
    -t exon \
    -p \
    -B \
    -C \
    "$ALIGN_DIR"/*.bam | tee "$LOG_FILE"

echo "✅ Quantification complete. Results saved in $COUNTS_DIR/gene_counts.txt"
echo "📄 Log saved to $LOG_FILE"
