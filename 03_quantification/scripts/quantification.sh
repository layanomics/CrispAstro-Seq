#!/bin/bash

# ---------- CONFIG ----------
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../../config/config.sh"

# Create output directory if it doesn't exist
mkdir -p "$COUNTS_DIR"

# Run featureCounts for paired-end RNA-seq data
featureCounts \
    -T "$THREADS" \
    -a "$GENOME_GTF" \
    -o "$COUNTS_DIR/gene_counts.txt" \
    -g gene_id \
    -t exon \
    -p \
    -B \
    -C \
    "$ALIGN_DIR"/*.bam

echo "✅ Quantification complete. Results saved in $COUNTS_DIR/gene_counts.txt"
