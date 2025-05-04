#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../../config/config.sh"

# -------- OUTPUT DIR --------
FASTP_OUT="$PROJECT_DIR/01_preprocessing/results/fastp_cleaned"
mkdir -p "$FASTP_OUT"

# -------- RUN FASTP --------
echo "✂️ Running fastp on all paired-end FASTQ files with $THREADS threads..."

for r1 in "$RAW_DIR"/*_1.fastq.gz; do
    base=$(basename "$r1" _1.fastq.gz)
    r2="$RAW_DIR/${base}_2.fastq.gz"

    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$FASTP_OUT/${base}_clean_R1.fastq.gz" \
        -O "$FASTP_OUT/${base}_clean_R2.fastq.gz" \
        -q 20 \
        -l "$MIN_LENGTH" \
        -w "$THREADS" \
        -h "$FASTP_OUT/${base}_fastp.html" \
        -j "$FASTP_OUT/${base}_fastp.json"
done

echo "✅ fastp completed for all samples."
