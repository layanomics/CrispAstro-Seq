#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../../config/config.sh"

# -------- OUTPUT DIR --------
FASTQC_OUT="$PROJECT_DIR/01_preprocessing/results/fastqc_raw"
mkdir -p "$FASTQC_OUT"

# -------- RUN FASTQC --------
echo "🔬 Running FastQC on all raw FASTQ files using $THREADS threads..."

fastqc -t "$THREADS" "$RAW_DIR"/*.fastq.gz --outdir="$FASTQC_OUT"

echo "✅ FastQC completed for all files."
