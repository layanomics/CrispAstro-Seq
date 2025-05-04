#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../../config/config.sh"

# -------- OUTPUT DIR --------
MULTIQC_DIR="$PROJECT_DIR/01_preprocessing/results/multiqc"
mkdir -p "$MULTIQC_DIR"

# -------- RUN MULTIQC --------
echo "📊 Running MultiQC on FastQC and fastp results..."

multiqc "$PROJECT_DIR/01_preprocessing/results" -o "$MULTIQC_DIR"

echo "✅ MultiQC report generated in $MULTIQC_DIR"
