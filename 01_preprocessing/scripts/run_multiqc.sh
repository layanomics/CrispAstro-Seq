#!/bin/bash
# run_multiqc.sh — Aggregates QC results using MultiQC and saves timestamped report

# =============================
# 1. Load Configuration & Set Up
# =============================
# shellcheck disable=SC1090
source ~/CrispAstro-Seq/config/config.sh

# =============================
# 2. Handle Optional --tag Argument with Error Checking
# =============================
# Default tag value
TAG="qc"

# Check if user passed --tag and a value
if [[ "$1" == "--tag" ]]; then
    if [[ -n "$2" ]]; then
        TAG="$2"
    else
        echo "❌ You used --tag but didn't give a tag name."
        echo "✅ Example: bash run_multiqc.sh --tag post_trim"
        exit 1
    fi
elif [[ -n "$1" ]]; then
    echo "❌ Unknown option: $1"
    echo "✅ Usage: bash run_multiqc.sh --tag your_label"
    exit 1
fi

# =============================
# 3. Set Output Directory with Timestamp and Tag
# =============================
# Start time for duration tracking
START_TIME=$(date +%s)

# Create timestamped output directory for MultiQC summary
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUT_DIR="$QC_DIR/multiqc_${TAG}_$TIMESTAMP"
mkdir -p "$OUT_DIR"

echo "📊 Running MultiQC..."
echo "➤ Tag             : $TAG"
echo "➤ Input Directory : $QC_DIR"
echo "➤ Output Directory: $OUT_DIR"

# =============================
# 4. Run MultiQC
# =============================
multiqc "$QC_DIR" -o "$OUT_DIR"

# =============================
# 5. Finalize and Report Runtime
# =============================
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo "✅ MultiQC Report saved to: $OUT_DIR"
echo "Total Runtime: $RUNTIME_FMT"
