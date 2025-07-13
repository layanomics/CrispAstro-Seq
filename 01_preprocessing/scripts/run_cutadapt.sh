#!/bin/bash
# run_cutadapt.sh — Trim 10bp from 5' end using Cutadapt for all paired-end samples
# Project: CrispAstro-Seq

set -euo pipefail

# ================================
# Load configuration
# ================================
# shellcheck disable=SC1090
source ~/CrispAstro-Seq/config/config.sh

# Timestamp and logging
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
LOG_FILE="$LOG_DIR/cutadapt_run_$TIMESTAMP.log"
mkdir -p "$TRIMMED_DIR" "$LOG_DIR"

echo -e "\n🚀 Starting Cutadapt 5′-trimming at $TIMESTAMP"
echo "Logging to: $LOG_FILE"
echo "Samples found: $NUM_SAMPLES"

# ================================
# Trim loop for paired-end samples
# ================================
START_TIME=$(date +%s)

for r1 in "$RAW_DIR"/*_1.fastq.gz; do
    sample=$(basename "$r1" | cut -d_ -f1)
    r2="$RAW_DIR/${sample}_2.fastq.gz"

    # Output setup
    SAMPLE_OUT="$TRIMMED_DIR/$sample"
    mkdir -p "$SAMPLE_OUT"
    out1="$SAMPLE_OUT/${sample}_clean_R1.fastq.gz"
    out2="$SAMPLE_OUT/${sample}_clean_R2.fastq.gz"
    report="$SAMPLE_OUT/${sample}_cutadapt_report.txt"

    echo "🧬 Trimming: $sample"

    cutadapt \
        -u 10 -U 10 \
        -o "$out1" \
        -p "$out2" \
        "$r1" "$r2" \
        >"$report" 2>>"$LOG_FILE"

    echo "✅ Done: $sample"
done

echo -e "\n🎉 Cutadapt trimming completed for all samples."
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo -e "\n⏱️  Total Runtime: $RUNTIME_FMT"
