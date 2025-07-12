#!/bin/bash
# run_fastqc_trimmed.sh — Runs FastQC on trimmed reads only

# =============================
# 1. Load Configuration & Set Up
# =============================
source ~/CrispAstro-Seq/config/config.sh

# Optimization
MAX_JOBS=6
THREADS_PER_JOB=$((TOTAL_THREADS / MAX_JOBS))
[ "$THREADS_PER_JOB" -lt 1 ] && THREADS_PER_JOB=1

echo -e "\n⚙️ FastQC TRIMMED:"
echo "  - MAX_JOBS: $MAX_JOBS"
echo "  - THREADS_PER_JOB: $THREADS_PER_JOB"

# Setup paths
FASTQC_TYPE="trimmed"
FASTQC_INPUT_DIR="$TRIMMED_DIR"
FASTQC_OUT="$QC_DIR/fastqc_$FASTQC_TYPE"
mkdir -p "$FASTQC_OUT"
LOG_FILE="$FASTQC_OUT/fastqc_${FASTQC_TYPE}_$(date +%Y%m%d_%H%M%S).log"

# Run FastQC
START_TIME=$(date +%s)
CURRENT_JOBS=0

echo -e "\n🔍 Starting FastQC on trimmed reads...\n"

for r1 in "$FASTQC_INPUT_DIR"/*/*_clean_R1.fastq.gz; do
    sample=$(basename "$r1" _clean_R1.fastq.gz)
    r2="$FASTQC_INPUT_DIR/$sample/${sample}_clean_R2.fastq.gz"

    echo "[$(date +'%F %T')] FastQC for $sample (trimmed)..." | tee -a "$LOG_FILE"
    mkdir -p "$FASTQC_OUT/$sample"

    fastqc -t "$THREADS_PER_JOB" "$r1" "$r2" \
        --outdir="$FASTQC_OUT/$sample" >>"$LOG_FILE" 2>&1 &

    ((CURRENT_JOBS++))
    [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ] && wait && CURRENT_JOBS=0
done
wait

END_TIME=$(date +%s)
echo -e "\n✅ FastQC (trimmed) done in $((END_TIME - START_TIME)) sec."
