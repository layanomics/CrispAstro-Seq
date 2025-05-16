#!/bin/bash
# run_fastqc.sh — Runs FastQC in parallel per sample on both raw and trimmed reads
#                 using dynamic threading and per-sample output organization

# =============================
# 1. Load Configuration & Set Up
# =============================
# Load project paths and resource logic (THREADS, MAX_JOBS, etc.)
# shellcheck disable=SC1090
source ~/CrispAstro-Seq/config/config.sh

# =============================
# [FASTQC OPTIMIZATION] Override Parallel Settings for FastQC
# =============================
# FastQC benefits from running many lightweight jobs

MAX_JOBS=6
THREADS_PER_JOB=$((TOTAL_THREADS / MAX_JOBS))
[ "$THREADS_PER_JOB" -lt 1 ] && THREADS_PER_JOB=1

echo -e "\n⚙️ FastQC override:"
echo "  - MAX_JOBS (FastQC)       : $MAX_JOBS"
echo "  - THREADS_PER_JOB (FastQC): $THREADS_PER_JOB"

# =============================
# 2. Pass 1 — Run FastQC on RAW Reads
# =============================
START_TIME=$(date +%s)
FASTQC_TYPE="raw"
FASTQC_INPUT_DIR="$RAW_DIR"
FASTQC_OUT="$QC_DIR/fastqc_$FASTQC_TYPE"
mkdir -p "$FASTQC_OUT"
LOG_FILE="$FASTQC_OUT/fastqc_${FASTQC_TYPE}_$(date +%Y%m%d_%H%M%S).log"

echo -e "\n🔍 Starting FastQC on $FASTQC_TYPE reads...\n"

CURRENT_JOBS=0

for r1 in "$FASTQC_INPUT_DIR"/*_1.fastq.gz; do
    sample=$(basename "$r1" | cut -d'_' -f1)
    r2="$FASTQC_INPUT_DIR/${sample}_2.fastq.gz"

    echo "[$(date +'%F %T')] FastQC for $sample ($FASTQC_TYPE)..." | tee -a "$LOG_FILE"
    mkdir -p "$FASTQC_OUT/$sample"

    fastqc -t "$THREADS_PER_JOB" \
        "$r1" "$r2" \
        --outdir="$FASTQC_OUT/$sample" >>"$LOG_FILE" 2>&1 &

    ((CURRENT_JOBS++))
    [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ] && wait && CURRENT_JOBS=0
done
wait

# =============================
# 3. Pass 2 — Run FastQC on TRIMMED Reads
# =============================
FASTQC_TYPE="trimmed"
FASTQC_INPUT_DIR="$TRIMMED_DIR"
FASTQC_OUT="$QC_DIR/fastqc_$FASTQC_TYPE"
mkdir -p "$FASTQC_OUT"
LOG_FILE="$FASTQC_OUT/fastqc_${FASTQC_TYPE}_$(date +%Y%m%d_%H%M%S).log"

echo -e "\n🔍 Starting FastQC on $FASTQC_TYPE reads...\n"

CURRENT_JOBS=0

for r1 in "$FASTQC_INPUT_DIR"/*/*_clean_R1.fastq.gz; do
    sample=$(basename "$r1" _clean_R1.fastq.gz)
    r2="$FASTQC_INPUT_DIR/$sample/${sample}_clean_R2.fastq.gz"

    echo "[$(date +'%F %T')] FastQC for $sample ($FASTQC_TYPE)..." | tee -a "$LOG_FILE"
    mkdir -p "$FASTQC_OUT/$sample"

    fastqc -t "$THREADS_PER_JOB" \
        "$r1" "$r2" \
        --outdir="$FASTQC_OUT/$sample" >>"$LOG_FILE" 2>&1 &

    ((CURRENT_JOBS++))
    [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ] && wait && CURRENT_JOBS=0
done
wait

# =============================
# 4. Finalize and Report Runtime
# =============================
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo -e "\n✅ FastQC runs (raw + trimmed) complete."
echo "Total Runtime: $RUNTIME_FMT"
