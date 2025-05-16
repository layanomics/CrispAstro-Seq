#!/bin/bash
# run_fastp.sh — Performs Fastp QC on paired-end samples with parallelism and logging

# =============================
# 1. Load Configuration & Set Up
# =============================
# shellcheck disable=SC1090
source ~/CrispAstro-Seq/config/config.sh
# =============================
# [FASTP OPTIMIZATION] Override Parallel Settings for Fastp
# =============================
# Fastp benefits more from multithreading within fewer jobs

# Manually override MAX_JOBS for Fastp
MAX_JOBS=4
THREADS_PER_JOB=2
[ "$THREADS_PER_JOB" -lt 1 ] && THREADS_PER_JOB=1

echo -e "\n⚙️ Fastp override:"
echo "  - MAX_JOBS (Fastp)       : $MAX_JOBS"
echo "  - THREADS_PER_JOB (Fastp): $THREADS_PER_JOB"

# Start time for runtime tracking
START_TIME=$(date +%s)

# Ensure output directory exists
mkdir -p "$TRIMMED_DIR"

# Optional: log file if you want to keep run logs
LOG_FILE="$QC_DIR/fastp_run_$(date +%Y%m%d_%H%M%S).log"

# =============================
# 2. Loop Through Samples and Run Fastp in Parallel
# =============================

CURRENT_JOBS=0

# Loop over all R1 reads to identify samples
for r1 in "$RAW_DIR"/*_1.fastq.gz; do
    sample=$(basename "$r1" | cut -d'_' -f1)
    r2="$RAW_DIR/${sample}_2.fastq.gz"

    # Validate file pair existence and non-zero size
    if [[ ! -f "$r2" || ! -s "$r1" || ! -s "$r2" ]]; then
        echo "❌ Skipping $sample — missing or empty read pair" | tee -a "$LOG_FILE"
        continue
    fi

    echo "[$(date +'%F %T')] Fastp for $sample..." | tee -a "$LOG_FILE"

    # Create subfolder for each sample
    SAMPLE_DIR="$TRIMMED_DIR/$sample"
    mkdir -p "$SAMPLE_DIR"

    # Run fastp with multithreading, logging to file
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$SAMPLE_DIR/${sample}_clean_R1.fastq.gz" \
        -O "$SAMPLE_DIR/${sample}_clean_R2.fastq.gz" \
        -h "$SAMPLE_DIR/${sample}_fastp.html" \
        -j "$SAMPLE_DIR/${sample}_fastp.json" \
        --thread "$THREADS_PER_JOB" >>"$LOG_FILE" 2>&1 &

    ((CURRENT_JOBS++))
    [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ] && wait && CURRENT_JOBS=0
done

# Final wait to ensure all background jobs complete
wait

# =============================
# 3. Finalize and Report Runtime
# =============================

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo "===== Done: $(date) =====" | tee -a "$LOG_FILE"
echo "Total Runtime: $RUNTIME_FMT" | tee -a "$LOG_FILE"
