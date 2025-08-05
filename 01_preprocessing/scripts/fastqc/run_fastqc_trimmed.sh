#!/bin/bash
set -euo pipefail

######################################################################
# 🧬 CrispAstro-Seq Pipeline
# Script        : run_fastqc_trimmed.sh
# Description   : Run FastQC on trimmed paired-end reads in parallel
# Author        : Layan Essam
# Tool Version  : FastQC v0.12.1
# Last Updated  : July 2025
######################################################################

# =============================
# Load Configuration
# =============================
source "$(dirname "$0")/../../config/config.sh"
source "$(dirname "$0")/../../config/utils.sh"

# =============================
# 🔧 Tool Version Check
# =============================
check_tool_version fastqc

# =============================
# 🕒 Timestamp Initialization
# =============================
initialize_timestamps

# =============================
# 📁 Output Paths Setup
# =============================
FASTQC_TYPE="trimmed"
FASTQC_INPUT_DIR="$TRIMMED_DIR"
FASTQC_OUTPUT_DIR="$FASTQC_TRIMMED_DIR"
mkdir -p "$FASTQC_OUTPUT_DIR"
LOG_FILE="$FASTQC_OUTPUT_DIR/fastqc_${FASTQC_TYPE}_${TIMESTAMP}.log"

# =============================
# Sample, Parallelism & Thread Logic
# =============================
TOTAL_THREADS=$(nproc)
NUM_SAMPLES=$(find "$FASTQC_INPUT_DIR" -name "*_clean_R1.fastq.gz" | wc -l)
MAX_JOBS=$(( NUM_SAMPLES > 0 ? ( NUM_SAMPLES < TOTAL_THREADS ? NUM_SAMPLES : TOTAL_THREADS ) : 1 ))

# =============================
# 📝 Session Log Header
# =============================
{
  echo ""
  echo "====================== 🧪 FASTQC ANALYSIS ======================="
  echo "⚙️ FastQC CONFIGURATION [$FASTQC_TYPE]"
  echo "🔧 Version               : $TOOL_VERSION"
  echo "🧬 Samples Detected      : $NUM_SAMPLES"
  echo "🖥️ Total Threads         : $TOTAL_THREADS"
  echo "⚡ Max Parallel Jobs     : $MAX_JOBS"
  echo "📂 Input Directory       : $FASTQC_INPUT_DIR"
  echo "📁 Output Directory      : $FASTQC_OUTPUT_DIR"
  echo "🕒 Script launched at    : $START_TIME_KSA | UTC: $START_TIME_UTC"
  echo "🧾 Log File              : $LOG_FILE"
  echo "================================================================"
} | tee -a "$LOG_FILE"

echo -e "\n🧬 Starting FastQC ($FASTQC_TYPE) for $NUM_SAMPLES samples..."

# =============================
# 🚀 FastQC Parallel Execution
# =============================
CURRENT_JOBS=0
SUCCESS_COUNT=0
START_TIME=$(date +%s)

for r1 in "$FASTQC_INPUT_DIR"/*/*_clean_R1.fastq.gz; do
    sample=$(basename "$r1" | sed -E 's/_clean_R1\.fastq\.gz$//')
    r2="${r1/_clean_R1.fastq.gz/_clean_R2.fastq.gz}"

    SAMPLE_OUT="$FASTQC_OUTPUT_DIR/$sample"
    mkdir -p "$SAMPLE_OUT"
    SAMPLE_LOG="$SAMPLE_OUT/fastqc_${sample}_${TIMESTAMP}.log"

    echo -e "\n[$(TZ='Asia/Riyadh' date +'%Y-%m-%d %H:%M')] 🚀 FastQC: $sample" | tee -a "$LOG_FILE"

    (
        if fastqc -t 1 "$r1" "$r2" --outdir="$SAMPLE_OUT"; then
            echo "✅ FastQC succeeded for $sample"
            echo "✅ $sample passed" >> "$LOG_FILE"
            touch "$SAMPLE_OUT/success.flag"
        else
            echo "❌ FastQC failed for $sample"
            echo "❌ $sample failed" >> "$LOG_FILE"
        fi
    ) >> "$SAMPLE_LOG" 2>&1 &

    # If the number of active jobs reaches the allowed maximum...
    # Then: wait for all to finish, and reset the counter to zero
    ((CURRENT_JOBS++))
    [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ] && wait && CURRENT_JOBS=0
done

wait

# =============================
# ✅ Final Summary & Runtime
# =============================
SUCCESS_COUNT=$(find "$FASTQC_OUTPUT_DIR" -name success.flag | wc -l)
log_runtime

{
  echo ""
  echo "====================== ✅ FASTQC COMPLETED ======================="
  echo "🎯 Successfully processed : $SUCCESS_COUNT of $NUM_SAMPLES samples"
  echo "📄 Reports directory      : $FASTQC_OUTPUT_DIR"
  echo "📄 Example report         : $(find "$FASTQC_OUTPUT_DIR" -name '*_fastqc.html' | head -n 1)"
  print_runtime_summary
  if [[ "$SUCCESS_COUNT" -eq "$NUM_SAMPLES" ]]; then
    echo "✅ Step completed successfully: run_fastqc_trimmed"
  else
    echo "⚠️  Step completed with some errors: $((NUM_SAMPLES - SUCCESS_COUNT)) failed sample(s)"
  fi
  echo "=================================================================="
} | tee -a "$LOG_FILE"

