#!/bin/bash
set -euo pipefail

######################################################################
# 🧬 CrispAstro-Seq Pipeline
# Script        : run_cutadapt.sh
# Description   : Trim 10bp from 5' end using Cutadapt for all paired-end samples
# Author        : Layan Essam
# Tool Version  : Cutadapt v4.6
# Last Updated  : July 2025
######################################################################

# =============================
# ⚙️ Load Configuration
# =============================
source "$HOME/CrispAstro-Seq/config/config.sh"
source "$PROJECT_DIR/config/utils.sh"

# =============================
# 🔧 Tool Version Check
# =============================
check_tool_version cutadapt

# =============================
# 🕒 Timestamp Initialization
# =============================
initialize_timestamps

# =============================
# ✂️ Trimming Strategy Tag
# =============================
CUTADAPT_INPUT_DIR="$RAW_DIR"
CUTADAPT_OUTPUT_DIR="$TRIMMED_DIR"
mkdir -p "$CUTADAPT_OUTPUT_DIR" "$CUTADAPT_LOG_DIR"
CUTADAPT_LOG_FILE="$CUTADAPT_LOG_DIR/cutadapt_trim10_${TIMESTAMP}.log"

# =============================
# 📝 Session Log Header
# =============================

NUM_SAMPLES=$(find "$CUTADAPT_INPUT_DIR" -name "*_1.fastq.gz" | wc -l)

{
  echo ""
  echo "====================== ✂️ CUTADAPT TRIMMING ======================="
  echo "⚙️ Cutadapt CONFIGURATION"
  echo "🔧 Version               : $TOOL_VERSION"
  echo "🏷️ Trimming Tag          : 10bp from 5′ ends"
  echo "🧬 Samples Detected      : $NUM_SAMPLES"
  echo "📂 Input Directory       : $CUTADAPT_INPUT_DIR"
  echo "📁 Output Directory      : $CUTADAPT_OUTPUT_DIR"
  echo "🕒 Script launched at    : $START_TIME_KSA | UTC: $START_TIME_UTC"
  echo "🧾 Log File              : $CUTADAPT_LOG_FILE"
  echo "==================================================================="
} | tee -a "$CUTADAPT_LOG_FILE"

# =============================
# 🚀 Run Cutadapt for Each Sample
# =============================
echo -e "\n✂️ Starting Cutadapt trimming for $NUM_SAMPLES samples..."

START_TIMESTAMP=$(date +%s)
SUCCESS_COUNT=0

for r1 in "$RAW_DIR"/*_1.fastq.gz; do
    sample=$(basename "$r1" | sed -E 's/_1\.fastq\.gz$//')
    r2="$CUTADAPT_INPUT_DIR/${sample}_2.fastq.gz"

    SAMPLE_OUT="$CUTADAPT_OUTPUT_DIR/$sample"
    mkdir -p "$SAMPLE_OUT"

    out1="$SAMPLE_OUT/${sample}_clean_R1.fastq.gz"
    out2="$SAMPLE_OUT/${sample}_clean_R2.fastq.gz"
    report="$CUTADAPT_LOG_DIR/${sample}_cutadapt_report.txt"
    sample_log="$CUTADAPT_LOG_DIR/cutadapt_${sample}_${TIMESTAMP}.log"

    echo -e "\n[$(TZ='Asia/Riyadh' date +'%Y-%m-%d %H:%M') KSA (+03:00)] ✂️ Trimming: $sample" | tee -a "$sample_log" "$CUTADAPT_LOG_FILE"
    
    if [[ ! -f "$r2" || ! -s "$r1" || ! -s "$r2" ]]; then
        echo "❌ Skipping $sample — missing or empty reads" | tee -a "$sample_log" "$CUTADAPT_LOG_FILE"
        continue
    fi
    
    cutadapt \
        -u 10 -U 10 \
        -o "$out1" \
        -p "$out2" \
        "$r1" "$r2" \
        > "$report" 2>&1

    # shellcheck disable=SC2181
    if [ $? -ne 0 ]; then
        echo "❌ Cutadapt failed for $sample. See log: $sample_log" | tee -a "$CUTADAPT_LOG_FILE"
    else
        echo "✅ Done: $sample" | tee -a "$CUTADAPT_LOG_FILE"
        ((SUCCESS_COUNT++))
    fi
done

# =============================
# ✅ Final Summary
# =============================
log_runtime

{
  echo ""
  echo "====================== ✅ CUTADAPT COMPLETED ======================="
  echo "🎯 Successfully trimmed   : $SUCCESS_COUNT of $NUM_SAMPLES samples"
  echo "📄 Cleaned FASTQs in      : $CUTADAPT_OUTPUT_DIR/<sample>/*_clean_R*.fastq.gz"
  echo "📄 Reports/logs in        : $CUTADAPT_LOG_DIR/"
  print_runtime_summary
  echo "✅ Step completed successfully: run_cutadapt"
  echo "==================================================================="
} | tee -a "$CUTADAPT_LOG_FILE"


