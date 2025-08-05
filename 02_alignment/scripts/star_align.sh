#!/bin/bash
set -euo pipefail
######################################################################
# 🧬 CrispAstro-Seq Pipeline
# Script        : star_align.sh
# Description   : STAR alignment for trimmed FASTQ reads (paired-end)
# Author        : Layan Essam
# Tool Version  : 2.6.0a
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
check_tool_version STAR

# =============================
# 🕒 Timestamp Initialization
# =============================
initialize_timestamps

# =============================
# 📁 Output Paths Setup
# =============================
ALIGN_INPUT_DIR="$TRIMMED_DIR"
ALIGN_OUTPUT_DIR="$ALIGN_DIR"
mkdir -p "$ALIGN_OUTPUT_DIR" "$ALIGN_LOG_DIR"
ALIGN_LOG="$ALIGN_LOG_DIR/star_align_cutadapt_trim10_${TIMESTAMP}.log"

# =============================
# 🧠 System Resource Detection
# =============================
TOTAL_THREADS=$(nproc); THREADS=$TOTAL_THREADS
# Memory in bytes for BAM sorting (80% of total RAM)
TOTAL_MEM_BYTES=$(grep MemTotal /proc/meminfo | awk '{print $2 * 1024}')
BAM_SORT_RAM=$((TOTAL_MEM_BYTES * 80 / 100))

# =============================
# 📝 Session Log Header
# =============================
NUM_SAMPLES=$(find "$ALIGN_INPUT_DIR" -name "*_clean_R1.fastq.gz" | wc -l)
{
  echo ""
  echo "===================== 🧬 STAR ALIGNMENT ====================="
  echo "⚙️ STAR CONFIGURATION"
  echo "🔢 STAR Version         : $TOOL_VERSION"
  echo "🧠 Threads              : $THREADS"
  echo "💾 BAM Sort RAM         : $BAM_SORT_RAM bytes"
  echo ""
  echo "📋 STAR alignment setup:"
  echo "🕒 Script launched at   : $START_TIME_KSA | UTC: $START_TIME_UTC"
  echo "🔢 Samples to align     : $NUM_SAMPLES"
  echo "📂 Input Directory      : $ALIGN_INPUT_DIR"
  echo "📁 Output Directory     : $ALIGN_OUTPUT_DIR"
  echo "📄 Logging to           : $ALIGN_LOG"
  echo "==============================================================="
} | tee -a "$ALIGN_LOG"

# =============================
# 🚀 STAR Alignment Execution
# =============================
SUCCESS_COUNT=0
START_TIMESTAMP=$(date +%s)

for r1 in "$ALIGN_INPUT_DIR"/*/*_clean_R1.fastq.gz; do
    sample=$(basename "$r1" | sed -E 's/_clean_R1\.fastq\.gz$//')
    r2="${r1/_clean_R1.fastq.gz/_clean_R2.fastq.gz}"

    SAMPLE_OUT_DIR="$ALIGN_OUTPUT_DIR/$sample"
    mkdir -p "$SAMPLE_OUT_DIR"
    SAMPLE_LOG="$ALIGN_LOG_DIR/star_align_${sample}_${TIMESTAMP}.log"

    {
        SAMPLE_START_KSA=$(TZ=Asia/Riyadh date +'%H:%M')
        echo "🕒 Starting STAR for $sample at $SAMPLE_START_KSA (KSA, UTC+3)"

        if [[ ! -f "$r2" || ! -s "$r1" || ! -s "$r2" ]]; then
            echo "❌ Skipping $sample — missing or empty reads" | tee -a "$ALIGN_LOG"
            continue
        fi

        SAMPLE_START_TS=$(date +%s)

        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX_DIR" \
            --readFilesIn "$r1" "$r2" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$SAMPLE_OUT_DIR/${sample}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --quantMode GeneCounts \
            --twopassMode Basic \
            --limitBAMsortRAM "$BAM_SORT_RAM" 
        
        STAR_EXIT_CODE=$?
        SAMPLE_END_TS=$(date +%s)


        echo "✅ Finished $sample in $(date -ud "@$((SAMPLE_END_TS - SAMPLE_START_TS))" +'%Hh %Mm')" | tee -a "$ALIGN_LOG"

        if [[ $STAR_EXIT_CODE -eq 0 ]]; then
            echo "✅ $sample passed — STAR alignment completed" | tee -a "$ALIGN_LOG"
            touch "$SAMPLE_OUT_DIR/success.flag"
            ((SUCCESS_COUNT++))
        else
            echo "❌ $sample failed — STAR encountered an error" | tee -a "$ALIGN_LOG"
        fi

        
        echo "--------------------------------------------------------"
        echo ""
    } >> "$SAMPLE_LOG" 2>&1
done

# =============================
# ✅ Alignment Output Validation
# =============================
{
  echo ""
  echo "==================== 🔍 ALIGNMENT OUTPUT VALIDATION ===================="
  echo ""

  FAILED_SAMPLES=0

  for sample_dir in "$ALIGN_DIR"/*; do
    [[ -d "$sample_dir" ]] || continue
    bam_file=$(find "$sample_dir" -name "*Aligned.sortedByCoord.out.bam" -type f)

    if [[ -s "$bam_file" ]]; then
      echo "✅ $(basename "$sample_dir"): BAM file present and non-empty."
    else
      echo "❌ $(basename "$sample_dir"): Missing or empty BAM file."
      ((FAILED_SAMPLES++))
    fi
  done

  echo ""

  if [[ $FAILED_SAMPLES -gt 0 ]]; then
    echo "❌ Alignment output check failed for $FAILED_SAMPLES sample(s)."
  else
    echo "✅ All aligned BAM files validated successfully."
  fi

} | tee -a "$ALIGN_LOG"

# =============================
# Finalize and Runtime Summary
# =============================
log_runtime

{
  echo ""
  echo "=================== ✅ STAR ALIGNMENT COMPLETED ==================="
  echo ""
  echo "🎯 Successfully aligned $SUCCESS_COUNT of $NUM_SAMPLES samples."
  print_runtime_summary
  echo "📄 Log file saved to : $ALIGN_LOG"
  echo "=============================================================="
} | tee -a "$ALIGN_LOG"
