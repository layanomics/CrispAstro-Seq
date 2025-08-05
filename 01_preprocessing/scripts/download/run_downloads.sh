#!/bin/bash
set -euo pipefail

######################################################################
# 🧬 CrispAstro-Seq Pipeline
# Script        : run_downloads.sh
# Description   : Orchestrate parallel ENA FASTQ downloads using joblist
# Author        : Layan Essam
# Last Updated  : August 2025
######################################################################

# =============================
# 📦 Load Configuration & Utils
# =============================
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../../../config/config.sh"
source "$SCRIPT_DIR/../../../../config/utils.sh"

# =============================
# 🕒 Timestamp Initialization
# =============================
initialize_timestamps

# =============================
# 📁 Path Setup
# =============================
mkdir -p "$RAW_DIR" "$PREPROC_LOG_DIR"
JOBLIST="$METADATA_DIR/joblist.txt"
LOG_FILE="$PREPROC_LOG_DIR/download_summary_${TIMESTAMP}.log"

# =============================
# ⚙️ Determine max jobs dynamically (cap at 6)
# =============================
TOTAL_CORES=$(nproc)
MAX_JOBS=$(( TOTAL_CORES / 2 ))
[[ $MAX_JOBS -lt 1 ]] && MAX_JOBS=1
[[ $MAX_JOBS -gt 6 ]] && MAX_JOBS=6

# =============================
# 📝 Generate download task list
bash "$(dirname "$0")/make_joblist.sh"

#=============================
# 📝 Session Log Header
# =============================
{
  echo ""
  echo "====================== 🌍 PARALLEL DOWNLOAD ======================="
  echo "📄 Joblist              : $JOBLIST"
  echo "🧵 Max Parallel Jobs    : $MAX_JOBS"
  echo "📦 Output Directory     : $RAW_DIR"
  echo "🕒 Script Launched at   : $START_TIME_KSA | UTC: $START_TIME_UTC"
  echo "📓 Log File             : $LOG_FILE"
  echo "==================================================================="
} | tee -a "$LOG_FILE"
#=============================
# 🚀 Launch parallel downloads
# =============================
START_TIMESTAMP=$(date +%s)

parallel -j "$MAX_JOBS" \
  --retries 3 \
  --tag --line-buffer \
  "$(dirname "$0")/download_one.sh" {1} {2} :::: "$JOBLIST" \
| tee -a "$LOG_FILE"

# =============================
# ✅ Final Summary & Runtime
# =============================
SUCCESS_COUNT=$(find "$RAW_DIR" -name "*.fastq.gz" -exec bash -c 'gunzip -t "{}" &>/dev/null && echo OK' \; | wc -l)
EXPECTED_COUNT=$(( $(wc -l < "$JOBLIST") ))

{
  echo ""
  echo "====================== ✅ DOWNLOAD COMPLETED ======================"
  echo "📂 Output directory     : $RAW_DIR"
  report_runtime
  echo "🎯 Successful files     : $SUCCESS_COUNT of $EXPECTED_COUNT"

  if [[ "$SUCCESS_COUNT" -eq "$EXPECTED_COUNT" ]]; then
    echo "✅ Step completed successfully: run_downloads.sh"
  else
    echo "⚠️  Step completed with errors: $((EXPECTED_COUNT - SUCCESS_COUNT)) failed download(s)"
  fi

  echo "==================================================================="
} | tee -a "$LOG_FILE"





