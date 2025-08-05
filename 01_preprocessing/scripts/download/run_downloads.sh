#!/bin/bash
set -euo pipefail

# 📁 Load config and utils
source "$(dirname "$0")/../../config/config.sh"
source "$(dirname "$0")/../../config/utils.sh"

# 🧱 Setup output folders
mkdir -p "$RAW_DIR" "$PREPROC_LOG_DIR"
JOBLIST="$METADATA_DIR/joblist.txt"
LOG_FILE="$PREPROC_LOG_DIR/download_summary.log"

# ⚙️ Determine max jobs dynamically (cap at 6)
TOTAL_CORES=$(nproc)
MAX_JOBS=$(( TOTAL_CORES / 2 ))
[[ $MAX_JOBS -lt 1 ]] && MAX_JOBS=1
[[ $MAX_JOBS -gt 6 ]] && MAX_JOBS=6

# 📝 Generate download task list
bash "$(dirname "$0")/make_joblist.sh"

# 🚀 Launch parallel downloads
echo "====================== 🌍 PARALLEL DOWNLOAD ======================="
echo "📄 Joblist        : $JOBLIST"
echo "🧵 Max Jobs       : $MAX_JOBS"
echo "📦 Output Folder  : $RAW_DIR"
echo "📓 Log File       : $LOG_FILE"
echo "==================================================================="

parallel -j "$MAX_JOBS" \
  --retries 3 \
  --tag --line-buffer \
  "$(dirname "$0")/download_one.sh" {1} {2} :::: "$JOBLIST" \
| tee "$LOG_FILE"
