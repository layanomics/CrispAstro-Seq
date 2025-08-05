#!/bin/bash
set -euo pipefail
######################################################################
# 🧬 CrispAstro-Seq Pipeline
# Script        : run_multiqc.sh
# Description   : Run MultiQC on raw and/or trimmed FastQC results
# Author        : Layan Essam
# Tool Version  : MultiQC v1.14
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
check_tool_version multiqc

# =============================
# 🏷️ Handle Optional Arguments
# =============================
TYPE=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --type)
      TYPE="$2"
      shift 2
      ;;
    *)
      echo "❌ Unknown option: $1"
      echo "✅ Usage: bash run_multiqc.sh --type raw|trimmed"
      exit 1
      ;;
  esac
done

if [[ -z "$TYPE" ]]; then
  echo "❌ ERROR: --type is required. Please specify --type raw or --type trimmed"
  exit 1
fi

# =============================
# 🕒 Timestamp Initialization
# =============================
initialize_timestamps

# =============================
# 📁 Input & Output Paths
# =============================
RAW_INPUT="$FASTQC_RAW_DIR"
TRIMMED_INPUT="$FASTQC_TRIMMED_DIR"

RAW_OUTPUT="$MULTIQC_RAW_DIR"
TRIMMED_OUTPUT="$MULTIQC_TRIMMED_DIR"

mkdir -p "$RAW_OUTPUT" "$TRIMMED_OUTPUT" "$PREPROC_RESULTS_DIR/multiqc"

# Select log file destination
if [[ "$TYPE" == "raw" ]]; then
  LOG_FILE="$RAW_OUTPUT/multiqc_raw_${TIMESTAMP}.log"
elif [[ "$TYPE" == "trimmed" ]]; then
  LOG_FILE="$TRIMMED_OUTPUT/multiqc_trimmed_${TIMESTAMP}.log"
else
  LOG_FILE="$PREPROC_RESULTS_DIR/multiqc/multiqc_all_${TIMESTAMP}.log"
fi

# =============================
# 📝 Session Log Header
# =============================
{
  echo ""
  echo "====================== 🧪 MULTIQC REPORTING ======================="
  echo "⚙️  MultiQC CONFIGURATION"
  echo "🔧 Version               : $TOOL_VERSION"
  echo "📦 Report Type           : $TYPE"
  echo "🕒 Script launched at    : $START_TIME_KSA | UTC: $START_TIME_UTC"
  echo ""
  echo "======================📂 File Structure ============================"
  echo "📂 Input Subdirectory    :"
  [[ "$TYPE" == "raw" || "$TYPE" == "all" ]] && echo "   └─ $RAW_INPUT"
  [[ "$TYPE" == "trimmed" || "$TYPE" == "all" ]] && echo "   └─ $TRIMMED_INPUT"
  echo ""
  echo "📁 MultiQC Output Folder(s):"
  [[ "$TYPE" == "raw" || "$TYPE" == "all" ]] && echo "   └─ $RAW_OUTPUT"
  [[ "$TYPE" == "trimmed" || "$TYPE" == "all" ]] && echo "   └─ $TRIMMED_OUTPUT"
  echo ""
  echo "📄 Logging to             : $LOG_FILE"
  echo "==================================================================="
} | tee -a "$LOG_FILE"


# =============================
# 🚀 Run MultiQC
# =============================
START_TIMESTAMP=$(date +%s)

if [[ "$TYPE" == "raw" || "$TYPE" == "all" ]]; then
    echo -e "\n📁 Running MultiQC for RAW..." | tee -a "$LOG_FILE"
    multiqc "$RAW_INPUT" -o "$RAW_OUTPUT" 2>&1 | tee -a "$LOG_FILE"
fi

if [[ "$TYPE" == "trimmed" || "$TYPE" == "all" ]]; then
    echo -e "\n📁 Running MultiQC for TRIMMED..." | tee -a "$LOG_FILE"
    multiqc "$TRIMMED_INPUT" -o "$TRIMMED_OUTPUT" 2>&1 | tee -a "$LOG_FILE"
fi


# =============================
# ✅ Final Summary
# =============================
log_runtime

{
  echo ""
  echo "==================== ✅ MULTIQC COMPLETED ====================="
  [[ "$TYPE" == "raw" || "$TYPE" == "all" ]] && echo "📄 Raw report saved to     : $RAW_OUTPUT/multiqc_report.html"
  [[ "$TYPE" == "trimmed" || "$TYPE" == "all" ]] && echo "📄 Trimmed report saved to : $TRIMMED_OUTPUT/multiqc_report.html"
  print_runtime_summary
  echo "✅ Step completed successfully: run_multiqc"
  echo "==================================================================="
} | tee -a "$LOG_FILE"
