#!/bin/bash
set -euo pipefail
######################################################################
# 🧬 CrispAstro-Seq Pipeline
# Script        : build_star_index.sh
# Description   : Build STAR genome index (GRCh38.p12) using GENCODE v38 GTF annotation 
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
# 📁 Path & Directory Setup
# =============================
mkdir -p "$STAR_INDEX_DIR" "$ALIGN_LOG_DIR"
LOG_FILE="$ALIGN_LOG_DIR/star_index_${TIMESTAMP}.log"

# =============================
# 🧠 Thread Detection (Self-contained)
# =============================
TOTAL_THREADS=$(nproc); THREADS=$TOTAL_THREADS

# =============================
#📝 Session Log Header
# =============================
{
  echo ""
  echo "==================== 🧬 STAR GENOME INDEXING ===================="
  echo "⚙️  STAR CONFIGURATION"
  echo "🔢 STAR Version         : $TOOL_VERSION"
  echo "🧠 Threads              : $THREADS"
  echo ""
  echo "📋 STAR genome indexing setup:"
  echo "🕒 Script launched at   : $START_TIME_KSA | UTC: $START_TIME_UTC"
  echo "🧬 Reference genome     : $GENOME_FASTA"
  echo "📘 Annotation GTF       : $GENOME_GTF"
  echo "📁 Output directory     : $STAR_INDEX_DIR"
  echo "📄 Logging to           : $LOG_FILE"
  echo "==============================================================="
} | tee -a "$LOG_FILE"

# =============================
# 🚀 Run STAR Indexing
# =============================
echo  "🚀 Starting STAR genome indexing..."

START_TIMESTAMP=$(date +%s)  
STAR \
    --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$GENOME_GTF" \
    >>"$LOG_FILE" 2>&1

# =============================
# ✅ Index Output Validation
# =============================
if compgen -G "$STAR_INDEX_DIR/*" > /dev/null; then
  echo "✅ Index files successfully generated." | tee -a "$LOG_FILE"
else
  echo "❌ Indexing failed: No files found in $STAR_INDEX_DIR" | tee -a "$LOG_FILE"
  exit 1
fi

# =============================
# Finalize and Runtime Summary
# =============================
log_runtime

{
  echo ""
  echo "=================== ✅ STAR INDEX COMPLETE ==================="
  echo "✅ STAR indexing completed successfully!"
  print_runtime_summary
  echo "📄 Log file saved to : $LOG_FILE"
  echo "=============================================================="
} | tee -a "$LOG_FILE"
