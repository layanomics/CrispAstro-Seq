#!/bin/bash
set -euo pipefail

#############################################################################
# 🧬 CrispAstro-Seq Pipeline 
# script: extract_metadata_and_ids.sh
# Description: Extracts metadata.csv and sra_ids.txt from SraRunTable.csv
# Author: Layan Essam
# Last Updated  : August 2025
#############################################################################

# ==============================================
# 📌 Load shared configuration and utilities
# ==============================================
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../../../.. && pwd)"
source "$PROJECT_DIR/config/config.sh"
source "$PROJECT_DIR/config/utils.sh"

# =============================
# 🕒 Timestamp Initialization
# =============================
initialize_timestamps

# =============================
# 📁 Paths Setup
# =============================
SRA_TABLE="$METADATA_DIR/SraRunTable.csv"
METADATA_OUT="$METADATA_DIR/metadata.csv"
SRA_IDS_OUT="$METADATA_DIR/sra_ids.txt"

# Ensure output directory exists
mkdir -p "$METADATA_DIR"

# Check input file exists
if [[ ! -f "$SRA_TABLE" ]]; then
    echo "❌ ERROR: SraRunTable.csv not found at: $SRA_TABLE"
    exit 1
fi

# =============================
# 📝 Session Log Header
# =============================
echo ""
echo "================== 🧬 METADATA EXTRACTION ==================="
echo "📥 Source Table    : $SRA_TABLE"
echo "📄 Metadata Output : $METADATA_OUT"
echo "🧾 SRA IDs Output  : $SRA_IDS_OUT"
echo "🕒 Started at      : $START_TIME_KSA | UTC: $START_TIME_UTC"
echo "============================================================="
echo ""

# =============================
# 🔍 Extract Metadata
# =============================

START_TIMESTAMP=$(date +%s)

# Extract metadata.csv: sample and condition
csvcut -c "Run","genotype" "$SRA_TABLE" | sed '1s/.*/sample,condition/' > "$METADATA_OUT"

# Extract SRA IDs only (no header)
csvcut -c "Run" "$SRA_TABLE" | tail -n +2 > "$SRA_IDS_OUT"

# =============================
# ✅ Final Summary & Runtime
# =============================
log_runtime
print_runtime_summary

echo "✅ metadata.csv and sra_ids.txt created under: $METADATA_DIR"
echo "============================================================="
