#!/bin/bash
set -euo pipefail

#############################################################################
# 🌐 CrispAstro-Seq Pipeline 
# script: download_fastq_from_ena.sh
# Description: Download paired-end FASTQ files from ENA with retry & integrity check
# Author: Layan Essam
# Last Updated: August 2025
#############################################################################

# ==============================================
# 📌 Load shared configuration and utilities
# ==============================================
source "$(dirname "$0")/../../config/config.sh"
source "$(dirname "$0")/../../config/utils.sh"

# =============================
# 🕒 Timestamp Initialization
# =============================
initialize_timestamps

# =============================
# 📁 Setup Input/Output & Log Paths
# =============================
SRA_LIST="$METADATA_DIR/sra_ids.txt"
mkdir -p "$RAW_DIR" "$PREPROC_LOG_DIR"
LOG_FILE="$PREPROC_LOG_DIR/download_${TIMESTAMP}.log"

# =============================
# 📑 SRA ID List Validation
# =============================
mapfile -t SRA_IDS < "$SRA_LIST"

if [[ ${#SRA_IDS[@]} -eq 0 ]]; then
  echo "❌ ERROR: SRA list is empty: $SRA_LIST" | tee -a "$LOG_FILE"
  exit 1
fi

# ==============================
# 📝 Session Log Header
# ==============================

{
  echo ""
  echo "====================== 🌍 FASTQ DOWNLOAD ======================="
  echo "📄 SRA IDs File       : $SRA_LIST"
  echo "📄 Samples to Download: ${#SRA_IDS[@]}"
  echo "📂 Output Directory   : $RAW_DIR"
  echo "🖥️ Log File           : $LOG_FILE"
  echo "🕒 Script Launched at : $START_TIME_KSA | UTC: $START_TIME_UTC"
  echo "================================================================"
  echo ""
} | tee -a "$LOG_FILE"

# ============================================
# ⚙️ Parallelism & Sample Thread Logic Section
# ============================================
TOTAL_CORES=$(nproc)
MAX_PARALLEL=$((TOTAL_CORES / 2))
[[ $MAX_PARALLEL -lt 1 ]] && MAX_PARALLEL=1
[[ $MAX_PARALLEL -gt 6 ]] && MAX_PARALLEL=6
MAX_JOBS=$MAX_PARALLEL

# ===========================================
# 🔁 Dynamic Parallel Download with retry and integrity check
# ===========================================
START_TIMESTAMP=$(date +%s)

JOB_COUNT=0
SUCCESS_COUNT=0

for SRA in "${SRA_IDS[@]}"; do
  echo -e "\n🧬 $SRA" | tee -a "$LOG_FILE"

  # 🔄 ENA block detection using API
  BLOCK_FOUND=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRA&result=read_run&fields=fastq_ftp&format=tsv" |
                grep -v fastq_ftp | head -n1 | awk -F '/' '{print $(NF-2)}')

  if [[ -z "$BLOCK_FOUND" ]]; then
    echo "🚫 ENA files not found for: $SRA" | tee -a "$LOG_FILE"
    continue
  fi

  for READ in 1 2; do
    (
      FILE="${SRA}_${READ}.fastq.gz"
      DEST="$RAW_DIR/$FILE"
      URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/${SRA:0:6}/${BLOCK_FOUND}/$SRA/$FILE"

      if ! [[ -f "$DEST" ]] || ! gunzip -t "$DEST" &>/dev/null; then
        ATTEMPT=1
        MAX_TRIES=3
        SUCCESS=0

        while [[ $ATTEMPT -le $MAX_TRIES ]]; do
          echo "🌐 Attempt $ATTEMPT: Downloading $FILE..." | tee -a "$LOG_FILE"

          if wget -q --show-progress --progress=bar:force -c "$URL" -O "$DEST"; then
            echo "📦 Downloaded: $FILE" | tee -a "$LOG_FILE"
            echo "🔍 Verifying..." | tee -a "$LOG_FILE"

            if gunzip -t "$DEST" &>/dev/null; then
              echo "✅ Verified: $FILE" | tee -a "$LOG_FILE"
              SUCCESS=1
              ((SUCCESS_COUNT++))
              break
            else
              echo "❌ Corrupted: $FILE — retrying..." | tee -a "$LOG_FILE"
              rm -f "$DEST"
            fi

          else
            echo "⚠️  wget failed: $FILE" | tee -a "$LOG_FILE"
          fi

          ((ATTEMPT++))
          sleep 2
        done

        if [[ $SUCCESS -ne 1 ]]; then
          echo "🚫 Failed after $MAX_TRIES attempts: $FILE" | tee -a "$LOG_FILE"
        fi

      else
        echo "✅ Already downloaded and verified: $FILE" | tee -a "$LOG_FILE"
      fi
    ) &

    ((JOB_COUNT++))

    # ✅ Real fix to prevent early prompt — handles both batch and final job
    if (( JOB_COUNT % MAX_JOBS == 0 )) || (( JOB_COUNT == ${#SRA_IDS[@]} * 2 )); then
      wait
    fi

  done  # End of READ loop

done  # End of SRA loop

wait  # 🔚 Final batch to ensure all background jobs finish


# =============================
# ✅ Final Summary & Runtime
# =============================

{
  echo ""
  echo "====================== ✅ DOWNLOAD COMPLETED ======================"
  echo "🎯 Successfully verified : $SUCCESS_COUNT of $(( ${#SRA_IDS[@]} * 2 )) FASTQ files"
  echo "📂 FASTQ files directory : $RAW_DIR"
  report_runtime
  echo ""

  echo "📦 Downloaded FASTQ Files:"
  printf "%-30s %10s\n" "Filename" "Size (MB)"
  echo "-----------------------------------------"
  for FILE in "$RAW_DIR"/*.fastq.gz; do
    SIZE=$(du -m "$FILE" | cut -f1)
    printf "%-30s %10s\n" "$(basename "$FILE")" "$SIZE"
  done
  echo ""


  if [[ "$SUCCESS_COUNT" -eq $(( ${#SRA_IDS[@]} * 2 )) ]]; then
    echo "✅ Step completed successfully: download_fastq_from_ena.sh"
  else
    echo "⚠️  Step completed with errors: $(( (${#SRA_IDS[@]} * 2) - SUCCESS_COUNT )) file(s) failed"
  fi
  echo "==================================================================="
} | tee -a "$LOG_FILE"













