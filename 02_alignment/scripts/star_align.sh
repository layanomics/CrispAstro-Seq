#!/bin/bash

# ---------- CONFIG ----------
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
source "$SCRIPT_DIR/../../config/config.sh"

# Detect if running on the VM and override THREADS accordingly
if [[ "$(hostname)" == "star-indexer" ]]; then
    THREADS=$(nproc) # Use all available VM threads
    echo "⭐ Running on VM: Using $THREADS threads."
else
    echo "💻 Running on local machine: Using $THREADS threads from config."
fi

# ---------- LOG SETUP ----------
LOG_DIR="$SCRIPT_DIR/../../logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/star_align_$(date +%Y%m%d_%H%M%S).log"

echo "🚀 Starting STAR Alignment at $(date)"
echo "Logging to: $LOG_FILE"
start_time=$(date +%s)

# ---------- ALIGNMENT ----------
mkdir -p "$ALIGN_DIR"

for R1_FILE in "$TRIMMED_DIR"/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1_FILE" | sed 's/_R1\.fastq\.gz//')
    R2_FILE="$TRIMMED_DIR/${SAMPLE}_R2.fastq.gz"

    echo "🧬 Processing Sample: $SAMPLE" | tee -a "$LOG_FILE"

    STAR --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$R1_FILE" "$R2_FILE" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$ALIGN_DIR/${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        &>>"$LOG_FILE"
done

# ---------- END ----------
end_time=$(date +%s)
runtime=$((end_time - start_time))

# shellcheck disable=SC2181
if [ $? -eq 0 ]; then
    echo "✅ STAR Alignment completed successfully in $((runtime / 60)) minutes and $((runtime % 60)) seconds."
else
    echo "❌ STAR Alignment failed. Check log: $LOG_FILE"
fi
