#!/bin/bash

# ---------- CONFIG ----------
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../../config/config.sh"

# ---------- LOG SETUP ----------
LOG_DIR="$SCRIPT_DIR/../../logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/build_star_index_$(date +%Y%m%d_%H%M%S).log"

# ---------- START ----------
echo "🌟 Starting STAR index build at $(date)"
echo "Saving log to $LOG_FILE"
start_time=$(date +%s)

# ---------- CREATE INDEX ----------
mkdir -p "$STAR_INDEX_DIR"
STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$GENOME_GTF" \
    --sjdbOverhang 149 \
    &>"$LOG_FILE"

# ---------- END ----------
end_time=$(date +%s)
runtime=$((end_time - start_time))

# shellcheck disable=SC2181
if [ $? -eq 0 ]; then
    echo "✅ STAR genome index built successfully in $((runtime / 60)) minutes and $((runtime % 60)) seconds."
    echo "📁 Index stored in: $STAR_INDEX_DIR"
else
    echo "❌ STAR index build failed. Check log: $LOG_FILE"
fi
