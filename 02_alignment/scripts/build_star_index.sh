#!/bin/bash
# run_star_index.sh — STAR genome indexing with GENCODE v38 annotation

# =============================
# 1. Load Configuration
# =============================
# shellcheck disable=SC1091
source "$HOME/CrispAstro-Seq/config/config.sh"

# ===========================================
# Timestamp Block (KSA-standard logging)
# ===========================================

# Set timezone to local machine setting (e.g., KSA)
TZ=$(timedatectl | grep "Time zone" | awk '{print $3}')
export TZ
# Human-readable full timestamp (for console logs)
START_TIME_HUMAN=$(date +'%Y-%m-%d %H:%M:%S %Z (%:z)')

# Compact timestamp (for folder/file naming)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# =============================
# Thread Detection (Self-contained)
# =============================
TOTAL_THREADS=$(nproc)
THREADS=$TOTAL_THREADS
echo -e "\n🧠 Threads available: $THREADS"

# =============================
# 2. Indexing Parameters
# =============================
mkdir -p "$STAR_INDEX_DIR"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/star_index_${TIMESTAMP}.log"

echo -e "\n🚀 Starting STAR genome indexing..."
echo "🕒 Started at      : $START_TIME_HUMAN"
echo "Reference genome  : $GENOME_FASTA"
echo "Annotation GTF    : $GENOME_GTF"
echo "Index output dir  : $STAR_INDEX_DIR"
echo "Threads           : $THREADS"

# =============================
# 3. Run STAR Indexing
# =============================
STAR \
    --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$GENOME_GTF" \
    --sjdbOverhang 140 \
    >>"$LOG_FILE" 2>&1

# =============================
# 4. Finalize
# =============================

END_TIMESTAMP=$(date +%s)
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo -e "\n✅ STAR indexing complete!"
echo "🕒 Finished at     : $(date +'%Y-%m-%d %H:%M:%S %Z (%:z)')"
echo "⏱️  Total Runtime  : $RUNTIME_FMT"
echo "📄 Log file     : $LOG_FILE"
