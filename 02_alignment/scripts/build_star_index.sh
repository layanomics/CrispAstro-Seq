#!/bin/bash
# run_star_index.sh — STAR genome indexing with GENCODE v38 annotation

# =============================
# 1. Load Configuration
# =============================
# shellcheck disable=SC1090
source ~/CrispAstro-Seq/config/config.sh

# Optional: Force date to use local timezone (e.g., KSA)
# shellcheck disable=SC2155
export TZ=$(timedatectl | grep "Time zone" | awk '{print $3}')

# Start time tracking
START_TIMESTAMP=$(date +%s)
START_TIME_HUMAN=$(date +'%Y-%m-%d %H:%M:%S %Z (%:z)')

# =============================
# 2. Indexing Parameters
# =============================
mkdir -p "$STAR_INDEX_DIR"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/star_index_$(date +%Y%m%d_%H%M%S).log"

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
echo "🕒 Finished at : $(date +'%Y-%m-%d %H:%M:%S %Z (%:z)')"
echo "⏱️  Total Runtime: $RUNTIME_FMT"
echo "📄 Log file     : $LOG_FILE"
