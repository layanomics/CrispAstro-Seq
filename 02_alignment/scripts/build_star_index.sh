#!/bin/bash
# run_star_index.sh — STAR genome indexing with GENCODE v38 annotation

# =============================
# 1. Load Configuration
# =============================
# shellcheck disable=SC1090
source ~/CrispAstro-Seq/config/config.sh

# Start time tracking
START_TIME=$(date +%s)

# =============================
# 2. Indexing Parameters
# =============================
mkdir -p "$STAR_INDEX_DIR"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/star_index_$(date +%Y%m%d_%H%M%S).log"

echo -e "\n🚀 Starting STAR genome indexing..."
echo "Reference genome : $GENOME_FASTA"
echo "Annotation GTF   : $GENOME_GTF"
echo "Index output dir : $STAR_INDEX_DIR"
echo "Threads          : $THREADS"

# =============================
# 3. Run STAR Indexing
# =============================
STAR \
    --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$GENOME_GTF" \
    --sjdbOverhang 100 \
    >>"$LOG_FILE" 2>&1

# =============================
# 4. Finalize
# =============================
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo -e "\n✅ STAR indexing complete!"
echo "Total Runtime: $RUNTIME_FMT"
echo "Log: $LOG_FILE"
