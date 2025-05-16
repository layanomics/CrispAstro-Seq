#!/bin/bash
# run_star_align.sh — Aligns trimmed FASTQ reads using STAR

# =============================
# 1. Load Configuration
# =============================
# shellcheck disable=SC1090
source ~/CrispAstro-Seq/config/config.sh

# Start time tracking
START_TIME=$(date +%s)

# Override parallel logic for STAR alignment
MAX_JOBS=1
THREADS_PER_JOB=$THREADS

echo -e "\n⚙️ STAR override:"
echo "  - MAX_JOBS (STAR Align)  : $MAX_JOBS"
echo "  - THREADS_PER_JOB        : $THREADS_PER_JOB"

mkdir -p "$ALIGN_DIR" "$LOG_DIR"

# =============================
# 2. Sample Loop
# =============================

CURRENT_JOBS=0
for r1 in "$TRIMMED_DIR"/*/*_clean_R1.fastq.gz; do
    sample=$(basename "$r1" | cut -d'_' -f1)
    r2="$TRIMMED_DIR/$sample/${sample}_clean_R2.fastq.gz"

    if [[ ! -f "$r2" || ! -s "$r1" || ! -s "$r2" ]]; then
        echo "❌ Skipping $sample — missing or empty reads" | tee -a "$LOG_DIR/star_align_errors.log"
        continue
    fi

    echo "[$(date +'%F %T')] STAR alignment for $sample..."

    SAMPLE_OUT_DIR="$ALIGN_DIR/$sample"
    mkdir -p "$SAMPLE_OUT_DIR"
    LOG_FILE="$LOG_DIR/star_align_${sample}.log"

    STAR \
        --runThreadN "$THREADS_PER_JOB" \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$SAMPLE_OUT_DIR/${sample}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode GeneCounts \
        --twopassMode Basic \
        >>"$LOG_FILE" 2>&1 &

    ((CURRENT_JOBS++))
    [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ] && wait && CURRENT_JOBS=0
done

wait

# =============================
# 3. Finalize
# =============================
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo -e "\n✅ All STAR alignments complete!"
echo "Total Runtime: $RUNTIME_FMT"
