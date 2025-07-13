#!/bin/bash
# star_align.sh — STAR alignment for trimmed FASTQ reads

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

echo -e "\n✅ Config Loaded:"
echo "🕒 Started at      : $START_TIME_HUMAN"
echo "- MACHINE TYPE             : Virtual Machine / Cloud Environment"
echo "- TOTAL CPU THREADS        : $THREADS"
echo "- MAX PARALLEL JOBS        : $MAX_JOBS"
echo "- THREADS PER JOB          : $THREADS_PER_JOB"
echo "- PAIRED-END SAMPLES FOUND : $(find "$TRIMMED_DIR" -name "*_clean_R1.fastq.gz" | wc -l)"
echo "- JOBS TO RUN IN PARALLEL  : $MAX_JOBS"
echo "- TOTAL SAMPLES TO PROCESS : $(find "$TRIMMED_DIR" -name "*_clean_R1.fastq.gz" | wc -l)"

echo -e "\n⚙️ STAR override:"
echo "  - MAX_JOBS (STAR Align)  : $MAX_JOBS"
echo "  - THREADS_PER_JOB        : $THREADS_PER_JOB"

mkdir -p "$ALIGN_DIR" "$LOG_DIR"

for r1 in "$TRIMMED_DIR"/*/*_clean_R1.fastq.gz; do
    sample=$(basename "$r1" | cut -d'_' -f1)
    r2="$TRIMMED_DIR/$sample/${sample}_clean_R2.fastq.gz"

    if [[ ! -f "$r2" || ! -s "$r1" || ! -s "$r2" ]]; then
        echo "❌ Skipping $sample — missing or empty reads" | tee -a "$LOG_DIR/star_align_errors_${TIMESTAMP}.log"

        continue
    fi

    echo "[$(date +'%F %T')] 🚀 STAR alignment for $sample..."

    SAMPLE_OUT_DIR="$ALIGN_DIR/$sample"
    mkdir -p "$SAMPLE_OUT_DIR"
    LOG_FILE="$LOG_DIR/star_align_${sample}_${TIMESTAMP}.log"

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
        --limitBAMsortRAM 30000000000 \
        >>"$LOG_FILE" 2>&1

    # shellcheck disable=SC2181
    if [ $? -ne 0 ]; then
        echo "❌ STAR failed for $sample. Check log: $LOG_FILE"
    else
        echo "✅ STAR finished successfully for $sample"
    fi

done

# =============================
# 3. Finalize
# =============================

END_TIMESTAMP=$(date +%s)
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo -e "\n✅ All STAR alignments complete!"
echo "🕒 Finished at     : $(date +'%Y-%m-%d %H:%M:%S %Z (%:z)')"
echo "⏱️  Total Runtime  : $RUNTIME_FMT"
echo "📄 Log file     : $LOG_FILE"
