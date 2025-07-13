#!/bin/bash
# run_fastqc_trimmed.sh — Run FastQC on trimmed paired-end reads

# =============================
# 1. Load Config & Set Up
# =============================
# shellcheck disable=SC1091
source "$HOME/CrispAstro-Seq/config/config.sh"

FASTQC_TYPE="trimmed"
FASTQC_INPUT_DIR="$TRIMMED_DIR"
FASTQC_OUT="$QC_DIR/fastqc_$FASTQC_TYPE"
mkdir -p "$FASTQC_OUT"
LOG_FILE="$FASTQC_OUT/fastqc_${FASTQC_TYPE}_$(date +%Y%m%d_%H%M%S).log"

# =============================
# 2. Parallelism Logic
# =============================
TOTAL_THREADS=$(nproc)
NUM_SAMPLES=$(find "$FASTQC_INPUT_DIR" -name "*_clean_R1.fastq.gz" | wc -l)

# Set MAX_JOBS = min(samples, threads)
if [ "$NUM_SAMPLES" -lt "$TOTAL_THREADS" ]; then
    MAX_JOBS=$NUM_SAMPLES
else
    MAX_JOBS=$TOTAL_THREADS
fi

[ "$MAX_JOBS" -lt 1 ] && MAX_JOBS=1
THREADS_PER_JOB=1

echo -e "\n⚙️ FastQC CONFIG (TRIMMED):"
echo "  - TOTAL THREADS        : $TOTAL_THREADS"
echo "  - TOTAL SAMPLES        : $NUM_SAMPLES"
echo "  - MAX PARALLEL JOBS    : $MAX_JOBS"
echo "  - THREADS PER JOB      : $THREADS_PER_JOB"

# =============================
# 3. Run FastQC in Parallel
# =============================
START_TIME=$(date +%s)
CURRENT_JOBS=0

echo -e "\n🔍 Running FastQC on trimmed reads...\n"

for r1 in "$FASTQC_INPUT_DIR"/*/*_clean_R1.fastq.gz; do
    sample=$(basename "$r1" _clean_R1.fastq.gz)
    r2="$FASTQC_INPUT_DIR/$sample/${sample}_clean_R2.fastq.gz"
    mkdir -p "$FASTQC_OUT/$sample"

    echo "[$(date +'%F %T')] FastQC for $sample (trimmed)" | tee -a "$LOG_FILE"
    fastqc -t "$THREADS_PER_JOB" "$r1" "$r2" \
        --outdir="$FASTQC_OUT/$sample" >>"$LOG_FILE" 2>&1 &

    ((CURRENT_JOBS++))
    [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ] && wait && CURRENT_JOBS=0
done
wait

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')

echo -e "\n✅ FastQC ($FASTQC_TYPE) finished in $RUNTIME_FMT"
