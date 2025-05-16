#!/bin/bash
# config.sh — Central configuration for CrispAstro-Seq project
# Used by all preprocessing and analysis scripts for consistency
# shellcheck disable=SC2034

# ============================================================
# 1. Define Project Directory Structure
# ============================================================

# Main project directory (root of all paths)
PROJECT_DIR="$HOME/CrispAstro-Seq"

# Raw and processed FASTQ data
RAW_DIR="$PROJECT_DIR/00_data/raw"
TRIMMED_DIR="$PROJECT_DIR/00_data/trimmed"
ALIGN_DIR="$PROJECT_DIR/00_data/aligned"
COUNTS_DIR="$PROJECT_DIR/00_data/counts"
REFERENCE_DIR="$PROJECT_DIR/00_data/reference"
METADATA_DIR="$PROJECT_DIR/00_data/metadata"

# Quality control and downstream visualization
QC_DIR="$PROJECT_DIR/01_preprocessing/results"
VISUALIZATION_DIR="$PROJECT_DIR/07_visualization/figures"

# Configs, conda environments, and logs
CONFIG_DIR="$PROJECT_DIR/config"
ENV_DIR="$PROJECT_DIR/environment"
LOG_DIR="$PROJECT_DIR/logs"

# ============================================================
# 2. Reference Files
# ============================================================

# Genome FASTA and annotation GTF used in STAR, featureCounts, etc.
GENOME_FASTA="$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="$REFERENCE_DIR/gencode.v38.annotation.gtf"

# Directory to hold STAR index files
STAR_INDEX_DIR="$REFERENCE_DIR/star_index"

# ============================================================
# 3. Threading & Parallelism Logic
# ============================================================

# Detect how many CPU threads are available
TOTAL_THREADS=$(nproc)

# Detect environment type (assumes battery = laptop)
MACHINE_TYPE=$([ -d "/sys/class/power_supply/BAT0" ] && echo "laptop" || echo "workstation_or_vm")

# Assign readable machine label
if [[ "$MACHINE_TYPE" == "laptop" ]]; then
    MACHINE_LABEL="Personal Laptop"
    MAX_JOBS=$((TOTAL_THREADS * 60 / 100)) # Limit jobs to preserve usability
else
    MACHINE_LABEL="Virtual Machine / Cloud Environment"
    MAX_JOBS=$TOTAL_THREADS # Use full power
fi

# Ensure at least one job is allowed
[ "$MAX_JOBS" -lt 1 ] && MAX_JOBS=1

# Distribute threads across jobs (auto-adjusted)
THREADS_PER_JOB=$((TOTAL_THREADS / MAX_JOBS))

# Export THREADS for downstream tools (e.g., STAR, R, featureCounts)
export THREADS=$TOTAL_THREADS

# ============================================================
# 4. Sample Counting for Reporting
# ============================================================

# Count number of R1 FASTQ files (paired-end assumption)
# shellcheck disable=SC2012
NUM_SAMPLES=$(ls "$RAW_DIR"/*_1.fastq.gz 2>/dev/null | wc -l)

# Estimate number of parallel jobs that will actually run
if [ "$NUM_SAMPLES" -lt "$MAX_JOBS" ]; then
    JOBS_IN_PARALLEL=$NUM_SAMPLES
else
    JOBS_IN_PARALLEL=$MAX_JOBS
fi

# ============================================================
# 5. Print Configuration Summary
# ============================================================

echo -e "\n✅ Config Loaded:
- MACHINE TYPE             : $MACHINE_LABEL
- TOTAL CPU THREADS        : $TOTAL_THREADS
- MAX PARALLEL JOBS        : $MAX_JOBS
- THREADS PER JOB          : $THREADS_PER_JOB
- PAIRED-END SAMPLES FOUND : $NUM_SAMPLES
- JOBS TO RUN IN PARALLEL  : $JOBS_IN_PARALLEL
- TOTAL SAMPLES TO PROCESS : $NUM_SAMPLES
"
