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
