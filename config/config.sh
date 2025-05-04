#!/bin/bash
# shellcheck disable=SC2034
# To Make sure it is executable run chmod +x config/config.sh
# ------------- Project Paths -------------

# Main project directory
PROJECT_DIR="$HOME/CrispAstro-Seq"

# Raw and processed data directories
RAW_DIR="$PROJECT_DIR/00_data/raw"
TRIMMED_DIR="$PROJECT_DIR/00_data/trimmed"
ALIGN_DIR="$PROJECT_DIR/00_data/aligned"
COUNTS_DIR="$PROJECT_DIR/00_data/counts"
REFERENCE_DIR="$PROJECT_DIR/00_data/reference"
METADATA_DIR="$PROJECT_DIR/00_data/metadata"

# QC and visualization results
QC_DIR="$PROJECT_DIR/01_qc/results"
VISUALIZATION_DIR="$PROJECT_DIR/07_visualization/figures"

# Configs and environment
CONFIG_DIR="$PROJECT_DIR/config"
ENV_DIR="$PROJECT_DIR/environment"

# ------------- Reference Files -------------

# Reference genome FASTA file
GENOME_FASTA="$REFERENCE_DIR/genome.fa"

# Reference genome annotation GTF file
GENOME_GTF="$REFERENCE_DIR/annotation.gtf"

# STAR genome index folder
STAR_INDEX_DIR="$REFERENCE_DIR/star_index"

# ------------- General Parameters -------------

# Number of threads for tools like STAR, FastQC, etc.
THREADS=12

# Other parameters (example: minimum read length after trimming)
MIN_LENGTH=30
