#!/bin/bash
# ===============================================================
# 📦 CrispAstro-Seq Central Configuration File
# Author: Layan Essam
# Last Updated: July 2025
# Purpose: Central configuration for CrispAstro-Seq pipeline-Define consistent paths and references across scripts 
# ===============================================================

# ========================================
# 📁 Main Project Directory
# ========================================
CONFIG_PATH="$(realpath "${BASH_SOURCE[0]}")"
PROJECT_DIR="$(realpath "$(dirname "$CONFIG_PATH")/..")"
# ========================================
# 📥 STEP 1: FASTQC on RAW READS
# ========================================
RAW_DIR="$PROJECT_DIR/00_data/raw"
REFERENCE_DIR="$PROJECT_DIR/00_data/reference"
PREPROC_RESULTS_DIR="$PROJECT_DIR/01_preprocessing/results"
FASTQC_RAW_DIR="$PREPROC_RESULTS_DIR/fastqc/raw"
PREPROC_LOG_DIR="$PROJECT_DIR/01_preprocessing/logs"

# ========================================
# ✂️ STEP 2: CUTADAPT TRIMMING
# ========================================
TRIMMED_DIR="$PROJECT_DIR/00_data/trimmed"
CUTADAPT_LOG_DIR="$PREPROC_RESULTS_DIR/cutadapt"

# ========================================
# 📊 STEP 3: FASTQC on TRIMMED READS
# ========================================
FASTQC_TRIMMED_DIR="$PREPROC_RESULTS_DIR/fastqc/trimmed"


# ========================================
# 🧪 STEP 4: MULTIQC REPORTING
# ========================================
MULTIQC_RAW_DIR="$PREPROC_RESULTS_DIR/multiqc/raw" 
MULTIQC_TRIMMED_DIR="$PREPROC_RESULTS_DIR/multiqc/trimmed"

# ========================================
# 🧬 STEP 5: STAR GENOME INDEXING
# ========================================
GENOME_FASTA="$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="$REFERENCE_DIR/gencode.v38.annotation.gtf"
STAR_INDEX_DIR="$REFERENCE_DIR/star_index"

# ========================================
# 🔭 STEP 6: STAR ALIGNMENT
# ========================================
ALIGN_DIR="$PROJECT_DIR/00_data/aligned"
ALIGN_LOG_DIR="$PROJECT_DIR/02_alignment/logs"

# ========================================
# 🧾 EXTRA: METADATA & FEATURE COUNTS
# ========================================
COUNTS_DIR="$PROJECT_DIR/00_data/counts"
METADATA_DIR="$PROJECT_DIR/00_data/metadata"

# ========================================

#!/bin/bash
# ===============================================================
# 📦 CrispAstro-Seq Central Configuration File
# Author: Layan Essam
# Last Updated: July 2025
# Purpose: Central configuration for CrispAstro-Seq pipeline-Define consistent paths and references across scripts 
# ===============================================================

# ========================================
# 📁 Main Project Directory
# ========================================
CONFIG_PATH="$(realpath "${BASH_SOURCE[0]}")"
PROJECT_DIR="$(realpath "$(dirname "$CONFIG_PATH")/..")"
# ========================================
# 📥 STEP 1: FASTQC on RAW READS
# ========================================
RAW_DIR="$PROJECT_DIR/00_data/raw"
REFERENCE_DIR="$PROJECT_DIR/00_data/reference"
PREPROC_RESULTS_DIR="$PROJECT_DIR/01_preprocessing/results"
FASTQC_RAW_DIR="$PREPROC_RESULTS_DIR/fastqc/raw"
PREPROC_LOG_DIR="$PROJECT_DIR/01_preprocessing/logs"

# ========================================
# ✂️ STEP 2: CUTADAPT TRIMMING
# ========================================
TRIMMED_DIR="$PROJECT_DIR/00_data/trimmed"
CUTADAPT_LOG_DIR="$PREPROC_RESULTS_DIR/cutadapt"

# ========================================
# 📊 STEP 3: FASTQC on TRIMMED READS
# ========================================
FASTQC_TRIMMED_DIR="$PREPROC_RESULTS_DIR/fastqc/trimmed"


# ========================================
# 🧪 STEP 4: MULTIQC REPORTING
# ========================================
MULTIQC_RAW_DIR="$PREPROC_RESULTS_DIR/multiqc/raw" 
MULTIQC_TRIMMED_DIR="$PREPROC_RESULTS_DIR/multiqc/trimmed"

# ========================================
# 🧬 STEP 5: STAR GENOME INDEXING
# ========================================
GENOME_FASTA="$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa"
GENOME_GTF="$REFERENCE_DIR/gencode.v38.annotation.gtf"
STAR_INDEX_DIR="$REFERENCE_DIR/star_index"

# ========================================
# 🔭 STEP 6: STAR ALIGNMENT
# ========================================
ALIGN_DIR="$PROJECT_DIR/00_data/aligned"
ALIGN_LOG_DIR="$PROJECT_DIR/02_alignment/logs"

# ========================================
# 🧾 EXTRA: METADATA & FEATURE COUNTS
# ========================================
COUNTS_DIR="$PROJECT_DIR/00_data/counts"
METADATA_DIR="$PROJECT_DIR/00_data/metadata"

# ========================================
