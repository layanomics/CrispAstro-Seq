#!/bin/bash
# merge_star_counts.sh — Final version for DESeq2: filters N_ rows, strips .version from GeneID

set -euo pipefail

# Paths
ALIGN_DIR=~/CrispAstro-Seq/00_data/aligned
OUT_MATRIX=~/CrispAstro-Seq/03_quantification/results/counts.tsv

mkdir -p "$(dirname "$OUT_MATRIX")"

echo "🧬 Generating DESeq2-ready count matrix (column 4: reverse-stranded)"
echo "Output file: $OUT_MATRIX"

HEADER="GeneID"
paste_files=()

for file in "$ALIGN_DIR"/SRR*/SRR*_ReadsPerGene.out.tab; do
    sample=$(basename "$file" | cut -d_ -f1)
    echo "➕ Processing sample: $sample"

    if [ ${#paste_files[@]} -eq 0 ]; then
        # First file: keep gene ID and counts, remove N_* lines, strip .version
        awk '$1 !~ /^N_/ { sub(/\..*/, "", $1); print $1 "\t" $4 }' "$file" >"tmp_${sample}.tsv"
    else
        # Others: just column 4, matching order
        awk '$1 !~ /^N_/ { print $4 }' "$file" >"tmp_${sample}.tsv"
    fi

    paste_files+=("tmp_${sample}.tsv")
    HEADER+="\t$sample"
done

paste "${paste_files[@]}" >tmp_combined.tsv
(echo -e "$HEADER" && cat tmp_combined.tsv) >"$OUT_MATRIX"

rm tmp_*.tsv tmp_combined.tsv

echo "✅ DESeq2-ready count matrix saved to: $OUT_MATRIX"
