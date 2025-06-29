# 🧬 DESeq2 Core Analysis — CrispAstro-Seq
# Description: Runs DESeq2 and saves all intermediate objects for downstream visualization
# Based on: paper methods — filter total count < 10, padj < 0.05, WT as reference

library(DESeq2)
library(dplyr) # For %>% and other dplyr functions
library(readr) # For read_tsv, read_csv
library(tibble) # For column_to_rownames
# 📁 Paths
count_file <- "/home/layanomics/CrispAstro-Seq/03_quantification/results/counts2.tsv" # Change to the new merged file
metadata_file <- "/home/layanomics/CrispAstro-Seq/00_data/metadata/metadata.csv" # Confirmed to use your existing metadata
output_dir <- "/home/layanomics/CrispAstro-Seq/04_differential_expression/results_DESeq2_1.24.0"


# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 📥 Load input files
counts <- read_tsv(count_file) %>% # Make sure 'show_col_types = FALSE' is NOT here
  column_to_rownames("GeneID")

metadata <- read_csv(metadata_file) %>%
  column_to_rownames("sample")

# 🔤 Recode condition names for clarity
metadata$condition <- recode(metadata$condition,
  "wild-type" = "WT",
  "MECP2-R270X" = "R270X"
)

# ✅ Set reference level BEFORE DESeqDataSet
metadata$condition <- factor(metadata$condition, levels = c("WT", "R270X"))

cat("✅ Loaded counts and metadata.\n")
print(dim(counts))
print(table(metadata$condition))

# 🔍 Filter genes with total count < 10
counts <- counts[rowSums(counts) >= 10, ]

# 🧬 Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition)

# 🚀 Run DESeq2
dds <- DESeq(dds)

# 📊 Extract results
res <- results(dds)
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

# 💾 Save results and objects
saveRDS(dds, file.path(output_dir, "dds.rds"))
saveRDS(res_df, file.path(output_dir, "deseq2_results_df.rds"))
write.csv(res_df, file.path(output_dir, "deseq2_results.csv"))

# 📊 Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
saveRDS(norm_counts, file.path(output_dir, "normalized_counts.rds"))
write.csv(norm_counts, file.path(output_dir, "normalized_counts.csv"))

cat("🎉 DESeq2 analysis complete. All outputs saved to:\n", output_dir, "\n")
