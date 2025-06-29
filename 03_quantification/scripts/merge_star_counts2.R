# 🧬 STAR Count Merger — CrispAstro-Seq Project
# Description: Merge STAR ReadsPerGene.out.tab files into a DESeq2-ready counts matrix

library(dplyr)
library(readr)
library(tibble)

# 📁 Paths
input_dir <- "/home/layanomics/CrispAstro-Seq/00_data/paper_counts" # <--- CHANGE THIS PATH
output_path <- "/home/layanomics/CrispAstro-Seq/03_quantification/results/counts2.tsv"

# 📄 Find STAR files
files <- list.files(input_dir, pattern = "_ReadsPerGene.out.tab$", full.names = TRUE)
sample_ids <- gsub("_ReadsPerGene.out.tab", "", basename(files))

cat("📦 Files detected:\n")
print(files)

cat("🧬 Sample IDs extracted:\n")
print(sample_ids)

cat("📊 Starting merge process...\n")

# 🧪 Read the first file and keep only gene counts (column 4), skipping first 4 STAR summary rows
df <- read_tsv(files[1], col_names = FALSE) %>% # REMOVE 'show_col_types = FALSE'
    slice(-1:-4) %>%
    select(GeneID = X1, !!sample_ids[1] := X2) # Keep 'X4' as confirmed by your tool

# 🔁 Loop through remaining files
for (i in 2:length(files)) {
    next_df <- read_tsv(files[i], col_names = FALSE) %>% # REMOVE 'show_col_types = FALSE'
        slice(-1:-4) %>%
        select(GeneID = X1, !!sample_ids[i] := X2) # Keep 'X4'

    df <- full_join(df, next_df, by = "GeneID")
}

cat("🔍 Preview before collapsing duplicates:\n")
print(head(df, 3))

# 🧼 Strip `.version` from gene IDs (e.g., ENSG00000123456.5 → ENSG00000123456)
df <- df %>% mutate(GeneID = gsub("\\..*", "", GeneID))

# 🧮 Collapse duplicated gene IDs (if multiple .version collapsed into one)
counts <- df %>%
    group_by(GeneID) %>%
    summarise(across(everything(), ~ sum(as.numeric(.), na.rm = TRUE))) %>%
    column_to_rownames("GeneID")

cat("✅ Final count matrix dimensions (genes x samples):\n")
print(dim(counts))

# 💾 Save as tab-separated file with GeneID column
write_tsv(counts %>% rownames_to_column("GeneID"), output_path)

cat("🎉 Done! Merged STAR gene count matrix saved to:\n", output_path, "\n")
