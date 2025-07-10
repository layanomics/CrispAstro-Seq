# 🧠 CrispAstro-Seq: Bulk RNA-seq Pipeline for CRISPR Knockouts in Neural Lineage

**CrispAstro-Seq** is a modular, reproducible pipeline for analyzing the transcriptional impact of CRISPR-mediated gene knockouts in human neural progenitors and astrocytes. The workflow processes raw RNA-seq data into biologically meaningful visualizations using cutting-edge tools and industry best practices.

> 📍 Based on GEO dataset [GSE185726](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185726)  
> ☁️ Cloud + Local: STAR alignment on GCP VM, DESeq2 downstream on local machine  
> 🧬 Genes visualized: **NES**, **GFAP**, **C3**

---

## 🎯 Project Objectives

- Reproduce and validate results from a published CRISPR-RNA-seq study
- Quantify transcriptional effects of specific gene knockouts on neural/astrocytic markers
- Build a scalable and transparent analysis pipeline suitable for GitHub portfolio

---

## 🧬 Workflow Overview

The pipeline is split into upstream (cloud) and downstream (local) stages:

```
📥 Download → 🧼 QC → 🧬 STAR Indexing → 🧲 Alignment
→ 📊 Count Merging → 📤 Transfer → 🧮 DESeq2 → 📈 Visualization
```

📄 Full breakdown: [`docs/steps.md`](docs/steps.md)  
📊 Workflow diagram: ![Pipeline](docs/crispastro_seq_pipeline.png)

---

## 🧪 Tools & Technologies

| Step                 | Tool(s) Used                        |
|----------------------|-------------------------------------|
| Data Acquisition     | `fasterq-dump` (SRA Toolkit)        |
| QC & Trimming        | `FastQC`, `Fastp`, `MultiQC`        |
| Alignment            | `STAR` (GRCh38.p12, GENCODE v38)    |
| Quantification       | STAR gene-level reverse counts      |
| DE Analysis          | `DESeq2` (R/Bioconductor)           |
| Plotting             | `ggplot2`, `ggsignif`               |
| Scripting            | Bash, R, Jupyter                    |
| Platform             | Google Cloud VM + Local machine     |

---

## 📁 Project Structure

```
CrispAstro-Seq/
├── 00_data/                    # raw, trimmed, aligned files
│   ├── raw/                    # raw FASTQs
│   ├── trimmed/                # fastp-cleaned reads
│   ├── aligned/               # STAR BAMs, Logs, ReadCounts
├── 01_qc/                      # FastQC and Fastp reports
├── 03_quantification/          # Merged gene-level count matrix
├── 04_differential_expression/
│   ├── analysis_deseq2.R
│   └── gene_barplots.R
├── notebooks/                  # Jupyter/Quarto Notebooks
├── scripts/                    # All shell scripts used in pipeline
├── results/                    # Plots, DE results, MultiQC
├── docs/                       # steps.md, pipeline diagram
├── config/                     # config.sh, path definitions
└── README.md
```

---

## 🔍 Key Results

Barplots of normalized expression for selected neurodevelopmental and astrocytic markers:

- **🧬 Nestin (NES)** — neural progenitor marker  
- **🔥 C3** — inflammation/reactive astrocyte marker  
- **🌟 GFAP** — mature astrocyte identity marker

| Gene  | Plot |
|-------|------|
| NES   | ![NES](results/plots/nes_barplot.png) |
| C3    | ![C3](results/plots/c3_barplot.png)   |
| GFAP  | ![GFAP](results/plots/gfap_barplot.png) |

---

## 🚀 How to Run This Project

### ⚙️ Requirements

- Conda or Mamba (for reproducible environments)
- Google Cloud (for VM-based alignment)
- R ≥ 4.2 + `DESeq2`, `ggplot2`, `ggsignif`

### 🧪 Example Commands

```bash
# Run FastQC on raw data
bash scripts/run_fastqc.sh raw

# Run Fastp trimming
bash scripts/run_fastp.sh

# STAR alignment (per sample)
bash scripts/run_star_alignment.sh

# Merge STAR counts
bash scripts/merge_star_counts.sh
```

In R (local):
```r
dds <- DESeqDataSetFromMatrix(countData, colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
plotCounts(dds, gene="ENSG00000132688", intgroup="condition")  # NES
```

---

## 📦 Reproducibility & Next Steps

- [x] Modular Bash scripts with dynamic thread handling  
- [x] Clean `.gitignore` and file organization  
- [x] Verified STAR alignment matches published data  
- [ ] Nextflow automation (in progress)  
- [ ] Volcano plot & PCA  
- [ ] Blog post or Quarto-based interactive report

---

## 📚 Dataset & Reference

- **GEO:** [GSE185726](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185726)  
- **Original paper:** [MECP2-deficient astrocytes alter neural development](#)

---

## 🤝 Acknowledgments

This project is part of my genomics data science portfolio, built to demonstrate hands-on expertise in:
- CRISPR knockout transcriptomics
- Reproducible cloud-scale RNA-seq analysis
- GitHub documentation & scientific storytelling

---

## 📫 Contact

Interested in collaborating or hiring?  
Find me on [GitHub](https://github.com/layanomics) or feel free to star this repo ⭐

---
