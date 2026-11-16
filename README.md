# Bulk RNA-seq Differential Expression & Functional Enrichment in Systemic Lupus Erythematosus (SLE)

## Project Overview

This project analyzes bulk RNA sequencing (RNA-seq) data to investigate differential gene expression and dysregulated biological pathways between patients with **Systemic Lupus Erythematosus (SLE)** and **healthy controls**.

Using publicly available sequencing data, transcript abundances were quantified with **StringTie** and converted into a gene-level count matrix via **`prepDE.py`**, followed by differential expression modeling in **DESeq2** and functional enrichment analysis with **g:Profiler (`gprofiler2`)**.

## Study Information

- **Disease:** Systemic Lupus Erythematosus (SLE)
- **Organism:** *Homo sapiens*
- **Data Type:** Bulk RNA-seq
- **BioProject:** PRJNA108950
- **Analysis:** Differential gene expression & functional over-representation analysis

## Experimental Design

| Group | Sample IDs | Number of Samples |
| :--- | :--- | :---: |
| **Control** | `SRR28387650`, `SRR28387651`, `SRR28387652` | 3 |
| **SLE** | `SRR28387657`, `SRR28387658`, `SRR28387660` | 3 |

## Tools and Technologies

- **Upstream Quantification:** StringTie, `prepDE.py`
- **Statistical & Enrichment Analysis:** R, `DESeq2`, `gprofiler2`
- **Visualization:** `ggplot2`, `pheatmap`

## Analysis Workflow

1. **Upstream RNA-seq Processing:** Transcript assembly and quantification using StringTie.
2. **Count Matrix Generation:** Extraction of raw gene-level counts via `prepDE.py`.
3. **Differential Gene Expression:** Normalization and statistical modeling (`~ condition`, with `control` as reference) using `DESeq2`.
4. **DEG Filtering:** Identification of significant genes (`padj < 0.05` and `|log2FoldChange| > 1`).
5. **Visualization & Functional Enrichment:** PCA, Volcano plot, Top 50 DEG Heatmap, and functional enrichment (`GO:BP`, `GO:MF`, `GO:CC`, `KEGG`, `Reactome`) via `gprofiler2`.

## Key Results

Differential expression analysis (`padj < 0.05`, `|log2FC| > 1`) identified **238 significant DEGs**:

- **173 upregulated genes** in SLE vs. Control
- **65 downregulated genes** in SLE vs. Control

### Biological Interpretation (g:Profiler Enrichment)
- **Upregulated Signatures (Immune & Type I/II Interferon Activation):** Strong enrichment of **Interferon alpha/beta signaling** (`R-HSA-909733`, FDR = 1.92 × 10⁻¹⁰), **Cytokine Signaling in Immune system** (`R-HSA-1280215`, FDR = 1.92 × 10⁻¹⁰), and **defense response to virus** (`GO:0051607`, FDR = 6.55 × 10⁻⁹), driven by classic SLE Interferon-Stimulated Genes (`MX1`, `OASL`, `IFI6`, `IFI27`, `IFI44L`, `IFIT5`, `USP18`) alongside `IRF7`, `SOCS1`, `SOCS3`, `IL10`, and `FOXP3`.
- **Downregulated Signatures (Chemokine Receptor Suppression):** Coordinated downregulation of **C-C chemokine receptor activity** (`GO:0016493`, FDR = 9.46 × 10⁻⁷) and **Chemokine signaling pathway** (`KEGG:04062`, FDR = 2.35 × 10⁻⁴), driven by `CCR2`, `CCR4`, `CCR5`, `CX3CR1`, `CXCR6`, and `CCL4`.

## Repository Structure

```text
├── data/               # Sample metadata and raw gene count matrix
├── scripts/            # R scripts for DESeq2 (01) and g:Profiler enrichment (02)
└── results/
    ├── tables/         # DESeq2 result tables, DEG gene lists, and g:Profiler CSVs
    └── figures/        # PCA, Volcano, Top 50 Heatmap, and enrichment dotplots
```

## How to Run

```bash
Rscript scripts/01_deseq2_analysis.R
Rscript scripts/02_go_enrichment.R
```
