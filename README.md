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

   ## Bioinformatics Workflow

The bulk RNA-seq analysis was performed using the following workflow:

1. **Data Download:** `prefetch` and `fasterq-dump` (SRA Toolkit) were used to download the paired-end raw sequencing data from NCBI SRA (`PRJNA108950`).
2. **Quality Control:** `FastQC` was used to assess the base quality scores, GC content, and adapter contamination of raw sequencing reads.
3. **Quality Trimming:** `Trimmomatic` (v0.39) was used with NexteraPE adapter sequences (`nextera_adapter.fa`) for adapter removal and quality trimming, followed by post-trimming `FastQC` verification.
4. **Reference Genome Alignment:** `HISAT2` was used to align high-quality paired-end reads to the human reference genome (`hg19`).
5. **SAM/BAM Processing:** `SAMtools` was used to convert SAM alignments to BAM format, coordinate-sort the reads, and generate BAM index files.
6. **Reference GTF Harmonization & Transcript Quantification:** Chromosome prefixes (`chr`) were standardized in `hg19.ncbiRefSeq.gtf` to match the `hg19_NoChr` index, and `StringTie` (`-e`) was used to assemble transcripts and estimate abundance levels.
7. **Count Matrix Generation:** `prepDE.py` was used to aggregate StringTie transcript estimates across all samples into a raw gene-level count matrix (`gene_count_matrix.csv`).
8. **Differential Expression Analysis:** `DESeq2` (R) was used to normalize counts and model differential gene expression (`~ condition`, with `control` as the reference group).
9. **DEG Filtering & Gene List Extraction:** Significantly differentially expressed genes were filtered (`padj < 0.05` and `|log2FoldChange| > 1`), identifying **173 upregulated** and **65 downregulated** genes in SLE compared to controls.
10. **Functional Enrichment & Visualization:** `gprofiler2` (R) was used for GO (`BP`, `MF`, `CC`), `KEGG`, and `Reactome` enrichment analyses, alongside `ggplot2` and `pheatmap` for PCA, Volcano, Heatmap, and enrichment dotplot visualizations.




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
Rscript scripts/01_deseq_script.R
Rscript scripts/02_FUNCTIONAL_ENRICHMENT.R
```
