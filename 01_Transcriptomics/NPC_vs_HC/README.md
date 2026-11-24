# NPC vs HC transcriptomic analyses

This directory contains R scripts for comparing Niemann–Pick disease type C (NPC) and healthy control (HC) samples in UDC-iNs and induced pluripotent stem cell-derived induced neurons (iPSC-iNs).

## Scripts

- `01_PCA.R`: joint PCA of UDC-iNs and iPSC-iNs using all batches, DESeq2 variance-stabilizing transformation (blind = TRUE), and all filtered nonconstant genes. No batch correction is applied. Ellipses represent 95% data ellipses by cell type and diagnosis.
- `02_DEG_GO.R`: differential gene expression, volcano plots, and GO Biological Process over-representation analysis, performed separately for each cell type. TMM normalization, limma-voom, donor blocking using duplicateCorrelation where repeated samples exist, and empirical Bayes moderation are used. The model includes diagnosis without age adjustment. DEGs are defined by BH-adjusted P < 0.05 and absolute log2 fold change > 1. GO analysis is performed separately for upregulated and downregulated genes using mapped tested genes as background and BH-adjusted P < 0.05.
- `03_DEG_Heatmap.R`: UDC-iN heatmaps of up to 50 significant DEGs ranked by adjusted P, using TMM log2CPM and row Z-scores. Differential expression is calculated independently using the same modeling approach.

All scripts retain genes with raw counts >= 10 in at least three samples (within each analyzed cell type for differential expression).

## Input data

Place the following files on the Desktop, or set the `RNASEQ_DATA_DIR` environment variable to their directory:

- `RNAseq_gene_count_matrix_26samples_GEO.csv`: raw gene counts, with a `Gene` column and sample names as remaining column headers.
- `RNA-seq_Donor information_2.xlsx`: sample metadata with columns `FASTQ name`, `Donor ID`, `Age`, `Condition`, `Batch`, and `Cell type`. Diagnosis values are HC and NPC; cell types are UDC-iN and iPSC-iN.

Sample identifiers must match between the input files. Input data are not included in this directory.

## Running the analyses

Run each script independently from this directory:

```sh
Rscript 01_PCA.R
Rscript 02_DEG_GO.R
Rscript 03_DEG_Heatmap.R
```

Required CRAN and Bioconductor packages are installed if missing. Each script writes results, figures and R session information to a separate output directory within the input directory.

## Interpretation

The differential expression contrast is NPC relative to HC. Positive log2 fold changes indicate higher expression in NPC; negative values indicate lower expression in NPC.
