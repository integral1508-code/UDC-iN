# UDC-iN vs UDC proteomic analyses

This directory contains R scripts for paired proteomic comparisons of urine-derived cells (UDCs) and directly induced neurons derived from UDCs (UDC-iNs).

## Scripts

- `01_PCA.R`: PCA of log2-transformed TMT scaled abundances using centered, unscaled proteins and no additional batch correction. Matched donor samples are connected, and empirical 95% data ellipses are shown.
- `02_DEP_GO.R`: paired differential protein abundance analysis using limma (`~ Donor + CellType`), empirical Bayes moderation and Benjamini–Hochberg correction. Significant differentially abundant proteins are defined by adjusted P < 0.05 and absolute log2 fold change > 1. Volcano plots and GO Biological Process over-representation analyses are generated.

Both scripts select samples classified as `1.HC` in the diagnosis mapping, excluding `6.YoungControl`. The differential abundance script requires nine matched donor pairs (18 samples). Proteins must have finite, strictly positive abundances in every selected sample and nonzero variance after log2 transformation. Duplicate accession rows are retained with unique suffixes; no imputation or additional batch correction is performed.

GO analysis is performed separately for increased and decreased proteins using unique Entrez IDs mapped from tested proteins as the background. UniProt accessions are mapped using `org.Hs.eg.db` (first mapping). GO terms are filtered at BH-adjusted P < 0.05; gene set sizes are 10–500. The top ten terms are selected by adjusted P and displayed in descending GeneRatio order.

## Input data

Place both files on the Desktop, or set the `PROTEOMICS_DATA_DIR` environment variable to their directory:

- `20230316_iN_batch_All_ys_Proteins_KM edited.xlsx`: TMT scaled protein abundance workbook. The first non-metadata sheet must contain `Accession` and sample abundance columns. The `metadata` sheet must contain `File Name`, `Type`, `ID2`, `Patient`, `Dx`, `Age`, and `batch`.
- `Proteome_対応表(2).xlsx`: diagnosis mapping workbook without a header; the first three columns are sample file name, diagnosis, and age.

Sample identifiers must match across the workbooks. Cell type values in the metadata are `UDC` and `iN`. Input data are not included.

## Running the analyses

Run each script independently:

```sh
Rscript 01_PCA.R
Rscript 02_DEP_GO.R
```

Required CRAN packages: readxl, ggplot2, ggrepel, dplyr, tibble, stringr. Required Bioconductor packages: limma, clusterProfiler, org.Hs.eg.db, AnnotationDbi. Missing packages are installed by the scripts.

Each script saves results, figures (including 600-dpi PNG and PDF), and R session information in a separate output directory within the input directory.

## Interpretation

Positive log2 fold changes indicate higher protein abundance in UDC-iNs relative to UDCs; negative values indicate lower abundance.
