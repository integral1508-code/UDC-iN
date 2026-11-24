# Disease-associated proteomic alterations in UDC-iNs

This directory contains R scripts for comparing healthy controls (HC), mild cognitive impairment due to Alzheimer's disease (MCI), and Alzheimer's disease (AD) in UDC-iNs.

## Input data

- `20230316_iN_batch_All_ys_Proteins_KM edited.xlsx`: TMT scaled abundance workbook with protein annotations (`Accession`, `Gene Symbol`) and a `metadata` worksheet (`File Name`, `Type`, `ID2`, `Sex`, `Age`, `batch`).
- `Proteome_対応表(2).xlsx`: diagnosis correspondence workbook without a header, containing sample names, diagnosis labels and age (third column, used by script 01).

Place both workbooks in `~/Desktop` or set `PROTEOME_DATA_DIR` to their directory. Use the original workbook contents with the filenames above. Scripts install missing CRAN and Bioconductor packages as needed. R 4.1 or later is required for native pipe syntax; use compatible Bioconductor packages.

## Donor exclusion

The 72-year-old donor originally classified as AD was excluded from the analyses because the updated clinical diagnosis indicated that the donor did not have AD. Both matched UDC and UDC-iN samples from this donor were excluded, leaving 11 AD UDC-iN donors.

## Scripts and execution

Run from this directory:

```sh
Rscript 01_DEP_GO.R
Rscript 02_GSEA_GO.R
Rscript 03_RNA_Splicing_GO_Figure.R
Rscript 04_RNA_Splicing_Top20.R
```

Scripts 01 and 02 independently fit their models. Scripts 03 and 04 require the outputs of script 02. Outputs are written under the input data directory.

### 01 — Differential protein abundance and GO over-representation analysis

UDC-iNs only; one sample per donor. The matched donor pair described above is excluded, leaving 11 AD UDC-iN donors. Diagnosis is classified by label suffix; age is obtained from the correspondence workbook. Positive log2 fold changes indicate higher abundance in the first-named group (MCI vs HC, AD vs HC, AD vs MCI). Positive finite abundances across all selected samples are log2-transformed; zero-variance proteins are removed without imputation. The limma model is `~ 0 + Dx + Age10 + Sex + Batch`, with centered age per decade and empirical Bayes moderation (`trend=TRUE`, `robust=FALSE`). BH adjustment is performed across proteins separately for each contrast. DEP threshold: adjusted P < 0.05, without a fold-change cutoff. GO BP ORA is performed separately for UP and DOWN proteins, using mapped tested proteins as background, BH FDR < 0.05, gene-set sizes 10–500 and at least five mapped input genes. The top 10 terms are selected by adjusted P and displayed by descending GeneRatio. SYMBOL mapping uses UniProt fallback where needed.

### 02 — Preranked GSEA and shared GO terms

The original diagnosis label rules and age from the proteome metadata are retained. The same specified donor pair is excluded. The model includes diagnosis and centered age per decade, plus sex and TMT batch when they have multiple levels; rank-deficient designs stop rather than silently dropping covariates. limma moderated t-statistics are averaged across mapped protein rows per Entrez gene, without DEP prefiltering. GO BP GSEA uses fgsea, nPermSimple=10000, eps=0, gene-set sizes 10–500 and seed 20261016. BH correction is applied globally to finite GO P values returned across all three contrasts; within-contrast adjusted P values are retained separately. Significant terms have global FDR < 0.05. Shared terms require the same GO ID and NES direction. Top 10 terms per direction are selected by adjusted P and displayed by absolute NES.

### 03 — RNA splicing GO figure

Reads the three significant GSEA tables from script 02. Displays seven predefined RNA processing/splicing GO terms with negative NES and global FDR < 0.05 in all three comparisons. Dot size represents absolute NES and color represents negative log10 adjusted P. No GSEA or statistical model is rerun.

### 04 — Leading-edge RNA splicing proteins

Selects the 20 most negative AD vs HC gene ranking scores from the union of leading-edge genes of the seven specified GO terms. MCI is not used in selection. Measured log2 protein abundances are averaged per mapped gene and donor. A separate gene-level linear model is fitted for visualization, retaining diagnosis and adjusting for age, sex and batch (when applicable). Adjusted donor values and estimated marginal means use the pooled observed covariate distribution. Figures show adjusted donor dots with mean lines (no confidence intervals) and centered adjusted group means. No additional significance tests are performed; these plots are not independent validation.

## Reproducibility notes

Scripts retain the original analysis-specific diagnosis and age sources rather than harmonizing them. Confirm that the selected samples and ages agree across scripts using the exported metadata tables. Protein-level BH in script 01 and global GO-level BH in script 02 are distinct procedures. Figures are exported as 600 dpi PNG and PDF. Session information is exported by the main analyses.
