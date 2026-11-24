# Age-dependent proteomic analyses

## Overview

Age-associated protein abundance changes in UDCs and UDC-iNs were assessed using nine paired adult healthy-control donors (18 samples), excluding YoungControl samples.

## Input files

- `20230316_iN_batch_All_ys_Proteins_KM edited.xlsx`
- `Proteome_対応表(2).xlsx`

Place both workbooks on the Desktop or set the `PROTEOMICS_DATA_DIR` environment variable to their directory. Outputs are written under that directory. Use the same workbook version for all scripts.

## Scripts and execution order

1. `01_Age_DEP_GO.R`: retains proteins with finite positive abundance in all 18 samples, log2-transforms abundance and excludes zero-variance proteins. Fits a joint limma model (`~ Age10 * CellType + Sex + Batch`) with donor blocking using `duplicateCorrelation`. Age is centered and expressed per decade. Three contrasts are derived from the same fit: UDC age slope, UDC-iN age slope and their difference (UDC-iN minus UDC). Empirical Bayes uses `trend=TRUE, robust=FALSE`. Protein-level BH FDR < 0.05, without an effect-size cutoff. GO BP over-representation analysis uses mapped tested Entrez genes as background, separately for positive and negative effects, with BH FDR < 0.05. All significant GO terms are plotted; top-ten tables are also saved.
2. `02_Age_Ribosomal_Proteins.R`: reads the original results from script 01 and selects all significant cytosolic RPL/RPS genes in the UDC-iN age contrast (BH FDR < 0.05), regardless of direction. Duplicate gene entries are resolved by largest absolute original slope, then smallest FDR. A visualization-only fit produces sex- and batch-adjusted abundance plots and regression lines; original significance values are retained. This is a significance-selected descriptive panel, not independent validation. The panel size is not fixed at 15.
3. `03_Age_GSEA.R`: reads the three original limma result tables without refitting or protein-level significance filtering. Maps SYMBOL to ENTREZID using `bitr()` and averages moderated t-statistics per Entrez gene. GO BP and Reactome GSEA use fgsea multilevel, `nPermSimple=10000`, `eps=0`, gene-set sizes 10–500, and BH FDR < 0.05. Top-ten terms per direction are selected by adjusted P and plotted using NES magnitude.

Run script 01 first; scripts 02 and 03 both consume its output directory `Proteome_HC_Age_JointModel_Sex_WithBatch`. Script 03 normally uses the SYMBOL column produced by script 01. Its original workbook fallback assumes worksheet row IDs, whereas script 01 assigns IDs after accession filtering; use the original CSVs with SYMBOL present.

```sh
Rscript 01_Age_DEP_GO.R
Rscript 02_Age_Ribosomal_Proteins.R
Rscript 03_Age_GSEA.R
```

Required R packages are listed and installed if missing within each script. The scripts retain their original analysis settings and output filenames. Figures are exported as PNG and/or PDF according to the original scripts.
