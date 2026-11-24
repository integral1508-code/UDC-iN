# UDC-iN vs UDC RNA-seq analysis

R scripts for comparing directly induced neurons derived from urine-derived cells (UDC-iNs) with parental urine-derived cells (UDCs). The input dataset comprises 22 samples from 11 paired donors.

## Requirements

R with the following packages (missing packages are installed by the scripts):

- CRAN: readxl, dplyr, tibble, ggplot2, ggrepel, stringr, pheatmap, BiocManager
- Bioconductor: DESeq2, clusterProfiler, org.Hs.eg.db, AnnotationDbi, edgeR, limma

Each script records or can record its R session information; package versions from the original analysis were not supplied. For exact reproduction, use the original session versions when available.

## Input data

Create a `data/` directory in this folder and place the following files in it:

| File | Description |
| --- | --- |
| `RNAseq_gene_count_matrix_22samples_GEO.csv` | Raw nonnegative integer counts; first identifier column named `Gene`, remaining columns named by sample ID. The heatmap requires unique gene symbols. |
| `RNA-seq_Donor information.xlsx` | Sample metadata for all three analyses. |

Metadata columns: `FASTQ name`, `Donor ID`, `Age`, `Condition`, `Cell type`. Cell types must be `UDC` and `UDC-iN`. PCA and DEG/GO require exact sample matching and one sample of each cell type per donor. The heatmap selects HC pairs from its metadata. All three scripts use the same metadata file. Ensure appropriate authorization for any donor-level information. Data are not included in this repository package.

## Usage

Run from this directory:

```sh
Rscript 01_PCA.R
Rscript 02_DEG_GO.R
Rscript 03_Marker_Heatmap.R
```

Scripts are independent. Outputs are written to `results/` in separate subdirectories. Override the input and output directories using the environment variables `UDC_INPUT_DIR` and `UDC_OUTPUT_DIR` if needed.

## Analysis overview

| Script | Methods | Outputs |
| --- | --- | --- |
| `01_PCA.R` | Raw count >=10 in >=3 samples; DESeq2 size factors and blind VST; centered, unscaled PCA of all retained nonconstant genes; 95% covariance-based data ellipses. | PCA coordinates, variance, gene list, VST matrix and PNG figures. |
| `02_DEG_GO.R` | Same count filter; paired DESeq2 Wald test (`~ Donor + CellType`); UDC-iN versus UDC; BH-adjusted P <0.05 and absolute log2 fold change >1. GO BP over-representation analysis separately for UP and DOWN genes; mapped genes with nonmissing adjusted P as background; BH-adjusted P <0.05; gene-set size 10–500. | DEG tables, volcano plot, GO tables, top-10 bar/dot plots, normalized counts and analysis settings. |
| `03_Marker_Heatmap.R` | TMM normalization using all genes, log2 CPM (prior count 2), row-wise Z scores of 10 curated markers; display saturation at +/-1.5 only. Paired limma-voom analysis is supplementary QC. | Marker expression and unclipped Z scores, QC tables, paired statistics, PDF and 600-dpi PNG heatmaps. |

Positive fold changes indicate higher expression in UDC-iNs. PCA group-label positions are fixed to match the original figure. The marker heatmap is descriptive and is not a significance test.

## Reproducibility notes

The original statistical settings and plotting parameters are retained. Local desktop paths and automatic macOS opening commands were replaced with portable paths and batch-friendly execution. Empty GO results are handled without failing the significance filter. Scripts have not been executed against the original input data in preparing this package.
