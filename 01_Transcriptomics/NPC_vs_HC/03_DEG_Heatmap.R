# ================================================================
# UDC-iN | NPC vs HC | Top 50 DEG Heatmap
# NATURE MAIN FIGURE | UNIFIED DEG ANALYSIS
#
# DEG analysis matches Volcano / GO:
#   Raw count >= 10 in >= 3 samples
#   TMM + voom + duplicateCorrelation + limma
#   Design: ~ Diagnosis
#   BH-adjusted P < 0.05 and |log2FC| > 1
#
# Heatmap:
#   Top 50 DEGs ranked by adjusted P
#   TMM log2CPM -> row Z-score
#   NPC-up first, NPC-down second
#   Preserve HC/NPC, donor, batch order
#
# Outputs: PNG 600 dpi, PDF, CSV
# ================================================================

# 0. Packages ----------------------------------------------------

cran <- c("readxl", "pheatmap")
bioc <- c("edgeR", "limma")

for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(readxl)
library(edgeR)
library(limma)
library(pheatmap)

# 1. Files and settings ------------------------------------------

desktop <- path.expand(Sys.getenv("RNASEQ_DATA_DIR", unset = "~/Desktop"))

count_file <- file.path(
  desktop, "RNAseq_gene_count_matrix_26samples_GEO.csv"
)

meta_file <- file.path(
  desktop, "RNA-seq_Donor information_2.xlsx"
)

outdir <- file.path(
  desktop, "UDCiN_NPC_vs_HC_Top50_DEG_Heatmap_Unified"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(count_file), file.exists(meta_file))

min_count <- 10
min_samples <- 3
padj_cutoff <- 0.05
lfc_cutoff <- 1
top_n <- 50L
z_limit <- 1.5

# 2. Metadata ----------------------------------------------------

meta_raw <- as.data.frame(readxl::read_excel(meta_file))

required <- c(
  "FASTQ name", "Donor ID", "Age",
  "Condition", "Cell type", "Batch"
)

if (!all(required %in% names(meta_raw))) {
  stop("Missing metadata columns: ",
       paste(setdiff(required, names(meta_raw)), collapse = ", "))
}

meta <- data.frame(
  Sample = trimws(as.character(meta_raw[["FASTQ name"]])),
  Donor = trimws(as.character(meta_raw[["Donor ID"]])),
  Age = as.numeric(meta_raw[["Age"]]),
  Diagnosis = toupper(trimws(as.character(meta_raw[["Condition"]]))),
  CellType = trimws(as.character(meta_raw[["Cell type"]])),
  Batch = trimws(as.character(meta_raw[["Batch"]])),
  stringsAsFactors = FALSE
)

meta <- meta[
  meta$CellType == "UDC-iN" &
    meta$Diagnosis %in% c("HC", "NPC"),
  ,
  drop = FALSE
]

if (nrow(meta) == 0) stop("No UDC-iN samples found.")
if (anyDuplicated(meta$Sample)) stop("Duplicate sample names.")
if (anyNA(meta[, c("Sample", "Donor", "Age", "Diagnosis")])) {
  stop("Missing sample metadata.")
}
if (!all(c("HC", "NPC") %in% meta$Diagnosis)) {
  stop("Both HC and NPC samples are required.")
}
if (any(tapply(meta$Diagnosis, meta$Donor,
               function(x) length(unique(x))) != 1L)) {
  stop("Donor appears in multiple diagnosis groups.")
}

# HC first, NPC second; group replicates by donor
meta <- meta[
  order(
    match(meta$Diagnosis, c("HC", "NPC")),
    meta$Age,
    meta$Donor,
    meta$Batch
  ),
  ,
  drop = FALSE
]

rownames(meta) <- NULL
sample_order <- meta$Sample

write.csv(
  meta,
  file.path(outdir, "Included_Samples.csv"),
  row.names = FALSE
)

# 3. Raw count matrix --------------------------------------------

raw <- read.csv(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Gene" %in% names(raw)) {
  stop("Gene column missing.")
}

raw$Gene <- trimws(as.character(raw$Gene))
raw <- raw[!is.na(raw$Gene) & nzchar(raw$Gene), , drop = FALSE]

if (length(setdiff(sample_order, names(raw))) > 0) {
  stop("Sample names missing from count matrix.")
}

# Match Volcano/GO workflow: sum duplicated gene identifiers
sample_cols <- sample_order

raw[sample_cols] <- lapply(
  raw[sample_cols],
  function(x) as.numeric(as.character(x))
)

if (anyNA(raw[, sample_cols, drop = FALSE])) {
  stop("Non-numeric count values.")
}

counts <- rowsum(
  as.matrix(raw[, sample_cols, drop = FALSE]),
  group = raw$Gene,
  reorder = FALSE
)

counts <- counts[, sample_order, drop = FALSE]

if (
  any(!is.finite(counts)) ||
  any(counts < 0) ||
  any(abs(counts - round(counts)) > 1e-8)
) {
  stop("Counts must be nonnegative integers.")
}

storage.mode(counts) <- "integer"

# 4. SAME filtering as Volcano / GO ------------------------------

keep <- rowSums(counts >= min_count) >= min_samples

counts_filtered <- counts[keep, , drop = FALSE]

if (nrow(counts_filtered) < 100) {
  stop("Too few genes after filtering.")
}

cat("Genes retained:", nrow(counts_filtered), "\n")

# 5. TMM normalization -------------------------------------------

y <- edgeR::DGEList(counts = counts_filtered)
y <- edgeR::calcNormFactors(y, method = "TMM")

# 6. limma-voom DEG analysis -------------------------------------

diagnosis <- factor(meta$Diagnosis, levels = c("HC", "NPC"))
donor <- factor(meta$Donor)

design <- model.matrix(~ diagnosis)
colnames(design) <- c("Intercept", "NPC_vs_HC")

v <- limma::voom(y, design, plot = FALSE)

has_repeats <- any(table(donor) > 1)

if (has_repeats) {

  corfit <- limma::duplicateCorrelation(
    v, design, block = donor
  )

  if (!is.finite(corfit$consensus.correlation)) {
    stop("Donor correlation estimation failed.")
  }

  v <- limma::voom(
    y,
    design,
    block = donor,
    correlation = corfit$consensus.correlation,
    plot = FALSE
  )

  corfit <- limma::duplicateCorrelation(
    v, design, block = donor
  )

  fit <- limma::lmFit(
    v,
    design,
    block = donor,
    correlation = corfit$consensus.correlation
  )

} else {

  fit <- limma::lmFit(v, design)

}

fit <- limma::eBayes(fit)

de <- limma::topTable(
  fit,
  coef = "NPC_vs_HC",
  number = Inf,
  sort.by = "P"
)

de$Gene <- rownames(de)

de <- de[
  ,
  c("Gene", setdiff(names(de), "Gene")),
  drop = FALSE
]

write.csv(
  de,
  file.path(outdir, "All_Genes_NPC_vs_HC.csv"),
  row.names = FALSE
)

# 7. Select top 50 DEGs ------------------------------------------

sig <- de[
  is.finite(de$adj.P.Val) &
    de$adj.P.Val < padj_cutoff &
    abs(de$logFC) > lfc_cutoff,
  ,
  drop = FALSE
]

sig <- sig[
  order(
    sig$adj.P.Val,
    -abs(sig$logFC),
    sig$Gene
  ),
  ,
  drop = FALSE
]

if (nrow(sig) == 0) {
  stop("No significant DEGs.")
}

chosen <- head(sig, top_n)

if (nrow(chosen) < top_n) {
  warning("Only ", nrow(chosen), " DEGs meet the thresholds.")
}

write.csv(
  chosen,
  file.path(outdir, "Selected_Top50_DEGs.csv"),
  row.names = FALSE
)

cat("Total significant DEGs:", nrow(sig), "\n")
cat("Selected genes:", nrow(chosen), "\n")

# 8. TMM log2CPM and row Z-score ---------------------------------

logcpm <- edgeR::cpm(
  y,
  log = TRUE,
  prior.count = 2
)

expr <- logcpm[
  chosen$Gene,
  sample_order,
  drop = FALSE
]

z <- t(scale(t(expr)))

if (any(!is.finite(z))) {
  stop("Selected genes contain zero-variance expression.")
}

write.csv(
  data.frame(Gene = rownames(expr), expr, check.names = FALSE),
  file.path(outdir, "Top50_TMM_log2CPM.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Gene = rownames(z), z, check.names = FALSE),
  file.path(outdir, "Top50_RowZscore.csv"),
  row.names = FALSE
)

# 9. Order rows: NPC-up then NPC-down -----------------------------

ord <- order(
  chosen$logFC < 0,
  chosen$adj.P.Val,
  -abs(chosen$logFC)
)

chosen <- chosen[ord, , drop = FALSE]
expr <- expr[chosen$Gene, , drop = FALSE]
z <- z[chosen$Gene, , drop = FALSE]

direction <- factor(
  ifelse(
    chosen$logFC > 0,
    "Higher in NPC",
    "Lower in NPC"
  ),
  levels = c("Higher in NPC", "Lower in NPC")
)

ann_row <- data.frame(Direction = direction)
rownames(ann_row) <- chosen$Gene

ann_col <- data.frame(
  Diagnosis = factor(meta$Diagnosis, levels = c("HC", "NPC")),
  Donor = factor(meta$Donor, levels = unique(meta$Donor))
)

rownames(ann_col) <- sample_order

# 10. Annotations and colors -------------------------------------

donor_palette <- setNames(
  grDevices::hcl.colors(length(unique(meta$Donor)), "Set 2"),
  unique(meta$Donor)
)

ann_colors <- list(
  Diagnosis = c(
    HC = "#A9B9C9",
    NPC = "#E58C83"
  ),
  Donor = donor_palette,
  Direction = c(
    "Higher in NPC" = "#C46A63",
    "Lower in NPC" = "#5D86B3"
  )
)

if (
  anyNA(meta$Batch) ||
  any(!nzchar(meta$Batch)) ||
  any(tolower(meta$Batch) %in% c("na", "nan", "null"))
) {
  stop("Batch metadata is missing.")
}

batch_label <- ifelse(
  grepl("^batch\\s*", meta$Batch, ignore.case = TRUE),
  meta$Batch,
  paste0("Batch ", meta$Batch)
)

labels <- paste0(
  meta$Diagnosis,
  ". ",
  format(meta$Age, trim = TRUE, scientific = FALSE),
  "y, ",
  batch_label
)

donor_sizes <- as.integer(
  table(factor(meta$Donor, levels = unique(meta$Donor)))
)

donor_gaps <- head(cumsum(donor_sizes), -1)
group_gap <- sum(meta$Diagnosis == "HC")

gaps_col <- sort(unique(c(donor_gaps, group_gap)))

up_n <- sum(chosen$logFC > 0)

gaps_row <- if (
  up_n > 0 && up_n < nrow(z)
) up_n else NULL

# 11. Heatmap color scales ---------------------------------------

heat_colors <- grDevices::colorRampPalette(
  c(
    "#174A8B",
    "#4385BE",
    "#D7E5F2",
    "#FFFFFF",
    "#F3D7D3",
    "#DA7368",
    "#A5232C"
  )
)(201)

heat_breaks <- seq(
  -z_limit,
  z_limit,
  length.out = 202
)

# 12. Plot function ----------------------------------------------

plot_heatmap <- function(mat, cols, brks, main, ticks) {

  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    labels_col = labels,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    fontsize = 11,
    fontsize_row = 9.5,
    fontsize_col = 8,
    angle_col = 45,
    cellwidth = 55,
    cellheight = 11.5,
    border_color = NA,
    color = cols,
    breaks = brks,
    legend_breaks = ticks,
    main = main,
    silent = TRUE
  )

}

save_plot <- function(hm, stem, width = 19.5, height = 13.5) {

  pdf(
    file.path(outdir, paste0(stem, ".pdf")),
    width = width,
    height = height,
    useDingbats = FALSE
  )

  grid::grid.newpage()
  grid::grid.draw(hm$gtable)
  dev.off()

  png(
    file.path(outdir, paste0(stem, "_600dpi.png")),
    width = width,
    height = height,
    units = "in",
    res = 600,
    bg = "white"
  )

  grid::grid.newpage()
  grid::grid.draw(hm$gtable)
  dev.off()

}

# 13. MAIN FIGURE: Row Z-score -----------------------------------

hm <- plot_heatmap(
  z,
  heat_colors,
  heat_breaks,
  "UDC-iN | NPC vs HC | Top 50 DEGs",
  c(-1.5, -1, 0, 1, 1.5)
)

save_plot(
  hm,
  "MainFigure_Top50_DEG_RowZscore"
)

# 14. SUPPLEMENT: Actual log2CPM -------------------------------

range_expr <- range(expr)

abs_colors <- grDevices::colorRampPalette(
  c(
    "#F7F7F7",
    "#D4E3EF",
    "#73A2C6",
    "#1D497C"
  )
)(201)

abs_breaks <- seq(
  range_expr[1],
  range_expr[2],
  length.out = 202
)

ticks <- pretty(range_expr, n = 5)

ticks <- ticks[
  ticks >= range_expr[1] &
    ticks <= range_expr[2]
]

hm_abs <- plot_heatmap(
  expr,
  abs_colors,
  abs_breaks,
  "UDC-iN | NPC vs HC | Top 50 DEGs | log2 CPM",
  ticks
)

save_plot(
  hm_abs,
  "Supplement_Top50_DEG_log2CPM"
)

# 15. Audit and summary ------------------------------------------

writeLines(
  c(
    "UDC-iN NPC vs HC Top 50 DEG Heatmap",
    "Filter: raw count >=10 in >=3 samples",
    "Normalization: TMM",
    "DEG: limma-voom",
    "Design: ~ Diagnosis",
    paste0("Donor correlation: ", has_repeats),
    "Test: empirical Bayes moderated t",
    "FDR: BH",
    "DEG: adjusted P <0.05 and |log2FC|>1",
    "Top 50: ranked by adjusted P",
    "Heatmap: TMM log2CPM, row Z-score"
  ),
  file.path(outdir, "Analysis_Settings.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_SessionInfo.txt")
)

cat("\n========================================\n")
cat("HEATMAP COMPLETED\n")
cat("========================================\n")
cat("HC samples:", sum(meta$Diagnosis == "HC"), "\n")
cat("NPC samples:", sum(meta$Diagnosis == "NPC"), "\n")
cat("HC donors:", length(unique(meta$Donor[meta$Diagnosis == "HC"])), "\n")
cat("NPC donors:", length(unique(meta$Donor[meta$Diagnosis == "NPC"])), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
cat("Significant DEGs:", nrow(sig), "\n")
cat("Heatmap genes:", nrow(chosen), "\n")
cat("Output:", outdir, "\n")
cat("========================================\n")

print(hm)
