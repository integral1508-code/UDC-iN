# ================================================================
# RNA-seq PCA | UDC-iN vs iPSC-iN
# NATURE MAIN FIGURE | ALL BATCHES | HC/NPC ELLIPSES
#
# Input:
# ~/Desktop/RNAseq_gene_count_matrix_26samples_GEO.csv
# ~/Desktop/RNA-seq_Donor information_2.xlsx
#
# - ALL batches for BOTH cell types
# - No batch correction
# - Count >= 10 in >= 3 samples
# - VST, blind = TRUE
# - PCA using ALL filtered genes
# - HC/NPC separate 95% data ellipses
# - Labels arranged in non-overlapping columns
# - No sample-to-sample lines
# - 600 dpi PNG
# ================================================================


# ================================================================
# 0. Packages
# ================================================================

cran_packages <- c(
  "readxl", "ggplot2", "dplyr", "tibble"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    "DESeq2",
    ask = FALSE,
    update = FALSE
  )
}

library(DESeq2)
library(readxl)
library(ggplot2)
library(dplyr)
library(tibble)

set.seed(123)


# ================================================================
# 1. Paths
# ================================================================

desktop <- path.expand(Sys.getenv("RNASEQ_DATA_DIR", unset = "~/Desktop"))

count_file <- file.path(
  desktop,
  "RNAseq_gene_count_matrix_26samples_GEO.csv"
)

meta_file <- file.path(
  desktop,
  "RNA-seq_Donor information_2.xlsx"
)

outdir <- file.path(
  desktop,
  "RNAseq_PCA_AllBatches_HC_NPC_Final"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(count_file),
  file.exists(meta_file)
)


# ================================================================
# 2. Read count matrix
# ================================================================

raw <- read.csv(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("Gene" %in% names(raw))

raw <- raw[
  !is.na(raw$Gene) &
    trimws(as.character(raw$Gene)) != "",
  ,
  drop = FALSE
]

raw$Gene <- trimws(as.character(raw$Gene))

sample_cols <- setdiff(names(raw), "Gene")

raw[sample_cols] <- lapply(
  raw[sample_cols],
  function(x) as.numeric(as.character(x))
)

count_check <- as.matrix(
  raw[, sample_cols, drop = FALSE]
)

if (
  anyNA(count_check) ||
  any(!is.finite(count_check)) ||
  any(count_check < 0) ||
  any(abs(count_check - round(count_check)) > 1e-8)
) {
  stop("Count matrix contains invalid values.")
}

# Merge duplicated gene IDs

counts <- raw %>%
  group_by(Gene) %>%
  summarise(
    across(all_of(sample_cols), sum),
    .groups = "drop"
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

if (any(counts > .Machine$integer.max)) {
  stop("Count exceeds integer limit.")
}

storage.mode(counts) <- "integer"


# ================================================================
# 3. Metadata
# ================================================================

meta_raw <- read_excel(meta_file)

required_cols <- c(
  "FASTQ name",
  "Donor ID",
  "Age",
  "Condition",
  "Batch",
  "Cell type"
)

if (!all(required_cols %in% names(meta_raw))) {
  stop(
    "Missing metadata columns: ",
    paste(
      setdiff(required_cols, names(meta_raw)),
      collapse = ", "
    )
  )
}

meta <- meta_raw %>%
  transmute(
    Sample = trimws(as.character(`FASTQ name`)),
    Donor = trimws(as.character(`Donor ID`)),
    Age = as.numeric(Age),
    Diagnosis = toupper(
      trimws(as.character(Condition))
    ),
    Batch = trimws(as.character(Batch)),
    CellType = trimws(as.character(`Cell type`))
  )

if (anyDuplicated(meta$Sample)) {
  stop("Duplicated sample names in metadata.")
}

if (!setequal(meta$Sample, colnames(counts))) {

  cat("\nOnly in count matrix:\n")
  print(setdiff(colnames(counts), meta$Sample))

  cat("\nOnly in metadata:\n")
  print(setdiff(meta$Sample, colnames(counts)))

  stop("Count matrix and metadata do not match.")
}

meta <- meta[
  match(colnames(counts), meta$Sample),
  ,
  drop = FALSE
]

stopifnot(
  identical(meta$Sample, colnames(counts))
)

if (
  anyNA(meta$Sample) ||
  anyNA(meta$Donor) ||
  anyNA(meta$Age) ||
  anyNA(meta$Diagnosis) ||
  anyNA(meta$Batch) ||
  anyNA(meta$CellType)
) {
  stop("Missing or invalid metadata.")
}

if (!all(meta$Diagnosis %in% c("HC", "NPC"))) {
  stop("Diagnosis must be HC or NPC.")
}

if (!all(meta$CellType %in% c("UDC-iN", "iPSC-iN"))) {
  stop("Unexpected Cell type.")
}

# ================================================================
# 4. Include ALL batches
# ================================================================

# IMPORTANT:
# No batch filtering is performed.
# Every sample in the count matrix is retained.

meta$BatchNumber <- sub(
  "^\\s*[Bb]atch\\s*",
  "",
  meta$Batch
)

meta$CellType <- factor(
  meta$CellType,
  levels = c("UDC-iN", "iPSC-iN")
)

meta$Diagnosis <- factor(
  meta$Diagnosis,
  levels = c("HC", "NPC")
)

meta$Donor <- factor(meta$Donor)

meta$Label <- paste0(
  meta$Diagnosis,
  ". ",
  format(meta$Age, trim = TRUE),
  "y, Batch ",
  meta$BatchNumber
)

write.csv(
  meta,
  file.path(outdir, "Included_Samples.csv"),
  row.names = FALSE
)

cat("\n====================================\n")
cat("ALL BATCHES INCLUDED\n")
cat("====================================\n")

cat("Samples:", nrow(meta), "\n")

print(
  table(
    meta$CellType,
    meta$Diagnosis,
    meta$BatchNumber
  )
)


# ================================================================
# 5. Low-expression filtering
# ================================================================

keep <- rowSums(
  counts >= 10
) >= 3

counts_filtered <- counts[
  keep,
  ,
  drop = FALSE
]

if (nrow(counts_filtered) < 2) {
  stop("Too few genes after filtering.")
}

cat(
  "Filtered genes:",
  nrow(counts_filtered),
  "\n"
)


# ================================================================
# 6. VST
# ================================================================

coldata <- as.data.frame(meta)

rownames(coldata) <- coldata$Sample

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)

vsd <- varianceStabilizingTransformation(
  dds,
  blind = TRUE
)

vst_mat <- assay(vsd)


# ================================================================
# 7. PCA
# ================================================================

pca_input <- t(vst_mat)

gene_var <- apply(
  pca_input,
  2,
  var
)

valid_genes <- is.finite(gene_var) &
  gene_var > 0

pca_input <- pca_input[
  ,
  valid_genes,
  drop = FALSE
]

if (ncol(pca_input) < 2) {
  stop("Too few genes for PCA.")
}

pca <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)

variance <- 100 *
  pca$sdev^2 /
  sum(pca$sdev^2)

pca_df <- as.data.frame(
  pca$x[, 1:2, drop = FALSE]
) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by = "Sample")


# ================================================================
# 8. Ellipse function
# ================================================================

make_ellipse <- function(
    df,
    level = 0.95,
    n = 300
) {

  xy <- as.matrix(
    df[, c("PC1", "PC2")]
  )

  mu <- colMeans(xy)
  sigma <- cov(xy)

  eig <- eigen(
    sigma,
    symmetric = TRUE
  )

  if (any(eig$values <= 0)) {
    return(NULL)
  }

  theta <- seq(
    0,
    2 * pi,
    length.out = n
  )

  circle <- rbind(
    cos(theta),
    sin(theta)
  )

  transform_mat <- eig$vectors %*%
    diag(
      sqrt(eig$values),
      nrow = 2
    )

  radius <- sqrt(
    qchisq(level, df = 2)
  )

  coords <- sweep(
    radius * transform_mat %*% circle,
    1,
    mu,
    "+"
  )

  data.frame(
    PC1 = coords[1, ],
    PC2 = coords[2, ]
  )
}


# ================================================================
# 9. HC / NPC ellipses
# ================================================================

ellipse_groups <- pca_df %>%
  group_by(CellType, Diagnosis) %>%
  group_split()

ellipse_list <- lapply(
  ellipse_groups,
  function(df) {

    if (nrow(df) < 3) {
      return(NULL)
    }

    ell <- make_ellipse(df)

    if (is.null(ell)) {
      return(NULL)
    }

    ell$CellType <- as.character(df$CellType[1])
    ell$Diagnosis <- as.character(df$Diagnosis[1])

    ell
  }
)

ellipse_df <- bind_rows(ellipse_list)

if (nrow(ellipse_df) == 0) {
  stop("No valid ellipse groups.")
}

ellipse_df$CellType <- factor(
  ellipse_df$CellType,
  levels = c("UDC-iN", "iPSC-iN")
)

ellipse_df$Diagnosis <- factor(
  ellipse_df$Diagnosis,
  levels = c("HC", "NPC")
)

ellipse_df$EllipseGroup <- interaction(
  ellipse_df$CellType,
  ellipse_df$Diagnosis,
  drop = TRUE
)


# ================================================================
# 10. Colors
# ================================================================

cell_colors <- c(
  "UDC-iN" = "#CC6859",
  "iPSC-iN" = "#3C78A6"
)

diagnosis_linetypes <- c(
  "HC" = "solid",
  "NPC" = "dashed"
)


# ================================================================
# 11. Deterministic label layout
# ================================================================

# Instead of relying on ggrepel, sample labels are
# assigned fixed columns and evenly spaced vertically.
#
# This prevents label-label overlap within each column.
#
# Layout:
#
# UDC-iN:
#   LEFT and RIGHT columns
#
# iPSC-iN:
#   LEFT and RIGHT columns
#
# All label positions are based on PCA coordinates.
# PCA itself is unchanged.

pc1_span <- diff(range(pca_df$PC1))
pc2_span <- diff(range(pca_df$PC2))

# Label vertical spacing in PCA coordinate units.
# Increase if labels appear too close.

label_step <- max(10, 0.045 * pc2_span)

# Label columns:
# Keep UDC and iPSC labels in separate horizontal regions.

udc_data <- pca_df %>%
  filter(CellType == "UDC-iN") %>%
  arrange(PC2, PC1) %>%
  mutate(
    Side = ifelse(
      row_number() %% 2 == 0,
      "left",
      "right"
    )
  )

ipsc_data <- pca_df %>%
  filter(CellType == "iPSC-iN") %>%
  arrange(PC2, PC1) %>%
  mutate(
    Side = ifelse(
      row_number() %% 2 == 0,
      "left",
      "right"
    )
  )

udc_left_x <- min(udc_data$PC1) -
  0.22 * pc1_span

udc_right_x <- max(udc_data$PC1) +
  0.12 * pc1_span

ipsc_left_x <- min(ipsc_data$PC1) -
  0.13 * pc1_span

ipsc_right_x <- max(ipsc_data$PC1) +
  0.18 * pc1_span


# ------------------------------------------------
# Helper: assign non-overlapping label Y positions
# ------------------------------------------------

spread_y <- function(y, min_gap) {

  ord <- order(y)
  result <- numeric(length(y))

  sorted <- y[ord]

  if (length(sorted) > 1) {
    for (i in 2:length(sorted)) {
      sorted[i] <- max(
        sorted[i],
        sorted[i - 1] + min_gap
      )
    }
  }

  # Recenter to preserve overall group location

  sorted <- sorted -
    mean(sorted) +
    mean(y)

  result[ord] <- sorted

  result
}


# ------------------------------------------------
# Assign label positions
# ------------------------------------------------

make_label_positions <- function(
    df,
    left_x,
    right_x,
    min_gap
) {

  left_df <- df %>%
    filter(Side == "left") %>%
    arrange(PC2)

  right_df <- df %>%
    filter(Side == "right") %>%
    arrange(PC2)

  if (nrow(left_df) > 0) {
    left_df$LabelX <- left_x
    left_df$LabelY <- spread_y(
      left_df$PC2,
      min_gap
    )
    left_df$Hjust <- 1
  }

  if (nrow(right_df) > 0) {
    right_df$LabelX <- right_x
    right_df$LabelY <- spread_y(
      right_df$PC2,
      min_gap
    )
    right_df$Hjust <- 0
  }

  bind_rows(
    left_df,
    right_df
  )
}

udc_labels <- make_label_positions(
  udc_data,
  udc_left_x,
  udc_right_x,
  label_step
)

ipsc_labels <- make_label_positions(
  ipsc_data,
  ipsc_left_x,
  ipsc_right_x,
  label_step
)

label_df <- bind_rows(
  udc_labels,
  ipsc_labels
)


# ================================================================
# 12. Group labels
# ================================================================

# Place cell-type labels above the corresponding
# HC/NPC ellipses and sample labels.

make_group_label <- function(cell_type) {

  sub_points <- pca_df %>%
    filter(CellType == cell_type)

  sub_ellipse <- ellipse_df %>%
    filter(CellType == cell_type)

  sub_labels <- label_df %>%
    filter(CellType == cell_type)

  top <- max(
    sub_points$PC2,
    sub_ellipse$PC2,
    sub_labels$LabelY
  )

  data.frame(
    CellType = cell_type,
    PC1 = mean(sub_points$PC1),
    PC2 = top + 0.13 * pc2_span
  )
}

group_labels <- bind_rows(
  make_group_label("UDC-iN"),
  make_group_label("iPSC-iN")
)

group_labels$CellType <- factor(
  group_labels$CellType,
  levels = c("UDC-iN", "iPSC-iN")
)


# ================================================================
# 13. Plot limits
# ================================================================

# Include points, ellipses, labels and group names.
# No clipping.

all_x <- c(
  pca_df$PC1,
  ellipse_df$PC1,
  label_df$LabelX,
  group_labels$PC1
)

all_y <- c(
  pca_df$PC2,
  ellipse_df$PC2,
  label_df$LabelY,
  group_labels$PC2
)

x_span <- diff(range(all_x))
y_span <- diff(range(all_y))

x_limits <- c(
  min(all_x) - 0.07 * x_span,
  max(all_x) + 0.07 * x_span
)

y_limits <- c(
  min(all_y) - 0.08 * y_span,
  max(all_y) + 0.08 * y_span
)


# ================================================================
# 14. Main Figure
# ================================================================

p_main <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2)
) +

  # HC/NPC ellipse fill
  geom_polygon(
    data = ellipse_df,
    aes(
      x = PC1,
      y = PC2,
      group = EllipseGroup,
      fill = CellType
    ),
    inherit.aes = FALSE,
    alpha = 0.10,
    color = NA
  ) +

  # HC/NPC ellipse outline
  geom_path(
    data = ellipse_df,
    aes(
      x = PC1,
      y = PC2,
      group = EllipseGroup,
      color = CellType,
      linetype = Diagnosis
    ),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +

  # Label leader lines
  geom_segment(
    data = label_df,
    aes(
      x = PC1,
      y = PC2,
      xend = LabelX,
      yend = LabelY
    ),
    inherit.aes = FALSE,
    linewidth = 0.22,
    color = "#C5C9CE"
  ) +

  # Sample points
  geom_point(
    aes(fill = CellType),
    shape = 21,
    size = 3.3,
    stroke = 0.55,
    color = "white"
  ) +

  # Labels
  geom_text(
    data = label_df,
    aes(
      x = LabelX,
      y = LabelY,
      label = Label,
      color = CellType,
      hjust = Hjust
    ),
    inherit.aes = FALSE,
    size = 2.8,
    family = "sans",
    show.legend = FALSE
  ) +

  # Group labels
  geom_text(
    data = group_labels,
    aes(
      x = PC1,
      y = PC2,
      label = CellType,
      color = CellType
    ),
    inherit.aes = FALSE,
    size = 4.6,
    fontface = "bold",
    family = "sans",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = cell_colors
  ) +

  scale_fill_manual(
    values = cell_colors
  ) +

  scale_linetype_manual(
    values = diagnosis_linetypes
  ) +

  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      variance[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      variance[2]
    )
  ) +

  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits,
    clip = "off"
  ) +

  theme_classic(base_size = 12) +

  theme(
    text = element_text(
      family = "sans",
      color = "#202020"
    ),

    axis.title = element_text(
      size = 13,
      color = "#202020"
    ),

    axis.text = element_text(
      size = 10.5,
      color = "#303030"
    ),

    axis.line = element_line(
      linewidth = 0.55,
      color = "#303030"
    ),

    axis.ticks = element_line(
      linewidth = 0.45,
      color = "#303030"
    ),

    legend.position = "none",
    panel.grid = element_blank(),

    plot.margin = margin(
      35, 45, 30, 45
    )
  )


# ================================================================
# 15. Save figures
# ================================================================

png_file <- file.path(
  outdir,
  "PCA_AllBatches_HC_NPC_600dpi.png"
)

ggsave(
  filename = png_file,
  plot = p_main,
  device = "png",
  width = 12.5,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

preview_file <- file.path(
  outdir,
  "PCA_AllBatches_HC_NPC_Preview_300dpi.png"
)

ggsave(
  filename = preview_file,
  plot = p_main,
  device = "png",
  width = 12.5,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ================================================================
# 16. Export results
# ================================================================

write.csv(
  pca_df,
  file.path(outdir, "PCA_Coordinates.csv"),
  row.names = FALSE
)

write.csv(
  label_df,
  file.path(outdir, "PCA_Label_Positions.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    PC = paste0("PC", seq_along(variance)),
    VariancePercent = variance
  ),
  file.path(outdir, "PCA_Variance.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Gene = colnames(pca_input)
  ),
  file.path(outdir, "PCA_GeneList.csv"),
  row.names = FALSE
)

write.csv(
  vst_mat,
  file.path(outdir, "RNAseq_VST_AllGenes.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_SessionInfo.txt")
)


# ================================================================
# 17. Display
# ================================================================

print(p_main)

if (Sys.info()[["sysname"]] == "Darwin") {
  system2(
    "open",
    args = preview_file,
    wait = FALSE
  )
}


# ================================================================
# 18. Summary
# ================================================================

cat("\n====================================\n")
cat("PCA COMPLETED | ALL BATCHES\n")
cat("====================================\n")

cat("Samples:", nrow(pca_df), "\n")
cat("Donors:", nlevels(meta$Donor), "\n")
cat("Genes:", ncol(pca_input), "\n")

cat(
  "PC1:",
  round(variance[1], 1),
  "%\n"
)

cat(
  "PC2:",
  round(variance[2], 1),
  "%\n"
)

cat("\nSample counts:\n")

print(
  table(
    pca_df$CellType,
    pca_df$Diagnosis
  )
)

cat("\nSaved:\n", png_file, "\n")
cat("====================================\n")
