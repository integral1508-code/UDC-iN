# ================================================================
# RNA-seq PCA | UDC-iN vs UDC
# Paired RNA-seq PCA and figure generation
#
# Changes:
# - Avoid overlap between UDC / UDC-iN and donor labels
# - Preserve PCA coordinates, colors and 95% ellipses
# - PNG 600 dpi
#
# Input files:
# RNAseq_gene_count_matrix_22samples_GEO.csv
# RNA-seq_Donor information.xlsx
# ================================================================

# ================================================================
# 0. Packages
# ================================================================

cran_packages <- c(
  "readxl",
  "ggplot2",
  "ggrepel",
  "dplyr",
  "tibble"
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
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
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
library(ggrepel)
library(dplyr)
library(tibble)

set.seed(123)

# ================================================================
# 1. Paths
# ================================================================

input_dir <- Sys.getenv("UDC_INPUT_DIR", unset = "data")
output_dir <- Sys.getenv("UDC_OUTPUT_DIR", unset = "results")

count_file <- file.path(
  input_dir,
  "RNAseq_gene_count_matrix_22samples_GEO.csv"
)

meta_file <- file.path(
  input_dir,
  "RNA-seq_Donor information.xlsx"
)

outdir <- file.path(
  output_dir,
  "RNAseq_PCA_Nature_MainFigure"
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
# 2. Count matrix
# ================================================================

counts_raw <- read.csv(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("Gene" %in% names(counts_raw))

counts_raw <- counts_raw[
  !is.na(counts_raw$Gene) &
    trimws(as.character(counts_raw$Gene)) != "",
  ,
  drop = FALSE
]

sample_cols <- setdiff(
  names(counts_raw),
  "Gene"
)

counts_raw[sample_cols] <- lapply(
  counts_raw[sample_cols],
  function(x) as.numeric(as.character(x))
)

count_values <- as.matrix(
  counts_raw[, sample_cols, drop = FALSE]
)

if (
  anyNA(count_values) ||
  any(!is.finite(count_values)) ||
  any(count_values < 0) ||
  any(abs(count_values - round(count_values)) > 1e-8)
) {
  stop("Counts must be nonnegative integers.")
}

counts <- counts_raw %>%
  group_by(Gene) %>%
  summarise(
    across(
      all_of(sample_cols),
      sum
    ),
    .groups = "drop"
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

if (any(counts > .Machine$integer.max)) {
  stop("Count values exceed integer limit.")
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
  "Cell type"
)

stopifnot(
  all(required_cols %in% names(meta_raw))
)

meta <- meta_raw %>%
  transmute(
    Sample = trimws(as.character(`FASTQ name`)),
    Donor = trimws(as.character(`Donor ID`)),
    Age = as.numeric(Age),
    Diagnosis = trimws(as.character(Condition)),
    CellType = trimws(as.character(`Cell type`))
  )

stopifnot(!anyDuplicated(meta$Sample))

# ================================================================
# 4. Match samples
# ================================================================

if (!setequal(meta$Sample, colnames(counts))) {
  stop("Sample names do not match.")
}

meta <- meta[
  match(colnames(counts), meta$Sample),
  ,
  drop = FALSE
]

stopifnot(
  identical(meta$Sample, colnames(counts))
)

meta$CellType <- factor(
  meta$CellType,
  levels = c("UDC", "UDC-iN")
)

meta$Donor <- factor(meta$Donor)

if (
  anyNA(meta$Age) ||
  anyNA(meta$Diagnosis) ||
  anyNA(meta$Donor) ||
  anyNA(meta$CellType)
) {
  stop("Invalid metadata.")
}

meta$Label <- paste0(
  meta$Diagnosis,
  ". ",
  format(meta$Age, trim = TRUE),
  "y"
)

# ================================================================
# 5. Paired donor validation
# ================================================================

pair_table <- table(
  meta$Donor,
  meta$CellType
)

if (!all(pair_table == 1)) {
  stop("Each donor must have one UDC and one UDC-iN.")
}

donor_check <- meta %>%
  group_by(Donor) %>%
  summarise(
    n_age = n_distinct(Age),
    n_diagnosis = n_distinct(Diagnosis),
    .groups = "drop"
  )

if (
  any(donor_check$n_age != 1) ||
  any(donor_check$n_diagnosis != 1)
) {
  stop("Age or diagnosis differs within donor pairs.")
}

# ================================================================
# 6. Low-expression filtering
# ================================================================

# Raw count >=10 in >=3 samples

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

# ================================================================
# 7. VST normalization
# ================================================================

coldata <- meta
rownames(coldata) <- coldata$Sample

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ Donor + CellType
)

dds <- estimateSizeFactors(dds)

vsd <- varianceStabilizingTransformation(
  dds,
  blind = TRUE
)

vst_mat <- assay(vsd)

# ================================================================
# 8. PCA using ALL expressed genes
# ================================================================

pca_input <- t(vst_mat)

gene_var <- apply(
  pca_input,
  2,
  var
)

valid_genes <- is.finite(gene_var) & gene_var > 0

pca_input <- pca_input[
  ,
  valid_genes,
  drop = FALSE
]

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
# 9. 95% data ellipses
# ================================================================

make_ellipse <- function(df, level = 0.95, n = 300) {

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
    stop("Ellipse covariance is singular.")
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
    diag(sqrt(eig$values), nrow = 2)

  radius <- sqrt(qchisq(level, df = 2))

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

ellipse_df <- bind_rows(
  lapply(
    levels(meta$CellType),
    function(ct) {

      tmp <- pca_df[
        pca_df$CellType == ct,
        ,
        drop = FALSE
      ]

      ell <- make_ellipse(tmp)
      ell$CellType <- ct
      ell
    }
  )
)

ellipse_df$CellType <- factor(
  ellipse_df$CellType,
  levels = c("UDC", "UDC-iN")
)

# ================================================================
# 10. FIXED group label positions
# ================================================================

# Group labels are positioned in clear regions
# above the main sample clusters.
#
# These positions do not change PCA coordinates.

group_labels <- data.frame(
  CellType = factor(
    c("UDC", "UDC-iN"),
    levels = c("UDC", "UDC-iN")
  ),
  PC1 = c(-113, 86),
  PC2 = c(155, 64)
)

# ================================================================
# 11. Colors
# ================================================================

cell_colors <- c(
  "UDC" = "#3C78A6",
  "UDC-iN" = "#CC6859"
)

# ================================================================
# 12. Main Figure
# ================================================================

# Use ggrepel with invisible obstacle labels.
# This helps donor labels avoid the group names.

obstacles <- data.frame(
  PC1 = group_labels$PC1,
  PC2 = group_labels$PC2,
  Label = group_labels$CellType
)

label_data <- rbind(
  data.frame(
    PC1 = pca_df$PC1,
    PC2 = pca_df$PC2,
    Label = pca_df$Label,
    CellType = as.character(pca_df$CellType),
    IsGroup = FALSE
  ),
  data.frame(
    PC1 = obstacles$PC1,
    PC2 = obstacles$PC2,
    Label = as.character(obstacles$Label),
    CellType = as.character(obstacles$Label),
    IsGroup = TRUE
  )
)

label_data$CellType <- factor(
  label_data$CellType,
  levels = c("UDC", "UDC-iN")
)

# Group names are drawn separately.
# Empty labels still participate as repulsion points.
label_data$PlotLabel <- ifelse(
  label_data$IsGroup,
  "",
  label_data$Label
)

p_main <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2)
) +

  # Paired donor connections
  geom_line(
    aes(group = Donor),
    color = "#C9CDD1",
    linewidth = 0.28,
    alpha = 0.45
  ) +

  # Ellipse fill
  geom_polygon(
    data = ellipse_df,
    aes(
      x = PC1,
      y = PC2,
      group = CellType,
      fill = CellType
    ),
    inherit.aes = FALSE,
    alpha = 0.045,
    color = NA
  ) +

  # Ellipse outline
  geom_path(
    data = ellipse_df,
    aes(
      x = PC1,
      y = PC2,
      group = CellType,
      color = CellType
    ),
    inherit.aes = FALSE,
    linewidth = 0.75
  ) +

  # Sample points
  geom_point(
    aes(fill = CellType),
    shape = 21,
    size = 3.3,
    stroke = 0.55,
    color = "white"
  ) +

  # Donor labels with group obstacles
  geom_text_repel(
    data = label_data,
    aes(
      x = PC1,
      y = PC2,
      label = PlotLabel,
      color = CellType
    ),
    inherit.aes = FALSE,
    size = 3.05,
    family = "sans",
    box.padding = 0.5,
    point.padding = 0.45,
    force = 3,
    force_pull = 0.4,
    max.overlaps = Inf,
    max.iter = 15000,
    max.time = 8,
    min.segment.length = 0,
    segment.color = "#B5B9BE",
    segment.size = 0.22,
    seed = 123,
    show.legend = FALSE
  ) +

  # Group labels drawn at fixed positions
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

  scale_color_manual(values = cell_colors) +
  scale_fill_manual(values = cell_colors) +

  labs(
    x = sprintf("PC1 (%.1f%%)", variance[1]),
    y = sprintf("PC2 (%.1f%%)", variance[2])
  ) +

  scale_x_continuous(
    expand = expansion(mult = 0.10)
  ) +

  scale_y_continuous(
    expand = expansion(mult = 0.12)
  ) +

  coord_cartesian(clip = "off") +

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
      25, 35, 25, 25
    )
  )

# ================================================================
# 13. Save PNG
# ================================================================

png_file <- file.path(
  outdir,
  "PCA_AllGenes_Nature_MainFigure_600dpi.png"
)

ggsave(
  filename = png_file,
  plot = p_main,
  device = "png",
  width = 8.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ================================================================
# 14. Preview PNG
# ================================================================

preview_file <- file.path(
  outdir,
  "PCA_AllGenes_Nature_Preview_300dpi.png"
)

ggsave(
  filename = preview_file,
  plot = p_main,
  device = "png",
  width = 8.5,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

# Display in RStudio
if (interactive()) print(p_main)


# ================================================================
# 15. Export PCA results
# ================================================================

write.csv(
  pca_df,
  file.path(outdir, "PCA_Coordinates.csv"),
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
  data.frame(Gene = colnames(pca_input)),
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
# 16. Summary
# ================================================================

cat("\n====================================\n")
cat("PCA COMPLETED\n")
cat("====================================\n")

cat("Samples:", nrow(pca_df), "\n")
cat("Donors:", nlevels(meta$Donor), "\n")
cat("Genes:", ncol(pca_input), "\n")
cat("PC1:", round(variance[1], 1), "%\n")
cat("PC2:", round(variance[2], 1), "%\n")

cat("\nPNG saved:\n", png_file, "\n")
cat("====================================\n")
