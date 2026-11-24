# ================================================================
# UDC-iN vs UDC | RNA-seq DEG + GO Biological Process
# Paired differential expression and GO enrichment
#
# Samples: 22 (11 paired donors)
#
# DESeq2:
#   design = ~ Donor + CellType
#   contrast = UDC-iN / UDC
#   padj < 0.05 & abs(log2FC) > 1
#
# GO:
#   Biological Process ORA
#   UP / DOWN analyzed separately
#   Top 10 significant terms
#
# Figures:
#   01 Volcano
#   02 GO UP Bar
#   03 GO DOWN Bar
#   04 GO UP Dot
#   05 GO DOWN Dot
#
# Output: PNG 600 dpi
# No GSEA / No svglite required
# ================================================================


# ================================================================
# 0. PACKAGES
# ================================================================

cran_pkgs <- c(
  "readxl",
  "dplyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "stringr"
)

bioc_pkgs <- c(
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

library(DESeq2)
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(123)


# ================================================================
# 1. PATHS AND PARAMETERS
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
  "RNAseq_UDCiN_MainFigure_DEG_GO"
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

padj_cutoff <- 0.05
lfc_cutoff <- 1
go_padj_cutoff <- 0.05
top_n <- 10

color_up <- "#CC6859"
color_down <- "#3C78A6"
color_ns <- "#C9CDD1"


# ================================================================
# 2. SAVE PNG FUNCTION
# ================================================================

save_png <- function(
    plot,
    filename,
    width = 8,
    height = 6
) {

  path <- file.path(outdir, filename)

  ggplot2::ggsave(
    filename = path,
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  cat("Saved:", path, "\n")
}


# ================================================================
# 3. READ COUNT MATRIX
# ================================================================

counts_raw <- read.csv(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Gene" %in% names(counts_raw)) {
  stop("Gene column not found.")
}

counts_raw <- counts_raw[
  !is.na(counts_raw$Gene) &
    trimws(as.character(counts_raw$Gene)) != "",
  ,
  drop = FALSE
]

counts_raw$Gene <- trimws(
  as.character(counts_raw$Gene)
)

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
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_cols),
      sum
    ),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

if (any(counts > .Machine$integer.max)) {
  stop("Count values exceed integer limit.")
}

storage.mode(counts) <- "integer"

cat("Original genes:", nrow(counts), "\n")


# ================================================================
# 4. READ METADATA
# ================================================================

meta_raw <- readxl::read_excel(meta_file)

required_cols <- c(
  "FASTQ name",
  "Donor ID",
  "Age",
  "Condition",
  "Cell type"
)

missing_cols <- setdiff(
  required_cols,
  names(meta_raw)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

meta <- meta_raw %>%
  dplyr::transmute(
    Sample = trimws(as.character(`FASTQ name`)),
    Donor = trimws(as.character(`Donor ID`)),
    Age = as.numeric(Age),
    Diagnosis = trimws(as.character(Condition)),
    CellType = trimws(as.character(`Cell type`))
  )

if (anyDuplicated(meta$Sample)) {
  stop("Duplicated sample IDs.")
}

if (!setequal(meta$Sample, colnames(counts))) {

  cat("\nCount-only samples:\n")
  print(setdiff(colnames(counts), meta$Sample))

  cat("\nMetadata-only samples:\n")
  print(setdiff(meta$Sample, colnames(counts)))

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

meta$Donor <- factor(meta$Donor)

meta$CellType <- factor(
  meta$CellType,
  levels = c("UDC", "UDC-iN")
)

if (
  anyNA(meta$Donor) ||
  anyNA(meta$CellType) ||
  anyNA(meta$Age) ||
  anyNA(meta$Diagnosis)
) {
  stop("Invalid metadata.")
}


# ================================================================
# 5. VALIDATE PAIRED DESIGN
# ================================================================

pair_table <- table(
  meta$Donor,
  meta$CellType
)

if (!all(pair_table == 1)) {
  print(pair_table)
  stop("Incomplete donor pairing.")
}

donor_check <- meta %>%
  dplyr::group_by(Donor) %>%
  dplyr::summarise(
    n_age = dplyr::n_distinct(Age),
    n_diagnosis = dplyr::n_distinct(Diagnosis),
    .groups = "drop"
  )

if (
  any(donor_check$n_age != 1) ||
  any(donor_check$n_diagnosis != 1)
) {
  stop("Donor metadata mismatch.")
}

cat("\nSamples:", nrow(meta), "\n")
cat("Donors:", nlevels(meta$Donor), "\n")


# ================================================================
# 6. LOW-EXPRESSION FILTERING
# ================================================================

# Raw count >=10 in >=3 samples

keep <- rowSums(counts >= 10) >= 3

counts_filtered <- counts[
  keep,
  ,
  drop = FALSE
]

if (nrow(counts_filtered) < 2) {
  stop("Too few genes after filtering.")
}

cat(
  "Genes retained:",
  nrow(counts_filtered),
  "\n"
)


# ================================================================
# 7. PAIRED DESEQ2 ANALYSIS
# ================================================================

coldata <- as.data.frame(meta)
rownames(coldata) <- coldata$Sample

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ Donor + CellType
)

dds <- DESeq2::DESeq(
  dds,
  test = "Wald"
)

res <- DESeq2::results(
  dds,
  contrast = c(
    "CellType",
    "UDC-iN",
    "UDC"
  ),
  alpha = 0.05
)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("Gene")


# ================================================================
# 8. DEG CLASSIFICATION
# ================================================================

res_df <- res_df %>%
  dplyr::mutate(
    Regulation = dplyr::case_when(

      !is.na(padj) &
        padj < padj_cutoff &
        log2FoldChange > lfc_cutoff ~ "UP",

      !is.na(padj) &
        padj < padj_cutoff &
        log2FoldChange < -lfc_cutoff ~ "DOWN",

      TRUE ~ "NS"
    ),

    negLog10Padj = -log10(
      pmax(
        padj,
        .Machine$double.xmin
      )
    )
  )

deg_up <- res_df %>%
  dplyr::filter(Regulation == "UP")

deg_down <- res_df %>%
  dplyr::filter(Regulation == "DOWN")

deg_all <- res_df %>%
  dplyr::filter(Regulation != "NS")

cat("\n====================================\n")
cat("DEG SUMMARY\n")
cat("====================================\n")

cat("UP:", nrow(deg_up), "\n")
cat("DOWN:", nrow(deg_down), "\n")
cat("TOTAL:", nrow(deg_all), "\n")


# ================================================================
# 9. ROBUST GENE ID MAPPING
# ================================================================

# Gene identifier mapping:
# Always create SYMBOL and ENTREZID columns explicitly.
# Avoid renaming SYMBOL to InputID.

gene_ids <- res_df$Gene

ensembl_fraction <- mean(
  grepl("^ENSG[0-9]+", gene_ids)
)

if (ensembl_fraction > 0.5) {

  id_type <- "ENSEMBL"

  clean_ids <- sub(
    "\\.[0-9]+$",
    "",
    gene_ids
  )

} else {

  id_type <- "SYMBOL"
  clean_ids <- gene_ids
}

cat("\nDetected ID type:", id_type, "\n")

gene_map_input <- data.frame(
  Gene = gene_ids,
  InputID = clean_ids,
  stringsAsFactors = FALSE
)

# Map Entrez IDs
entrez_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(clean_ids),
  column = "ENTREZID",
  keytype = id_type,
  multiVals = "first"
)

# Map gene symbols
symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(clean_ids),
  column = "SYMBOL",
  keytype = id_type,
  multiVals = "first"
)

# Construct mapping explicitly
gene_map <- gene_map_input %>%
  dplyr::mutate(
    ENTREZID = unname(
      entrez_map[InputID]
    ),
    SYMBOL = unname(
      symbol_map[InputID]
    )
  )

# For SYMBOL input, preserve original gene symbol
if (id_type == "SYMBOL") {
  gene_map$SYMBOL <- gene_map$InputID
}

gene_map$ENTREZID <- as.character(
  gene_map$ENTREZID
)

gene_map$SYMBOL <- as.character(
  gene_map$SYMBOL
)

# Export all mapping results, including unmapped genes
write.csv(
  gene_map,
  file.path(
    outdir,
    "Gene_ID_Mapping_All.csv"
  ),
  row.names = FALSE
)

# GO mapping
gene_map_go <- gene_map %>%
  dplyr::filter(
    !is.na(ENTREZID),
    ENTREZID != ""
  ) %>%
  dplyr::distinct(
    Gene,
    ENTREZID,
    .keep_all = TRUE
  )

cat(
  "Mapped genes:",
  nrow(gene_map_go),
  "\n"
)


# ================================================================
# 10. EXPORT DEG
# ================================================================

res_df_export <- res_df %>%
  dplyr::left_join(
    gene_map %>%
      dplyr::select(
        Gene,
        SYMBOL,
        ENTREZID
      ),
    by = "Gene"
  )

write.csv(
  res_df_export,
  file.path(
    outdir,
    "DESeq2_AllGenes.csv"
  ),
  row.names = FALSE
)

write.csv(
  res_df_export %>%
    dplyr::filter(Regulation != "NS"),
  file.path(
    outdir,
    "DESeq2_DEG_All.csv"
  ),
  row.names = FALSE
)

write.csv(
  res_df_export %>%
    dplyr::filter(Regulation == "UP"),
  file.path(
    outdir,
    "DESeq2_DEG_UP.csv"
  ),
  row.names = FALSE
)

write.csv(
  res_df_export %>%
    dplyr::filter(Regulation == "DOWN"),
  file.path(
    outdir,
    "DESeq2_DEG_DOWN.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 11. VOLCANO PLOT
# ================================================================

volcano_df <- res_df_export %>%
  dplyr::filter(
    is.finite(log2FoldChange),
    !is.na(padj),
    is.finite(negLog10Padj)
  )

if (nrow(volcano_df) == 0) {
  stop("No valid genes for volcano plot.")
}

label_genes <- volcano_df %>%
  dplyr::filter(Regulation != "NS") %>%
  dplyr::group_by(Regulation) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs(log2FoldChange)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    DisplayGene = dplyr::if_else(
      !is.na(SYMBOL) & SYMBOL != "",
      SYMBOL,
      Gene
    )
  )

volcano_colors <- c(
  "UP" = color_up,
  "DOWN" = color_down,
  "NS" = color_ns
)

xmax <- max(
  2,
  ceiling(
    max(abs(volcano_df$log2FoldChange))
  )
)

p_volcano <- ggplot2::ggplot(
  volcano_df,
  ggplot2::aes(
    x = log2FoldChange,
    y = negLog10Padj
  )
) +

  ggplot2::geom_point(
    data = subset(
      volcano_df,
      Regulation == "NS"
    ),
    color = color_ns,
    size = 0.85,
    alpha = 0.4
  ) +

  ggplot2::geom_point(
    data = subset(
      volcano_df,
      Regulation != "NS"
    ),
    ggplot2::aes(color = Regulation),
    size = 1.15,
    alpha = 0.8
  ) +

  ggplot2::geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey55"
  ) +

  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey55"
  ) +

  ggrepel::geom_text_repel(
    data = label_genes,
    ggplot2::aes(
      label = DisplayGene,
      color = Regulation
    ),
    size = 3.0,
    max.overlaps = Inf,
    max.iter = 10000,
    max.time = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey65",
    segment.size = 0.25,
    seed = 123,
    show.legend = FALSE
  ) +

  ggplot2::scale_color_manual(
    values = volcano_colors
  ) +

  ggplot2::scale_x_continuous(
    limits = c(-xmax, xmax),
    expand = ggplot2::expansion(mult = 0.04)
  ) +

  ggplot2::labs(
    x = expression(
      log[2]~"fold change (UDC-iN / UDC)"
    ),
    y = expression(
      -log[10]~"(adjusted P)"
    )
  ) +

  ggplot2::theme_classic(base_size = 12) +

  ggplot2::theme(
    axis.title = ggplot2::element_text(
      size = 13,
      color = "black"
    ),
    axis.text = ggplot2::element_text(
      size = 10.5,
      color = "black"
    ),
    axis.line = ggplot2::element_line(
      linewidth = 0.55
    ),
    axis.ticks = ggplot2::element_line(
      linewidth = 0.45
    ),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      20, 25, 20, 20
    )
  )

save_png(
  p_volcano,
  "01_Volcano_UDCiN_vs_UDC.png",
  width = 8,
  height = 6.5
)


# ================================================================
# 12. GO ORA BACKGROUND
# ================================================================

tested_genes <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::pull(Gene)

universe_ids <- gene_map_go %>%
  dplyr::filter(Gene %in% tested_genes) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

up_ids <- gene_map_go %>%
  dplyr::filter(Gene %in% deg_up$Gene) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

down_ids <- gene_map_go %>%
  dplyr::filter(Gene %in% deg_down$Gene) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

up_ids <- intersect(
  up_ids,
  universe_ids
)

down_ids <- intersect(
  down_ids,
  universe_ids
)

cat("\nGO universe:", length(universe_ids), "\n")
cat("GO UP genes:", length(up_ids), "\n")
cat("GO DOWN genes:", length(down_ids), "\n")


# ================================================================
# 13. GO BIOLOGICAL PROCESS ORA
# ================================================================

run_go <- function(ids) {

  if (length(ids) < 5) {
    return(NULL)
  }

  clusterProfiler::enrichGO(
    gene = ids,
    universe = universe_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
}

go_up <- run_go(up_ids)
go_down <- run_go(down_ids)


# ================================================================
# 14. EXTRACT GO RESULTS
# ================================================================

extract_go <- function(object, direction) {

  if (is.null(object)) {
    return(data.frame(p.adjust = numeric(0)))
  }

  df <- as.data.frame(object)

  if (nrow(df) == 0) {
    return(data.frame(p.adjust = numeric(0)))
  }

  df %>%
    dplyr::mutate(
      Direction = direction,
      negLog10Padj = -log10(
        pmax(
          p.adjust,
          .Machine$double.xmin
        )
      )
    )
}

go_up_all <- extract_go(go_up, "UP")
go_down_all <- extract_go(go_down, "DOWN")

go_up_sig <- go_up_all %>%
  dplyr::filter(p.adjust < go_padj_cutoff)

go_down_sig <- go_down_all %>%
  dplyr::filter(p.adjust < go_padj_cutoff)

write.csv(
  go_up_all,
  file.path(outdir, "GO_BP_UP_All.csv"),
  row.names = FALSE
)

write.csv(
  go_down_all,
  file.path(outdir, "GO_BP_DOWN_All.csv"),
  row.names = FALSE
)

write.csv(
  go_up_sig,
  file.path(outdir, "GO_BP_UP_Significant.csv"),
  row.names = FALSE
)

write.csv(
  go_down_sig,
  file.path(outdir, "GO_BP_DOWN_Significant.csv"),
  row.names = FALSE
)


# ================================================================
# 15. SELECT TOP 10 GO TERMS
# ================================================================

select_top_go <- function(df) {

  if (nrow(df) == 0) {
    return(df)
  }

  df %>%
    dplyr::arrange(
      p.adjust,
      dplyr::desc(Count),
      ID
    ) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Term = stringr::str_wrap(
        Description,
        width = 43
      ),
      negLog10Padj = -log10(
        pmax(
          p.adjust,
          .Machine$double.xmin
        )
      )
    )
}

go_up_top <- select_top_go(go_up_sig)
go_down_top <- select_top_go(go_down_sig)

write.csv(
  go_up_top,
  file.path(outdir, "GO_BP_UP_Top10.csv"),
  row.names = FALSE
)

write.csv(
  go_down_top,
  file.path(outdir, "GO_BP_DOWN_Top10.csv"),
  row.names = FALSE
)


# ================================================================
# 16. GO BAR PLOT
# ================================================================

make_go_bar <- function(
    df,
    direction,
    color,
    filename
) {

  if (nrow(df) == 0) {
    message("No significant GO terms: ", direction)
    return(NULL)
  }

  df$TermID <- factor(
    df$ID,
    levels = rev(df$ID)
  )

  term_labels <- setNames(
    df$Term,
    df$ID
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = negLog10Padj,
      y = TermID
    )
  ) +

    ggplot2::geom_col(
      fill = color,
      width = 0.72,
      alpha = 0.92
    ) +

    ggplot2::scale_y_discrete(
      labels = term_labels
    ) +

    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(0, 0.08)
      )
    ) +

    ggplot2::labs(
      x = expression(
        -log[10]~"(BH-adjusted P)"
      ),
      y = NULL
    ) +

    ggplot2::theme_classic(base_size = 12) +

    ggplot2::theme(
      axis.title.x = ggplot2::element_text(
        size = 12,
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = 10,
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10.5,
        color = "black"
      ),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(
        linewidth = 0.55
      ),
      plot.margin = ggplot2::margin(
        15, 20, 15, 15
      )
    )

  save_png(
    p,
    filename,
    width = 10,
    height = 5.8
  )

  return(p)
}

p_go_up_bar <- make_go_bar(
  go_up_top,
  "UP",
  color_up,
  "02_GO_BP_UP_Top10_Bar.png"
)

p_go_down_bar <- make_go_bar(
  go_down_top,
  "DOWN",
  color_down,
  "03_GO_BP_DOWN_Top10_Bar.png"
)


# ================================================================
# 17. GENE RATIO CONVERSION
# ================================================================

ratio_to_numeric <- function(x) {

  vapply(
    strsplit(as.character(x), "/"),
    function(z) {
      as.numeric(z[1]) / as.numeric(z[2])
    },
    numeric(1)
  )
}

if (nrow(go_up_top) > 0) {
  go_up_top$GeneRatioNumeric <- ratio_to_numeric(
    go_up_top$GeneRatio
  )
}

if (nrow(go_down_top) > 0) {
  go_down_top$GeneRatioNumeric <- ratio_to_numeric(
    go_down_top$GeneRatio
  )
}


# ================================================================
# 18. GO DOT PLOT
# ================================================================

make_go_dot <- function(
    df,
    direction,
    filename
) {

  if (nrow(df) == 0) {
    message("No significant GO terms: ", direction)
    return(NULL)
  }

  df$TermID <- factor(
    df$ID,
    levels = rev(df$ID)
  )

  term_labels <- setNames(
    df$Term,
    df$ID
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = GeneRatioNumeric,
      y = TermID
    )
  ) +

    ggplot2::geom_point(
      ggplot2::aes(
        size = Count,
        color = p.adjust
      ),
      alpha = 0.95
    ) +

    ggplot2::scale_y_discrete(
      labels = term_labels
    ) +

    ggplot2::scale_size_continuous(
      range = c(3.5, 9),
      name = "Gene count"
    ) +

    ggplot2::scale_color_gradient(
      low = "#B9473D",
      high = "#F0C6B8",
      limits = c(0, go_padj_cutoff),
      name = "Adjusted P"
    ) +

    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(0.05, 0.12)
      )
    ) +

    ggplot2::labs(
      x = "Gene ratio",
      y = NULL
    ) +

    ggplot2::theme_classic(base_size = 12) +

    ggplot2::theme(
      axis.title.x = ggplot2::element_text(
        size = 12,
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = 10,
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10.5,
        color = "black"
      ),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(
        linewidth = 0.55
      ),
      legend.title = ggplot2::element_text(
        size = 10
      ),
      legend.text = ggplot2::element_text(
        size = 9
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(
        15, 20, 15, 15
      )
    )

  save_png(
    p,
    filename,
    width = 10,
    height = 6
  )

  return(p)
}

p_go_up_dot <- make_go_dot(
  go_up_top,
  "UP",
  "04_GO_BP_UP_Top10_Dot.png"
)

p_go_down_dot <- make_go_dot(
  go_down_top,
  "DOWN",
  "05_GO_BP_DOWN_Top10_Dot.png"
)


# ================================================================
# 19. NORMALIZED COUNTS
# ================================================================

normalized_counts <- DESeq2::counts(
  dds,
  normalized = TRUE
)

write.csv(
  normalized_counts,
  file.path(
    outdir,
    "DESeq2_Normalized_Counts.csv"
  )
)


# ================================================================
# 20. SUMMARY TABLE
# ================================================================

summary_df <- data.frame(
  Metric = c(
    "Samples",
    "Donors",
    "Genes after filtering",
    "Genes with non-NA padj",
    "DEG UP",
    "DEG DOWN",
    "Total DEG",
    "GO UP significant",
    "GO DOWN significant"
  ),

  Value = c(
    nrow(meta),
    nlevels(meta$Donor),
    nrow(counts_filtered),
    sum(!is.na(res_df$padj)),
    nrow(deg_up),
    nrow(deg_down),
    nrow(deg_all),
    nrow(go_up_sig),
    nrow(go_down_sig)
  )
)

write.csv(
  summary_df,
  file.path(
    outdir,
    "Analysis_Summary.csv"
  ),
  row.names = FALSE
)

print(summary_df)


# ================================================================
# 21. ANALYSIS SETTINGS
# ================================================================

settings <- c(
  "UDC-iN vs UDC RNA-seq",
  "Paired design: ~ Donor + CellType",
  "Donor included as a fixed effect",
  "Age and diagnosis constant within donor",
  "Contrast: UDC-iN / UDC",
  "Low-expression filter: count >=10 in >=3 samples",
  "DEG: padj <0.05 and abs(log2FC)>1",
  "DEG test: DESeq2 Wald test",
  "P adjustment: BH",
  "GO: Biological Process ORA",
  "GO UP and DOWN analyzed separately",
  "GO background: mapped genes with non-NA padj",
  "GO significance: BH-adjusted P <0.05",
  "GO Top 10 selected by adjusted P",
  "GO bar X-axis: -log10(adjusted P)",
  "GO dot X-axis: GeneRatio",
  "GO dot color: adjusted P",
  "GO dot size: DEG gene count",
  "No GSEA",
  "All figures: PNG 600 dpi"
)

writeLines(
  settings,
  file.path(
    outdir,
    "Analysis_Settings.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "R_SessionInfo.txt"
  )
)


# ================================================================
# 22. DISPLAY FIGURES
# ================================================================

if (interactive()) print(p_volcano)

if (interactive() && !is.null(p_go_up_bar)) {
  print(p_go_up_bar)
}

if (interactive() && !is.null(p_go_down_bar)) {
  print(p_go_down_bar)
}

if (interactive() && !is.null(p_go_up_dot)) {
  print(p_go_up_dot)
}

if (interactive() && !is.null(p_go_down_dot)) {
  print(p_go_down_dot)
}





# ================================================================
# 24. FINAL SUMMARY
# ================================================================

cat("\n====================================\n")
cat("DEG + GO BP ANALYSIS COMPLETED\n")
cat("====================================\n")

cat("Samples:", nrow(meta), "\n")
cat("Donors:", nlevels(meta$Donor), "\n")
cat("Genes:", nrow(counts_filtered), "\n")

cat("DEG UP:", nrow(deg_up), "\n")
cat("DEG DOWN:", nrow(deg_down), "\n")
cat("Total DEG:", nrow(deg_all), "\n")

cat("GO UP:", nrow(go_up_sig), "\n")
cat("GO DOWN:", nrow(go_down_sig), "\n")

cat("\nOutput directory:\n")
cat(outdir, "\n")

cat("====================================\n")
