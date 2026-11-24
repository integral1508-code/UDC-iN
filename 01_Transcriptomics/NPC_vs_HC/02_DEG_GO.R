# ============================================================
# RNA-seq | NPC vs HC
# UDC-iN / iPSC-iN
#
# STABLE FULL ANALYSIS
# TMM + voom + duplicateCorrelation + limma
#
# No dream
# No voomWithDreamWeights
# No parallel processing
#
# Input (Desktop):
# RNAseq_gene_count_matrix_26samples_GEO.csv
# RNA-seq_Donor information_2.xlsx
#
# Output:
# RNAseq_HC_NPC_Stable_Final
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

cran_packages <- c(
  "readxl",
  "dplyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "stringr",
  "scales"
)

bioc_packages <- c(
  "edgeR",
  "limma",
  "BiocParallel",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

for (pkg in cran_packages) {

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

for (pkg in bioc_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )

  }

}

library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(scales)
library(edgeR)
library(limma)
library(BiocParallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Prevent oversubscription
options(mc.cores = 1)

BiocParallel::register(
  BiocParallel::SerialParam()
)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

set.seed(123)


# ============================================================
# 1. SETTINGS
# ============================================================

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
  "RNAseq_HC_NPC_Stable_Final"
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
log2fc_cutoff <- 1

go_padj_cutoff <- 0.05

min_count <- 10
min_samples <- 3

top_n <- 10

color_up <- "#CC6859"
color_down <- "#3C78A6"
color_ns <- "#C9CDD1"

cell_types <- c(
  "UDC-iN",
  "iPSC-iN"
)


# ============================================================
# 2. LOGGING
# ============================================================

log_file <- file.path(
  outdir,
  "Progress_Log.txt"
)

log_step <- function(x) {

  msg <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    x
  )

  cat(msg, "\n")

  cat(
    msg,
    "\n",
    file = log_file,
    append = TRUE
  )

  flush.console()

}

log_step("ANALYSIS STARTED")


# ============================================================
# 3. READ COUNT MATRIX
# ============================================================

log_step("Reading counts")

raw <- read.csv(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  "Gene" %in% colnames(raw)
)

raw$Gene <- trimws(
  as.character(raw$Gene)
)

raw <- raw[
  !is.na(raw$Gene) &
    raw$Gene != "",
  ,
  drop = FALSE
]

sample_cols <- setdiff(
  colnames(raw),
  "Gene"
)

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
  stop("Count matrix contains invalid raw counts")
}

counts <- raw %>%
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

storage.mode(counts) <- "numeric"

log_step(
  paste(
    "Count matrix:",
    nrow(counts),
    "genes x",
    ncol(counts),
    "samples"
  )
)


# ============================================================
# 4. READ METADATA
# ============================================================

log_step("Reading metadata")

meta_raw <- readxl::read_excel(
  meta_file
)

required_columns <- c(
  "FASTQ name",
  "Donor ID",
  "Age",
  "Condition",
  "Batch",
  "Cell type"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_raw)
)

if (length(missing_columns) > 0) {

  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )

}

meta <- meta_raw %>%
  transmute(
    Sample = trimws(
      as.character(`FASTQ name`)
    ),
    Donor = trimws(
      as.character(`Donor ID`)
    ),
    Age = as.numeric(Age),
    Diagnosis = toupper(
      trimws(as.character(Condition))
    ),
    Batch = trimws(
      as.character(Batch)
    ),
    CellType = trimws(
      as.character(`Cell type`)
    )
  )

if (anyDuplicated(meta$Sample)) {
  stop("Duplicated sample names in metadata")
}

if (!setequal(meta$Sample, colnames(counts))) {

  print(
    setdiff(
      colnames(counts),
      meta$Sample
    )
  )

  print(
    setdiff(
      meta$Sample,
      colnames(counts)
    )
  )

  stop("Sample names do not match")
}

meta <- meta[
  match(
    colnames(counts),
    meta$Sample
  ),
  ,
  drop = FALSE
]

if (
  anyNA(meta$Sample) ||
  anyNA(meta$Donor) ||
  anyNA(meta$Diagnosis) ||
  anyNA(meta$CellType)
) {
  stop("Missing essential metadata")
}

if (!all(meta$Diagnosis %in% c("HC", "NPC"))) {
  stop("Diagnosis must be HC or NPC")
}

if (!all(meta$CellType %in% cell_types)) {
  stop("Unexpected cell type")
}

if (
  any(
    meta %>%
      group_by(Donor) %>%
      summarise(
        n_groups = n_distinct(Diagnosis),
        .groups = "drop"
      ) %>%
      pull(n_groups) > 1
  )
) {
  stop("One donor is assigned to multiple diagnoses")
}

meta$Diagnosis <- factor(
  meta$Diagnosis,
  levels = c("HC", "NPC")
)

meta$Donor <- factor(
  meta$Donor
)

write.csv(
  meta,
  file.path(
    outdir,
    "Included_Samples.csv"
  ),
  row.names = FALSE
)

donor_structure <- meta %>%
  group_by(
    CellType,
    Diagnosis,
    Donor
  ) %>%
  summarise(
    Samples = n(),
    Batches = paste(
      Batch,
      collapse = "; "
    ),
    .groups = "drop"
  )

write.csv(
  donor_structure,
  file.path(
    outdir,
    "Donor_Batch_Structure.csv"
  ),
  row.names = FALSE
)

log_step("Metadata completed")


# ============================================================
# 5. GENE ANNOTATION
# ============================================================

log_step("Gene annotation")

genes <- rownames(counts)

clean_genes <- sub(
  "\\.[0-9]+$",
  "",
  genes
)

id_type <- if (
  mean(grepl("^ENSG[0-9]+$", clean_genes)) > 0.5
) {
  "ENSEMBL"
} else {
  "SYMBOL"
}

unique_ids <- unique(clean_genes)

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique_ids,
  column = "SYMBOL",
  keytype = id_type,
  multiVals = "first"
)

entrez_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique_ids,
  column = "ENTREZID",
  keytype = id_type,
  multiVals = "first"
)

annotation_df <- data.frame(
  Gene = genes,
  SYMBOL = unname(
    symbol_map[clean_genes]
  ),
  ENTREZID = unname(
    entrez_map[clean_genes]
  ),
  stringsAsFactors = FALSE
)

annotation_df$SYMBOL[
  is.na(annotation_df$SYMBOL)
] <- annotation_df$Gene[
  is.na(annotation_df$SYMBOL)
]

write.csv(
  annotation_df,
  file.path(
    outdir,
    "Gene_Annotation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. FIGURE THEME
# ============================================================

theme_pub <- function() {

  theme_classic(
    base_size = 12
  ) +
    theme(
      text = element_text(
        family = "sans",
        color = "black"
      ),
      plot.title = element_text(
        size = 13,
        face = "bold"
      ),
      axis.title = element_text(
        size = 12
      ),
      axis.text = element_text(
        size = 10,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.55
      ),
      panel.grid = element_blank()
    )

}


# ============================================================
# 7. VOLCANO PLOT
# ============================================================

make_volcano <- function(
    df,
    cell_type,
    filename
) {

  plot_df <- df %>%
    filter(
      is.finite(logFC),
      !is.na(adj.P.Val)
    ) %>%
    mutate(
      Score = -log10(
        pmax(
          adj.P.Val,
          1e-300
        )
      )
    )

  labels <- plot_df %>%
    filter(
      Status != "NS"
    ) %>%
    group_by(Status) %>%
    arrange(
      adj.P.Val,
      .by_group = TRUE
    ) %>%
    slice_head(n = 8) %>%
    ungroup()

  p <- ggplot(
    plot_df,
    aes(
      x = logFC,
      y = Score
    )
  ) +

    geom_point(
      data = plot_df %>%
        filter(Status == "NS"),
      color = color_ns,
      size = 1,
      alpha = 0.5
    ) +

    geom_point(
      data = plot_df %>%
        filter(Status != "NS"),
      aes(color = Status),
      size = 1.2,
      alpha = 0.85
    ) +

    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      color = "grey50"
    ) +

    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey50"
    ) +

    ggrepel::geom_text_repel(
      data = labels,
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3,
      seed = 123
    ) +

    scale_color_manual(
      values = c(
        UP = color_up,
        DOWN = color_down
      )
    ) +

    labs(
      title = paste(
        cell_type,
        "| NPC vs HC"
      ),
      x = expression(log[2]~"fold change"),
      y = expression(-log[10]~"(adjusted P)")
    ) +

    theme_pub() +

    theme(
      legend.position = "none"
    )

  ggsave(
    filename,
    p,
    width = 7.5,
    height = 6,
    dpi = 600,
    bg = "white"
  )

}


# ============================================================
# 8. GO ANALYSIS
# ============================================================

run_go <- function(ids, universe) {

  ids <- unique(
    na.omit(as.character(ids))
  )

  universe <- unique(
    na.omit(as.character(universe))
  )

  ids <- intersect(
    ids,
    universe
  )

  if (length(ids) < 5) {
    return(data.frame(p.adjust = numeric(0)))
  }

  ego <- clusterProfiler::enrichGO(
    gene = ids,
    universe = universe,
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

  if (is.null(ego)) return(data.frame(p.adjust = numeric(0)))
  as.data.frame(ego)

}


# ============================================================
# 9. GO BAR PLOT
# ============================================================

make_go_bar <- function(
    df,
    direction,
    cell_type,
    filename
) {

  if (nrow(df) == 0) {
    return(NULL)
  }

  df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Term = stringr::str_wrap(
        Description,
        width = 43
      ),
      Score = -log10(
        pmax(p.adjust, 1e-300)
      )
    )

  df$TermID <- factor(
    df$ID,
    levels = rev(df$ID)
  )

  term_labels <- setNames(
    df$Term,
    df$ID
  )

  fill_color <- if (
    direction == "UP"
  ) {
    color_up
  } else {
    color_down
  }

  p <- ggplot(
    df,
    aes(
      x = Score,
      y = TermID
    )
  ) +

    geom_col(
      fill = fill_color,
      width = 0.72
    ) +

    scale_y_discrete(
      labels = term_labels
    ) +

    scale_x_continuous(
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +

    labs(
      title = paste(
        cell_type,
        "|",
        direction
      ),
      x = expression(-log[10]~"(adjusted P)"),
      y = NULL
    ) +

    theme_pub() +

    theme(
      axis.text.y = element_text(size = 9),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  ggsave(
    filename,
    p,
    width = 9,
    height = 5.5,
    dpi = 600,
    bg = "white"
  )

}


# ============================================================
# 10. GO DOT PLOT
# ============================================================

make_go_dot <- function(
    df,
    direction,
    cell_type,
    filename
) {

  if (nrow(df) == 0) {
    return(NULL)
  }

  df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Term = stringr::str_wrap(
        Description,
        width = 43
      )
    )

  df$GeneRatioNumeric <- vapply(
    strsplit(df$GeneRatio, "/"),
    function(x) {
      as.numeric(x[1]) /
        as.numeric(x[2])
    },
    numeric(1)
  )

  df$TermID <- factor(
    df$ID,
    levels = rev(df$ID)
  )

  term_labels <- setNames(
    df$Term,
    df$ID
  )

  p <- ggplot(
    df,
    aes(
      x = GeneRatioNumeric,
      y = TermID
    )
  ) +

    geom_point(
      aes(
        size = Count,
        color = p.adjust
      ),
      alpha = 0.95
    ) +

    scale_y_discrete(
      labels = term_labels
    ) +

    scale_size_continuous(
      range = c(3.5, 9),
      name = "Gene count"
    ) +

    scale_color_gradient(
      low = "#CC6859",
      high = "#3C78A6",
      limits = c(0, 0.05),
      oob = scales::squish,
      name = "Adjusted P"
    ) +

    labs(
      title = paste(
        cell_type,
        "|",
        direction
      ),
      x = "Gene ratio",
      y = NULL
    ) +

    theme_pub() +

    theme(
      axis.text.y = element_text(size = 9),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  ggsave(
    filename,
    p,
    width = 9,
    height = 5.5,
    dpi = 600,
    bg = "white"
  )

}


# ============================================================
# 11. MAIN ANALYSIS FUNCTION
# ============================================================

analyze_cell_type <- function(cell_type) {

  log_step(
    paste("START:", cell_type)
  )

  cell_dir <- file.path(
    outdir,
    gsub("-", "", cell_type)
  )

  dir.create(
    cell_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ----------------------------------------------------------
  # Samples
  # ----------------------------------------------------------

  meta_sub <- meta %>%
    filter(
      CellType == cell_type
    )

  meta_sub$Diagnosis <- droplevels(
    meta_sub$Diagnosis
  )

  meta_sub$Donor <- droplevels(
    meta_sub$Donor
  )

  counts_sub <- counts[
    ,
    meta_sub$Sample,
    drop = FALSE
  ]

  if (nlevels(meta_sub$Diagnosis) != 2) {
    stop("Both HC and NPC are required")
  }

  donor_n <- meta_sub %>%
    group_by(Diagnosis) %>%
    summarise(
      N = n_distinct(Donor),
      .groups = "drop"
    )

  if (any(donor_n$N < 2)) {
    stop("Too few independent donors")
  }

  # ----------------------------------------------------------
  # Filter
  # ----------------------------------------------------------

  keep <- rowSums(
    counts_sub >= min_count
  ) >= min_samples

  counts_sub <- counts_sub[
    keep,
    ,
    drop = FALSE
  ]

  log_step(
    paste(
      cell_type,
      "genes:",
      nrow(counts_sub)
    )
  )

  # ----------------------------------------------------------
  # TMM
  # ----------------------------------------------------------

  y <- edgeR::DGEList(
    counts = counts_sub
  )

  y <- edgeR::calcNormFactors(
    y,
    method = "TMM"
  )

  meta_model <- as.data.frame(
    meta_sub
  )

  rownames(meta_model) <- meta_model$Sample

  meta_model$Diagnosis <- factor(
    meta_model$Diagnosis,
    levels = c("HC", "NPC")
  )

  meta_model$Donor <- factor(
    meta_model$Donor
  )

  design <- model.matrix(
    ~ Diagnosis,
    data = meta_model
  )

  # ----------------------------------------------------------
  # Standard voom
  # ----------------------------------------------------------

  log_step(
    paste(
      cell_type,
      "voom START"
    )
  )

  v <- limma::voom(
    y,
    design = design,
    plot = FALSE,
    save.plot = FALSE
  )

  log_step(
    paste(
      cell_type,
      "voom FINISHED"
    )
  )

  # ----------------------------------------------------------
  # Donor correlation
  # ----------------------------------------------------------

  has_repeats <- any(
    table(meta_model$Donor) > 1
  )

  if (has_repeats) {

    log_step(
      paste(
        cell_type,
        "duplicateCorrelation START"
      )
    )

    corfit <- limma::duplicateCorrelation(
      v,
      design = design,
      block = meta_model$Donor
    )

    consensus_cor <- corfit$consensus.correlation

    if (!is.finite(consensus_cor)) {
      stop("Invalid donor correlation estimate")
    }

    log_step(
      paste(
        cell_type,
        "Donor correlation:",
        round(consensus_cor, 4)
      )
    )

    # Re-estimate voom weights using donor correlation

    v <- limma::voom(
      y,
      design = design,
      block = meta_model$Donor,
      correlation = consensus_cor,
      plot = FALSE
    )

    # Re-estimate correlation

    corfit <- limma::duplicateCorrelation(
      v,
      design = design,
      block = meta_model$Donor
    )

    consensus_cor <- corfit$consensus.correlation

    log_step(
      paste(
        cell_type,
        "duplicateCorrelation FINISHED"
      )
    )

    fit <- limma::lmFit(
      v,
      design = design,
      block = meta_model$Donor,
      correlation = consensus_cor
    )

    model_name <- "voom_duplicateCorrelation"

  } else {

    fit <- limma::lmFit(
      v,
      design = design
    )

    consensus_cor <- NA_real_

    model_name <- "voom_standard"

  }

  # ----------------------------------------------------------
  # Empirical Bayes
  # ----------------------------------------------------------

  log_step(
    paste(
      cell_type,
      "eBayes START"
    )
  )

  fit <- limma::eBayes(
    fit,
    robust = FALSE
  )

  log_step(
    paste(
      cell_type,
      "eBayes FINISHED"
    )
  )

  # ----------------------------------------------------------
  # DEG
  # ----------------------------------------------------------

  res <- limma::topTable(
    fit,
    coef = "DiagnosisNPC",
    number = Inf,
    sort.by = "P"
  )

  res$Gene <- rownames(res)

  res_df <- res %>%
    left_join(
      annotation_df,
      by = "Gene"
    )

  res_df$Status <- "NS"

  res_df$Status[
    !is.na(res_df$adj.P.Val) &
      res_df$adj.P.Val < padj_cutoff &
      res_df$logFC > log2fc_cutoff
  ] <- "UP"

  res_df$Status[
    !is.na(res_df$adj.P.Val) &
      res_df$adj.P.Val < padj_cutoff &
      res_df$logFC < -log2fc_cutoff
  ] <- "DOWN"

  deg_up <- res_df %>%
    filter(Status == "UP")

  deg_down <- res_df %>%
    filter(Status == "DOWN")

  deg_all <- res_df %>%
    filter(Status != "NS")

  write.csv(
    res_df,
    file.path(
      cell_dir,
      "All_Genes_NPC_vs_HC.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    deg_up,
    file.path(
      cell_dir,
      "DEG_UP.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    deg_down,
    file.path(
      cell_dir,
      "DEG_DOWN.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    deg_all,
    file.path(
      cell_dir,
      "Significant_DEGs.csv"
    ),
    row.names = FALSE
  )

  log_step(
    paste(
      cell_type,
      "DEG UP:",
      nrow(deg_up),
      "DOWN:",
      nrow(deg_down)
    )
  )

  # ----------------------------------------------------------
  # Volcano
  # ----------------------------------------------------------

  make_volcano(
    res_df,
    cell_type,
    file.path(
      cell_dir,
      "01_Volcano.png"
    )
  )

  log_step(
    paste(
      cell_type,
      "Volcano saved"
    )
  )

  # ----------------------------------------------------------
  # GO
  # ----------------------------------------------------------

  universe <- res_df %>%
    filter(
      !is.na(adj.P.Val),
      !is.na(ENTREZID)
    ) %>%
    pull(ENTREZID) %>%
    unique()

  go_up <- run_go(
    deg_up$ENTREZID,
    universe
  )

  go_down <- run_go(
    deg_down$ENTREZID,
    universe
  )

  go_up_sig <- go_up %>%
    filter(
      !is.na(p.adjust),
      p.adjust < go_padj_cutoff
    )

  go_down_sig <- go_down %>%
    filter(
      !is.na(p.adjust),
      p.adjust < go_padj_cutoff
    )

  write.csv(
    go_up,
    file.path(
      cell_dir,
      "GO_BP_UP_All.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    go_down,
    file.path(
      cell_dir,
      "GO_BP_DOWN_All.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    go_up_sig,
    file.path(
      cell_dir,
      "GO_BP_UP_Significant.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    go_down_sig,
    file.path(
      cell_dir,
      "GO_BP_DOWN_Significant.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    head(go_up_sig, top_n),
    file.path(
      cell_dir,
      "GO_BP_UP_Top10.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    head(go_down_sig, top_n),
    file.path(
      cell_dir,
      "GO_BP_DOWN_Top10.csv"
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # GO figures
  # ----------------------------------------------------------

  make_go_bar(
    go_up_sig,
    "UP",
    cell_type,
    file.path(
      cell_dir,
      "02_GO_BP_UP_Bar.png"
    )
  )

  make_go_bar(
    go_down_sig,
    "DOWN",
    cell_type,
    file.path(
      cell_dir,
      "03_GO_BP_DOWN_Bar.png"
    )
  )

  make_go_dot(
    go_up_sig,
    "UP",
    cell_type,
    file.path(
      cell_dir,
      "04_GO_BP_UP_Dot.png"
    )
  )

  make_go_dot(
    go_down_sig,
    "DOWN",
    cell_type,
    file.path(
      cell_dir,
      "05_GO_BP_DOWN_Dot.png"
    )
  )

  log_step(
    paste(
      cell_type,
      "GO analysis completed"
    )
  )

  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------

  summary_df <- data.frame(
    CellType = cell_type,
    Model = model_name,
    Samples = nrow(meta_sub),
    Donors = n_distinct(meta_sub$Donor),
    HC_Donors = n_distinct(
      meta_sub$Donor[
        meta_sub$Diagnosis == "HC"
      ]
    ),
    NPC_Donors = n_distinct(
      meta_sub$Donor[
        meta_sub$Diagnosis == "NPC"
      ]
    ),
    GenesTested = nrow(res_df),
    DonorCorrelation = consensus_cor,
    DEG_UP = nrow(deg_up),
    DEG_DOWN = nrow(deg_down),
    GO_UP = nrow(go_up_sig),
    GO_DOWN = nrow(go_down_sig)
  )

  write.csv(
    summary_df,
    file.path(
      cell_dir,
      "Analysis_Summary.csv"
    ),
    row.names = FALSE
  )

  log_step(
    paste(
      cell_type,
      "COMPLETED"
    )
  )

  summary_df

}


# ============================================================
# 12. RUN BOTH CELL TYPES
# ============================================================

summary_list <- lapply(
  cell_types,
  analyze_cell_type
)

summary_all <- dplyr::bind_rows(
  summary_list
)

write.csv(
  summary_all,
  file.path(
    outdir,
    "Summary_UDCiN_iPSCiN.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. ANALYSIS SETTINGS
# ============================================================

settings <- c(
  "RNA-seq NPC vs HC",
  "UDC-iN and iPSC-iN analyzed separately",
  "All differentiation batches included",
  "TMM normalization",
  "Standard limma voom",
  "duplicateCorrelation for repeated donors",
  "Common within-donor correlation",
  "No dream",
  "No voomWithDreamWeights",
  "No parallel processing",
  "Count >=10 in >=3 samples",
  "DEG BH adjusted P <0.05",
  "DEG absolute log2FC >1",
  "GO Biological Process ORA",
  "GO UP and DOWN separately",
  "GO background: mapped tested genes",
  "GO adjusted P <0.05",
  "GO top 10",
  "Age not adjusted",
  "PNG 600 dpi"
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


# ============================================================
# 14. FINAL SUMMARY
# ============================================================

log_step("ALL ANALYSES COMPLETED")

print(summary_all)

cat("\n====================================\n")
cat("FULL ANALYSIS COMPLETED\n")
cat("====================================\n")

cat("\nOutput directory:\n")
cat(outdir, "\n")

if (Sys.info()[["sysname"]] == "Darwin") {

  system2(
    "open",
    args = shQuote(outdir),
    wait = FALSE
  )

}
