# ============================================================
# HC Proteome Age GSEA
#
# UDC_Age / UDCiN_Age / Interaction
#
# SAME METHOD AS ORIGINAL INTERACTION GSEA
#
# - Existing limma moderated t-statistics
# - No model refitting
# - No protein-level significance prefilter
# - SYMBOL -> ENTREZID via bitr()
# - Mean moderated t across protein rows per ENTREZID
# - GO Biological Process + Reactome
# - fgsea multilevel
# - nPermSimple = 10000
# - eps = 0
# - minGSSize = 10, maxGSSize = 500
# - BH FDR < 0.05
#
# Interaction = UDC-iN age slope - UDC age slope
# ============================================================
options(
  stringsAsFactors = FALSE,
  repos = c(CRAN = "https://cloud.r-project.org")
)
# ============================================================
# 1. Packages
# ============================================================
cran <- c(
  "readxl",
  "readr",
  "dplyr",
  "ggplot2",
  "stringr",
  "scales"
)
bioc <- c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
)
for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(
      p,
      ask = FALSE,
      update = FALSE
    )
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
})
set.seed(20261016)
# ============================================================
# 2. Paths
# ============================================================
D <- path.expand(Sys.getenv("PROTEOMICS_DATA_DIR", unset = "~/Desktop"))
root <- file.path(
  D,
  "Proteome_HC_Age_JointModel_Sex_WithBatch"
)
if (!dir.exists(root)) {
  stop("Proteome result folder not found: ", root)
}
contrasts <- c(
  "UDC_Age",
  "UDCiN_Age",
  "Interaction"
)
# New common output folder
# Original Interaction GSEA remains untouched.
out_root <- file.path(
  root,
  "GSEA_All_Contrasts_Same_Method_nPerm10000"
)
dir.create(
  out_root,
  recursive = TRUE,
  showWarnings = FALSE
)
# ============================================================
# 3. Settings
# ============================================================
FDR_CUTOFF <- 0.05
MIN_GS <- 10
MAX_GS <- 500
N_PERM <- 10000
SEED <- 20261016
TOP_N <- 10
# ============================================================
# 4. Output helpers
# ============================================================
W <- function(x, name, folder) {
  readr::write_csv(
    as.data.frame(x),
    file.path(folder, name),
    na = ""
  )
}
save_plot <- function(
    p,
    name,
    folder,
    height = 6,
    width = 9.5
) {
  ggplot2::ggsave(
    filename = file.path(
      folder,
      paste0(name, ".png")
    ),
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = file.path(
      folder,
      paste0(name, ".pdf")
    ),
    plot = p,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}
# ============================================================
# 5. Find original proteome workbook
#
# Only needed when SYMBOL is absent from All_Proteins.csv
# ============================================================
find_workbook <- function() {
  exact <- file.path(
    D,
    "20230316_iN_batch_All_ys_Proteins_KM edited.xlsx"
  )
  if (file.exists(exact)) {
    return(exact)
  }
  hits <- list.files(
    D,
    pattern = "^20230316_iN_batch_All_ys_Proteins_KM edited.*\\.xlsx$",
    full.names = TRUE
  )
  if (length(hits) != 1) {
    stop(
      "Cannot uniquely identify proteome workbook.\n",
      "Matches: ",
      paste(basename(hits), collapse = ", ")
    )
  }
  hits
}
# ============================================================
# 6. Load original workbook annotation
#
# ProteinID = original worksheet row number
# BEFORE nonempty-Accession filtering.
#
# Never renumber rows after filtering.
# ============================================================
load_workbook_annotation <- function() {
  workbook <- find_workbook()
  sheets <- readxl::excel_sheets(workbook)
  sheet <- setdiff(sheets, "metadata")[1]
  if (is.na(sheet)) {
    stop("No usable worksheet found.")
  }
  raw <- as.data.frame(
    readxl::read_excel(
      workbook,
      sheet = sheet,
      .name_repair = "minimal"
    )
  )
  required <- c(
    "Accession",
    "Gene Symbol"
  )
  if (!all(required %in% names(raw))) {
    stop(
      "Workbook missing columns: ",
      paste(
        setdiff(required, names(raw)),
        collapse = ", "
      )
    )
  }
  data.frame(
    ProteinID = as.character(seq_len(nrow(raw))),
    GeneSymbol = trimws(
      as.character(raw[["Gene Symbol"]])
    ),
    WorkbookAccession = trimws(
      as.character(raw[["Accession"]])
    ),
    stringsAsFactors = FALSE
  )
}
# ============================================================
# 7. Obtain Gene Symbols
#
# Same logic as original Interaction script:
#
# 1. Prefer SYMBOL already in All_Proteins.csv
# 2. Otherwise match original worksheet rows
# ============================================================
prepare_proteins <- function(
    input,
    output
) {
  if (!file.exists(input)) {
    stop("Missing input: ", input)
  }
  res <- as.data.frame(
    readr::read_csv(
      input,
      show_col_types = FALSE
    )
  )
  required <- c("ProteinID", "t")
  if (!all(required %in% names(res))) {
    stop(
      "Missing columns in ",
      input,
      ": ",
      paste(
        setdiff(required, names(res)),
        collapse = ", "
      )
    )
  }
  res$ProteinID <- as.character(res$ProteinID)
  res$t <- suppressWarnings(
    as.numeric(res$t)
  )
  if (anyDuplicated(res$ProteinID)) {
    stop(
      "Duplicated ProteinID: ",
      input
    )
  }
  if ("SYMBOL" %in% names(res)) {
    res$GeneSymbol <- trimws(
      as.character(res$SYMBOL)
    )
    qc <- data.frame(
      Source = "Original CSV SYMBOL",
      OriginalRows = nrow(res),
      WithGeneSymbol = sum(
        !is.na(res$GeneSymbol) &
          nzchar(res$GeneSymbol)
      ),
      WithFiniteT = sum(
        is.finite(res$t)
      )
    )
  } else {
    anno <- load_workbook_annotation()
    res <- dplyr::left_join(
      res,
      anno,
      by = "ProteinID"
    )
    matched <- !is.na(res$WorkbookAccession) &
      nzchar(res$WorkbookAccession)
    accession_consistent <- rep(
      NA,
      nrow(res)
    )
    if ("Accession" %in% names(res)) {
      accession_consistent <- matched &
        trimws(
          as.character(res$Accession)
        ) == res$WorkbookAccession
    }
    qc <- data.frame(
      Source = "Original worksheet row ID",
      OriginalRows = nrow(res),
      MatchedWorkbookRows = sum(matched),
      MatchingAccession = sum(
        accession_consistent,
        na.rm = TRUE
      ),
      WithFiniteT = sum(
        is.finite(res$t)
      )
    )
    W(
      qc,
      "00_Input_Match_QC.csv",
      output
    )
    if (sum(matched) < 0.9 * nrow(res)) {
      stop(
        "ProteinID-to-workbook matching failed.\n",
        "Inspect 00_Input_Match_QC.csv"
      )
    }
    if (
      "Accession" %in% names(res) &&
      any(
        matched & !accession_consistent,
        na.rm = TRUE
      )
    ) {
      stop(
        "Accession mismatch between CSV and workbook.\n",
        "Check workbook version."
      )
    }
  }
  W(
    qc,
    "00_Input_Match_QC.csv",
    output
  )
  W(
    res,
    "01_Proteins_With_GeneSymbols.csv",
    output
  )
  res
}
# ============================================================
# 8. SYMBOL -> ENTREZID
#
# EXACT same principle as original Interaction code
#
# - Exclude ambiguous Gene Symbols
# - bitr mapping
# - Mean moderated t per ENTREZID
# ============================================================
build_ranking <- function(
    res,
    output
) {
  gs <- res |>
    dplyr::filter(
      !is.na(GeneSymbol),
      nzchar(GeneSymbol),
      !grepl("[,;/|]", GeneSymbol),
      is.finite(t)
    )
  valid_symbols <- AnnotationDbi::keys(
    org.Hs.eg.db::org.Hs.eg.db,
    keytype = "SYMBOL"
  )
  valid <- intersect(
    unique(gs$GeneSymbol),
    valid_symbols
  )
  qc <- data.frame(
    Metric = c(
      "FiniteT_and_single_SYMBOL_entries",
      "Distinct_input_SYMBOLs",
      "Valid_orgHs_SYMBOLs"
    ),
    Value = c(
      nrow(gs),
      length(unique(gs$GeneSymbol)),
      length(valid)
    )
  )
  W(
    qc,
    "02_SYMBOL_Mapping_QC.csv",
    output
  )
  if (length(valid) < 100) {
    stop(
      "Only ", length(valid),
      " valid Gene Symbols.\n",
      "Inspect 02_SYMBOL_Mapping_QC.csv"
    )
  }
  mapped <- suppressMessages(
    clusterProfiler::bitr(
      valid,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  )
  if (
    is.null(mapped) ||
    nrow(mapped) < 100
  ) {
    stop(
      "SYMBOL-to-ENTREZID mapping failed."
    )
  }
  # Remove existing ENTREZID before join.
  # This prevents ENTREZID.x / ENTREZID.y.
  if ("ENTREZID" %in% names(gs)) {
    gs <- dplyr::select(
      gs,
      -dplyr::all_of("ENTREZID")
    )
  }
  mapped <- mapped |>
    dplyr::transmute(
      SYMBOL = as.character(SYMBOL),
      ENTREZID = as.character(ENTREZID)
    ) |>
    dplyr::distinct(
      SYMBOL,
      ENTREZID
    )
  mapped_proteins <- dplyr::inner_join(
    gs,
    mapped,
    by = c("GeneSymbol" = "SYMBOL"),
    relationship = "many-to-many"
  ) |>
    dplyr::distinct(
      ENTREZID,
      ProteinID,
      .keep_all = TRUE
    )
  if (!"ENTREZID" %in% names(mapped_proteins)) {
    stop("ENTREZID missing after mapping.")
  }
  ranked <- mapped_proteins |>
    dplyr::group_by(ENTREZID) |>
    dplyr::summarise(
      RankScore = mean(t),
      ProteinRows = dplyr::n(),
      GeneSymbols = paste(
        sort(unique(GeneSymbol)),
        collapse = ";"
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(RankScore),
      ENTREZID
    )
  ranks <- ranked$RankScore
  names(ranks) <- ranked$ENTREZID
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  if (
    length(ranks) < 100 ||
    anyDuplicated(names(ranks))
  ) {
    stop("Invalid ranked gene vector.")
  }
  W(
    mapped_proteins,
    "03_Mapped_Protein_Entries.csv",
    output
  )
  W(
    ranked,
    "04_Ranked_Genes_Mean_ModeratedT.csv",
    output
  )
  list(
    ranks = ranks,
    ranked = ranked
  )
}
# ============================================================
# 9. Plot GSEA results
# ============================================================
plot_gsea <- function(
    z,
    contrast,
    db,
    direction,
    output
) {
  if (!nrow(z)) {
    return(invisible(NULL))
  }
  z$PlotNES <- if (
    direction == "UP"
  ) {
    z$NES
  } else {
    -z$NES
  }
  z$Term <- stringr::str_wrap(
    z$Description,
    width = 43
  )
  z <- z[
    order(
      -z$PlotNES,
      z$p.adjust,
      z$ID
    ),
    ,
    drop = FALSE
  ]
  z$TermID <- factor(
    z$ID,
    levels = rev(z$ID)
  )
  lab <- stats::setNames(
    z$Term,
    z$ID
  )
  cols <- if (direction == "UP") {
    c("#B9473D", "#F0C6B8")
  } else {
    c("#235C91", "#BFD8ED")
  }
  p <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = PlotNES,
      y = TermID
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = setSize,
        color = p.adjust
      ),
      alpha = 0.95
    ) +
    ggplot2::scale_y_discrete(
      labels = lab
    ) +
    ggplot2::scale_size_continuous(
      range = c(3.5, 9),
      name = "Gene-set size"
    ) +
    ggplot2::scale_color_gradient(
      low = cols[1],
      high = cols[2],
      limits = c(0, 0.05),
      oob = scales::squish,
      name = "BH FDR"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(0.04, 0.12)
      )
    ) +
    ggplot2::labs(
      title = paste0(
        contrast,
        " | ",
        db,
        " | ",
        direction
      ),
      x = if (
        direction == "UP"
      ) {
        "NES"
      } else {
        "-NES"
      },
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 12
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14
      ),
      axis.text.y = ggplot2::element_text(
        size = 9
      ),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  save_plot(
    p,
    paste0(
      "08_",
      db,
      "_",
      direction,
      "_NES_Dot"
    ),
    output,
    height = max(
      5.5,
      2.7 + 0.36 * nrow(z)
    )
  )
}
# ============================================================
# 10. Process GSEA results
# ============================================================
process_gsea <- function(
    obj,
    contrast,
    db,
    output,
    ranked_genes
) {
  tab <- if (is.null(obj)) {
    data.frame()
  } else {
    as.data.frame(obj)
  }
  if (!nrow(tab)) {
    W(
      tab,
      paste0("05_", db, "_ALL.csv"),
      output
    )
    W(
      tab,
      paste0("06_", db, "_FDR005.csv"),
      output
    )
    return(
      data.frame(
        Contrast = contrast,
        Database = db,
        RankedGenes = ranked_genes,
        TestedTerms = 0,
        Significant = 0,
        UP = 0,
        DOWN = 0
      )
    )
  }
  tab$Direction <- ifelse(
    tab$NES > 0,
    "UP",
    "DOWN"
  )
  W(
    tab,
    paste0("05_", db, "_ALL.csv"),
    output
  )
  sig <- tab[
    is.finite(tab$p.adjust) &
      tab$p.adjust < FDR_CUTOFF,
    ,
    drop = FALSE
  ]
  W(
    sig,
    paste0("06_", db, "_FDR005.csv"),
    output
  )
  summary <- data.frame(
    Contrast = contrast,
    Database = db,
    RankedGenes = ranked_genes,
    TestedTerms = nrow(tab),
    Significant = nrow(sig),
    UP = sum(sig$NES > 0),
    DOWN = sum(sig$NES < 0)
  )
  for (direction in c("UP", "DOWN")) {
    z <- sig[
      sig$Direction == direction,
      ,
      drop = FALSE
    ]
    z <- z[
      order(
        z$p.adjust,
        -abs(z$NES),
        z$ID
      ),
      ,
      drop = FALSE
    ]
    z <- head(
      z,
      TOP_N
    )
    W(
      z,
      paste0(
        "07_",
        db,
        "_",
        direction,
        "_Top10.csv"
      ),
      output
    )
    if (!nrow(z)) next
    plot_gsea(
      z = z,
      contrast = contrast,
      db = db,
      direction = direction,
      output = output
    )
  }
  summary
}
# ============================================================
# 11. Main GSEA loop
# ============================================================
all_summaries <- list()
for (contrast in contrasts) {
  message("\n")
  message("====================================")
  message("GSEA: ", contrast)
  message("====================================")
  input <- file.path(
    root,
    contrast,
    "All_Proteins.csv"
  )
  output <- file.path(
    out_root,
    contrast
  )
  dir.create(
    output,
    recursive = TRUE,
    showWarnings = FALSE
  )
  # ----------------------------------------------------------
  # Original limma results
  # ----------------------------------------------------------
  res <- prepare_proteins(
    input = input,
    output = output
  )
  # ----------------------------------------------------------
  # SYMBOL -> ENTREZID -> mean moderated t
  # ----------------------------------------------------------
  ranking <- build_ranking(
    res = res,
    output = output
  )
  ranks <- ranking$ranks
  message(
    "Ranked genes: ",
    length(ranks)
  )
  # ----------------------------------------------------------
  # GO BP GSEA
  # ----------------------------------------------------------
  set.seed(SEED)
  go <- clusterProfiler::gseGO(
    geneList = ranks,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = MIN_GS,
    maxGSSize = MAX_GS,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    by = "fgsea",
    seed = TRUE,
    nPermSimple = N_PERM
  )
  # ----------------------------------------------------------
  # Reactome GSEA
  # ----------------------------------------------------------
  set.seed(SEED)
  reactome <- ReactomePA::gsePathway(
    geneList = ranks,
    organism = "human",
    minGSSize = MIN_GS,
    maxGSSize = MAX_GS,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    by = "fgsea",
    seed = TRUE,
    nPermSimple = N_PERM
  )
  # ----------------------------------------------------------
  # Save results
  # ----------------------------------------------------------
  summary_go <- process_gsea(
    obj = go,
    contrast = contrast,
    db = "GO_BP",
    output = output,
    ranked_genes = length(ranks)
  )
  summary_reactome <- process_gsea(
    obj = reactome,
    contrast = contrast,
    db = "Reactome",
    output = output,
    ranked_genes = length(ranks)
  )
  summary_df <- dplyr::bind_rows(
    summary_go,
    summary_reactome
  )
  W(
    summary_df,
    "09_GSEA_Summary.csv",
    output
  )
  all_summaries[[contrast]] <- summary_df
  # ----------------------------------------------------------
  # Translation elongation enrichment plot
  # Same as original Interaction code
  # ----------------------------------------------------------
  react_tab <- if (is.null(reactome)) {
    data.frame()
  } else {
    as.data.frame(reactome)
  }
  if (nrow(react_tab)) {
    elong <- react_tab[
      grepl(
        "Eukaryotic Translation Elongation",
        react_tab$Description,
        ignore.case = TRUE
      ),
      ,
      drop = FALSE
    ]
    W(
      elong,
      "10_Eukaryotic_Translation_Elongation.csv",
      output
    )
    if (nrow(elong)) {
      eg <- tryCatch(
        enrichplot::gseaplot2(
          reactome,
          geneSetID = elong$ID[1],
          title = paste0(
            elong$Description[1],
            " | NES=",
            sprintf("%.2f", elong$NES[1]),
            ", FDR=",
            format(
              elong$p.adjust[1],
              digits = 2,
              scientific = TRUE
            )
          )
        ),
        error = function(e) {
          warning(
            "Enrichment plot unavailable: ",
            conditionMessage(e)
          )
          NULL
        }
      )
      if (!is.null(eg)) {
        tryCatch(
          save_plot(
            eg,
            "10_Reactome_Translation_Elongation_Enrichment",
            output,
            height = 6,
            width = 9
          ),
          error = function(e) {
            warning(
              "Could not save enrichment plot: ",
              conditionMessage(e)
            )
          }
        )
      }
    }
  }
  # ----------------------------------------------------------
  # Analysis settings
  # ----------------------------------------------------------
  writeLines(
    c(
      paste("Contrast:", contrast),
      "",
      "Input: existing limma All_Proteins.csv.",
      "No model refitting.",
      "No protein significance prefilter.",
      "",
      "Ranking:",
      "Unambiguous SYMBOL -> ENTREZID via bitr.",
      "Mean original moderated t per ENTREZID.",
      "",
      "GO BP: gseGO.",
      "Reactome: gsePathway.",
      "Algorithm: fgsea.",
      "nPermSimple = 10000.",
      "eps = 0.",
      "minGSSize = 10.",
      "maxGSSize = 500.",
      "pvalueCutoff = 1.",
      "pAdjustMethod = BH.",
      "Significance: FDR < 0.05.",
      "",
      "BH adjustment is performed separately",
      "within each database and contrast.",
      "",
      "UDC_Age: age slope in UDC.",
      "UDCiN_Age: age slope in UDC-iN.",
      "Interaction: UDC-iN slope minus UDC slope."
    ),
    file.path(
      output,
      "Analysis_Settings.txt"
    )
  )
  capture.output(
    sessionInfo(),
    file = file.path(
      output,
      "SessionInfo.txt"
    )
  )
  print(summary_df)
}
# ============================================================
# 12. Combined summary
# ============================================================
combined <- dplyr::bind_rows(
  all_summaries
)
W(
  combined,
  "GSEA_All_Contrasts_Summary.csv",
  out_root
)
print(combined)
message("\n====================================")
message("ALL GSEA ANALYSES COMPLETED")
message("====================================")
message("Output: ", out_root)
