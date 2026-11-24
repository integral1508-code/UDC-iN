# =====================================================================
# UDC-iN vs UDC | HC paired donors | Publication-ready marker heatmap
# Ten curated markers,
# paired design validation, TMM log2CPM, UNCLIPPED row Z-score export.
# Paired differential-expression QC via limma-voom
# =====================================================================

# --- Packages --------------------------------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")
for (p in c("readxl","pheatmap"))
  if (!requireNamespace(p,quietly=TRUE))
    install.packages(p,repos="https://cloud.r-project.org")
for (p in c("edgeR","limma"))
  if (!requireNamespace(p,quietly=TRUE))
    BiocManager::install(p,ask=FALSE,update=FALSE)
library(readxl)
library(edgeR)
library(limma)
library(pheatmap)

# --- Files and plot options ------------------------------------------
input_dir <- Sys.getenv("UDC_INPUT_DIR", unset = "data")
output_dir <- Sys.getenv("UDC_OUTPUT_DIR", unset = "results")
count_file <- file.path(input_dir,"RNAseq_gene_count_matrix_22samples_GEO.csv")
meta_file  <- file.path(input_dir,"RNA-seq_Donor information.xlsx")
outdir <- file.path(output_dir,"UDCiN_vs_UDC_MainFigure_10Markers")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

# Fixed, biologically curated marker list. Do NOT re-select genes based
# on the resulting heatmap or differential-expression P values.
markers <- list(
  "Epithelial"    = c("EPCAM","KRT8"),
  "Somatic / ECM" = c("VIM","COL1A1"),
  "Proliferation" = c("MKI67"),
  "Neuronal"      = c("MAP2","MAPT","SYN1","SYP","SNAP25")
)
genes <- unlist(markers,use.names=FALSE)
stopifnot(length(genes)==10L,!anyDuplicated(genes))
color_limit <- 1.5   # display only: values outside +/-1.5 are saturated

# --- Metadata: exactly matched UDC / UDC-iN pairs -------------------
if (!file.exists(count_file)) stop("Missing: ",count_file)
if (!file.exists(meta_file)) stop("Missing: ",meta_file)
meta <- as.data.frame(readxl::read_excel(meta_file))
required <- c("FASTQ name","Donor ID","Age","Condition","Cell type")
if (!all(required %in% names(meta)))
  stop("Metadata missing: ",paste(setdiff(required,names(meta)),collapse=", "))

meta <- meta[!is.na(meta[["FASTQ name"]]) &
               meta[["Condition"]]=="HC" &
               meta[["Cell type"]] %in% c("UDC","UDC-iN"),,drop=FALSE]
meta$Sample <- as.character(meta[["FASTQ name"]])
meta$Donor  <- as.character(meta[["Donor ID"]])
meta$AgeNum <- suppressWarnings(as.numeric(meta[["Age"]]))
if (anyNA(meta$Sample) || anyNA(meta$Donor) ||
    anyNA(meta$AgeNum) || anyDuplicated(meta$Sample))
  stop("Invalid sample/donor/age metadata")

pairtab <- table(meta$Donor,meta[["Cell type"]])
if (!all(c("UDC","UDC-iN") %in% colnames(pairtab)))
  stop("UDC or UDC-iN absent from metadata")
if (any(pairtab[,"UDC"]!=1L | pairtab[,"UDC-iN"]!=1L))
  stop("Every included HC donor must have exactly 1 UDC and 1 UDC-iN")
donor_age <- tapply(meta$AgeNum,meta$Donor,function(v) {
  if (length(unique(v))!=1L) stop("Paired age mismatch")
  v[1]
})
donors <- names(donor_age)[order(as.numeric(donor_age),names(donor_age))]
get_samples <- function(ct) vapply(donors,function(d) {
  x <- meta$Sample[meta$Donor==d & meta[["Cell type"]]==ct]
  if (length(x)!=1L) stop("Unpaired donor: ",d)
  x
},character(1))
# vapply() attaches donor IDs as names; remove them explicitly.
# identical() compares attributes as well as values.
sample_order <- unname(c(get_samples("UDC"), get_samples("UDC-iN")))
if (anyDuplicated(sample_order)) stop("Duplicated sample IDs in paired order")
matched_idx <- match(sample_order, meta$Sample)
if (anyNA(matched_idx)) stop("Some paired samples could not be matched to metadata")
meta <- meta[matched_idx, , drop=FALSE]
rownames(meta) <- NULL
if (!identical(unname(as.character(meta$Sample)), sample_order))
  stop("Metadata sample order does not match the expected paired order")
# Verify the expected donor order on both sides of the heatmap.
if (!identical(unname(as.character(meta$Donor)), rep(donors, 2L)))
  stop("UDC and UDC-iN donor pair order mismatch")
if (!identical(unname(as.character(meta[["Cell type"]])),
               rep(c("UDC", "UDC-iN"), each=length(donors))))
  stop("Cell type order mismatch")
if (length(donors)!=11L)
  warning("Expected 11 donors; found ",length(donors))
labels <- paste0("HC. ",format(meta$AgeNum,trim=TRUE,scientific=FALSE),"y")
meta$HeatmapLabel <- labels
meta$PlotOrder <- seq_len(nrow(meta))

# --- Count matrix / TMM normalization --------------------------------
raw <- read.csv(count_file,check.names=FALSE,stringsAsFactors=FALSE)
if (!("Gene" %in% names(raw))) stop("Missing Gene column")
if (anyNA(raw$Gene) || anyDuplicated(raw$Gene))
  stop("Gene symbols missing or duplicated; resolve before analysis")
missing_samples <- setdiff(sample_order,names(raw))
if (length(missing_samples))
  stop("Count matrix missing samples: ",paste(missing_samples,collapse=", "))
counts <- as.matrix(raw[,sample_order,drop=FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- raw$Gene
if (anyNA(counts) || any(!is.finite(counts)) ||
    any(counts<0) || any(abs(counts-round(counts))>1e-6))
  stop("Expected finite, nonnegative integer raw counts")

missing_genes <- setdiff(genes,rownames(counts))
if (length(missing_genes))
  stop("Missing marker genes: ",paste(missing_genes,collapse=", "))

# TMM estimated using ALL genes; not just the selected markers.
y <- edgeR::DGEList(counts=counts)
y <- edgeR::calcNormFactors(y,method="TMM")
logcpm <- edgeR::cpm(y,log=TRUE,prior.count=2)
expr <- logcpm[genes,,drop=FALSE]

# Raw row-wise Z scores; NEVER clip values in the exported matrix.
z <- t(scale(t(expr)))
bad <- !is.finite(rowSums(z))
if (any(bad)) {
  warning("Constant-expression markers: ",paste(rownames(z)[bad],collapse=", "))
  z[bad,] <- 0
}

# Low expression is important to flag, even if Z-score colors look strong.
cpm_linear <- edgeR::cpm(y,log=FALSE)
qc <- data.frame(
  Gene=genes,
  Category=rep(names(markers),lengths(markers)),
  Mean_CPM_UDC=rowMeans(cpm_linear[genes,seq_along(donors),drop=FALSE]),
  Mean_CPM_UDCiN=rowMeans(cpm_linear[genes,length(donors)+seq_along(donors),drop=FALSE]),
  Mean_log2CPM_UDC=rowMeans(expr[,seq_along(donors),drop=FALSE]),
  Mean_log2CPM_UDCiN=rowMeans(expr[,length(donors)+seq_along(donors),drop=FALSE]),
  Mean_paired_log2CPM_change=rowMeans(
    expr[,length(donors)+seq_along(donors),drop=FALSE] -
    expr[,seq_along(donors),drop=FALSE]
  ),
  stringsAsFactors=FALSE
)
qc$Low_expression_both_groups <- qc$Mean_CPM_UDC<1 & qc$Mean_CPM_UDCiN<1
if (any(qc$Low_expression_both_groups))
  warning("Markers with mean CPM <1 in BOTH groups: ",
          paste(qc$Gene[qc$Low_expression_both_groups],collapse=", "))

write.csv(meta,file.path(outdir,"Sample_Metadata_Paired.csv"),row.names=FALSE)
write.csv(data.frame(Gene=genes,expr,check.names=FALSE),
          file.path(outdir,"Marker_TMM_log2CPM.csv"),row.names=FALSE)
write.csv(data.frame(Gene=genes,z,check.names=FALSE),
          file.path(outdir,"Marker_RowZscore_UNCLIPPED.csv"),row.names=FALSE)
write.csv(qc,file.path(outdir,"Marker_Expression_QC.csv"),row.names=FALSE)

# --- OPTIONAL statistical check: paired limma-voom across ALL genes --
# FDR is calculated across all tested genes, NOT just selected markers.
# These results are supplemental QC; the heatmap itself does not imply
# statistical significance, and marker choice is descriptive.
keep <- edgeR::filterByExpr(y,group=meta[["Cell type"]])
if (sum(keep)<100) stop("Too few genes passed expression filtering")
yf <- y[keep,,keep.lib.sizes=FALSE]
yf <- edgeR::calcNormFactors(yf,method="TMM")
cell_type <- factor(meta[["Cell type"]],levels=c("UDC","UDC-iN"))
donor <- factor(meta$Donor,levels=donors)
design <- model.matrix(~ donor + cell_type)
colnames(design) <- make.names(colnames(design))
coef_name <- grep("^cell_typeUDC.iN$",colnames(design),value=TRUE)
if (length(coef_name)!=1L) stop("Cannot find UDC-iN contrast")
v <- limma::voom(yf,design,plot=FALSE)
fit <- limma::eBayes(limma::lmFit(v,design))
de <- limma::topTable(fit,coef=coef_name,number=Inf,sort.by="none")
de$Gene <- rownames(de)
write.csv(de,file.path(outdir,"Paired_limma_voom_AllGenes.csv"),row.names=FALSE)
de_marker <- de[match(genes,de$Gene),,drop=FALSE]
de_marker <- cbind(Gene=genes,de_marker[,setdiff(names(de_marker),"Gene"),drop=FALSE])
write.csv(de_marker,file.path(outdir,"Paired_limma_voom_MarkerStats.csv"),row.names=FALSE)

# --- Plot annotations -----------------------------------------------
ann_col <- data.frame("Cell type"=cell_type,check.names=FALSE)
rownames(ann_col) <- sample_order
ann_row <- data.frame(
  Category=factor(rep(names(markers),lengths(markers)),levels=names(markers))
)
rownames(ann_row) <- genes
ann_colors <- list(
  "Cell type"=c("UDC"="#C8A1D6","UDC-iN"="#E88980"),
  Category=c("Epithelial"="#CBA26D","Somatic / ECM"="#8F9C73",
             "Proliferation"="#A27B9E","Neuronal"="#648BB4")
)
row_gaps <- head(cumsum(lengths(markers)),-1L)

# Diverging palette: dark blue -> white -> dark red.
heat_colors <- grDevices::colorRampPalette(
  c("#174A8B","#4385BE","#D7E5F2","#FFFFFF",
    "#F3D7D3","#DA7368","#A5232C")
)(201)
heat_breaks <- seq(-color_limit,color_limit,length.out=202)

make_hm <- function(matrix_data,colors,breaks,title,legend_breaks,legend_labels) {
  pheatmap::pheatmap(
    matrix_data,
    cluster_rows=FALSE,cluster_cols=FALSE,
    labels_col=labels,
    annotation_col=ann_col,annotation_row=ann_row,
    annotation_colors=ann_colors,
    gaps_col=length(donors),gaps_row=row_gaps,
    fontsize=12,fontsize_row=14,fontsize_col=10,
    angle_col=45,
    cellwidth=34,cellheight=31,
    show_colnames=TRUE,show_rownames=TRUE,
    border_color=NA,
    color=colors,breaks=breaks,
    legend_breaks=legend_breaks,legend_labels=legend_labels,
    main=title,silent=TRUE
  )
}

# Wider device and an explicit outer margin avoid clipping first/last labels.
save_hm <- function(hm,stem,width=15.8,height=9.5) {
  pdf_file <- file.path(outdir,paste0(stem,".pdf"))
  png_file <- file.path(outdir,paste0(stem,"_600dpi.png"))
  grDevices::pdf(pdf_file,width=width,height=height,
                 onefile=FALSE,useDingbats=FALSE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x=0.51,y=0.52,width=0.93,height=0.90))
  grid::grid.draw(hm$gtable)
  grid::popViewport()
  grDevices::dev.off()
  grDevices::png(png_file,width=width,height=height,
                 units="in",res=600,bg="white")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x=0.51,y=0.52,width=0.93,height=0.90))
  grid::grid.draw(hm$gtable)
  grid::popViewport()
  grDevices::dev.off()
  message("Saved: ",pdf_file," | ",png_file)
}

main_hm <- make_hm(
  z,heat_colors,heat_breaks,"UDC-iN vs UDC",
  c(-1.5,-1,0,1,1.5),c("-1.5","-1","0","1","1.5")
)
save_hm(main_hm,"MainFigure_10Markers_RowZscore")

# Absolute-expression companion (not row-scaled).
range_expr <- range(expr,finite=TRUE)
if (diff(range_expr)==0) range_expr <- range_expr+c(-0.5,0.5)
abs_colors <- grDevices::colorRampPalette(
  c("#F7F7F7","#D4E3EF","#73A2C6","#1D497C")
)(201)
abs_breaks <- seq(range_expr[1],range_expr[2],length.out=202)
ticks <- pretty(range_expr,n=5)
ticks <- ticks[ticks>=range_expr[1] & ticks<=range_expr[2]]
abs_hm <- make_hm(
  expr,abs_colors,abs_breaks,"UDC-iN vs UDC | log2 CPM",
  ticks,format(ticks,trim=TRUE)
)
save_hm(abs_hm,"Supplement_10Markers_log2CPM")

cat("\nDONE:",length(donors),"paired HC donors,",
    length(genes),"markers.\nOutput:",outdir,"\n")
cat("Note: row Z-score colors saturate at +/-",color_limit,
    "; exported Z-scores are UNCLIPPED.\n")
cat("Important: selected markers are illustrative; check paired DE",
    " statistics and expression QC before interpreting.\n")

capture.output(sessionInfo(), file = file.path(outdir, "R_SessionInfo.txt"))
