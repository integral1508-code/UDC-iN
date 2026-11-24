# ============================================================================
# UDC-iN HC donors | age-associated ribosomal proteins
# Sex- and batch-adjusted donor points + cell-type-specific age-regression lines (no CI)
# Same joint paired limma model as supplied age-volcano analysis:
#   ~ Age10 * CellType + Sex + Batch; block = Donor (duplicateCorrelation)
# Selection from ORIGINAL saved UDCiN_Age/All_Proteins.csv (no re-testing)
# UDC-iN age-associated proteins (FDR < 0.05), irrespective of slope direction.
# No GO term membership filter; this is an exploratory selected panel.
# This is exploratory, significance-based selection (not an a priori panel).
# ============================================================================
cran <- c('readxl','readr','dplyr','ggplot2','patchwork','tidyr')
for (p in cran) if (!requireNamespace(p,quietly=TRUE)) install.packages(p,repos='https://cloud.r-project.org')
if (!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org')
for (p in c('limma')) if (!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,ask=FALSE,update=FALSE)
suppressPackageStartupMessages({library(readxl);library(readr);library(dplyr);library(ggplot2);library(limma)})
D <- path.expand(Sys.getenv("PROTEOMICS_DATA_DIR", unset = "~/Desktop"))
# Specify these manually if your filenames or folder locations differ.
proteome_file <- file.path(D,'20230316_iN_batch_All_ys_Proteins_KM edited.xlsx')
mapping_file <- file.path(D,'Proteome_対応表(2).xlsx')
# GO table produced by the original HC age analysis.
# Change this path if the folder is named differently.
output <- file.path(D,'UDCiN_Age_Ribosomal_Figures_From_Original_JointModel')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
for (p in c(proteome_file,mapping_file)) if (!file.exists(p)) stop('Missing input: ',p)
write_csv_out <- function(x,name) readr::write_csv(as.data.frame(x),file.path(output,name),na='')
clean <- function(x) trimws(as.character(x))
# ---------------------------------------------------------------------------
# 1. Reconstruct the same 9 paired HC donors and the same log2 matrix.
# ---------------------------------------------------------------------------
m <- as.data.frame(readxl::read_excel(mapping_file,col_names=FALSE,.name_repair='minimal'))
if(ncol(m)<3) stop('Mapping requires Sample / Diagnosis / Age columns.')
map <- data.frame(Sample=clean(m[[1]]),Dx_map=clean(m[[2]]),Age_map=suppressWarnings(as.numeric(m[[3]])))
map <- map[!is.na(map$Sample)&nzchar(map$Sample),,drop=FALSE]
if(anyDuplicated(map$Sample)) stop('Duplicated mapping sample.')
md <- as.data.frame(readxl::read_excel(proteome_file,sheet='metadata',.name_repair='minimal'))
needed <- c('File Name','Type','ID2','Dx','Age','Sex','batch')
if(!all(needed %in% names(md))) stop('Missing metadata: ',paste(setdiff(needed,names(md)),collapse=', '))
meta <- data.frame(Sample=clean(md[['File Name']]),CellType=clean(md[['Type']]),
                   Donor=clean(md[['ID2']]),Sex=clean(md[['Sex']]),Batch=clean(md[['batch']]))
meta <- dplyr::left_join(meta,map,by='Sample')
meta <- meta[meta$Dx_map=='1.HC' & meta$CellType %in% c('UDC','iN'),,drop=FALSE]
meta$CellType[meta$CellType=='iN'] <- 'UDC-iN'
meta <- meta[order(meta$Age_map,meta$Donor,meta$CellType),,drop=FALSE]
if(nrow(meta)!=18 || length(unique(meta$Donor))!=9) stop('Expected 9 paired HC donors (18 samples).')
if(!all(table(meta$Donor,meta$CellType)==1)) stop('Incomplete donor pairing.')
if(anyNA(meta$Age_map)||anyNA(meta$Sex)||any(!nzchar(meta$Sex))||
   anyNA(meta$Batch)||any(!nzchar(meta$Batch))) stop('Missing covariates.')
if(any(tapply(meta$Age_map,meta$Donor,function(x) length(unique(x)))!=1)) stop('Age differs within donor.')
if(any(tapply(meta$Sex,meta$Donor,function(x) length(unique(x)))!=1)) stop('Sex differs within donor.')
main_sheet <- setdiff(readxl::excel_sheets(proteome_file),'metadata')[1]
prot <- as.data.frame(readxl::read_excel(proteome_file,sheet=main_sheet,.name_repair='minimal'))
if(!all(c('Accession','Gene Symbol') %in% names(prot))) stop('Missing Accession / Gene Symbol.')
if(!all(meta$Sample %in% names(prot))) stop('Missing proteome sample columns.')
valid <- !is.na(prot$Accession)&nzchar(clean(prot$Accession))
prot <- prot[valid,,drop=FALSE]
protein_id <- as.character(seq_len(nrow(prot)))
ab <- as.matrix(data.frame(lapply(prot[meta$Sample],function(x) suppressWarnings(as.numeric(as.character(x)))),check.names=FALSE))
rownames(ab) <- protein_id
colnames(ab) <- meta$Sample
keep <- rowSums(is.finite(ab)&ab>0)==ncol(ab)
mat <- log2(ab[keep,,drop=FALSE])
variance <- apply(mat,1,var)
mat <- mat[is.finite(variance)&variance>0,,drop=FALSE]
# ---------------------------------------------------------------------------
# 2. Identify cytosolic ribosomal proteins only (RPL/RPS; exclude MRPL/MRPS).
# ---------------------------------------------------------------------------
syms <- trimws(as.character(prot[['Gene Symbol']][as.integer(rownames(mat))]))
syms[is.na(syms)] <- ''
# ---------------------------------------------------------------------------
# 3. Read original saved limma results. Refit coefficients ONLY for plotting nuisance-adjusted abundance.
# ---------------------------------------------------------------------------
# Read the ORIGINAL three contrasts, never recalculate FDR or select from a new fit.
analysis_root <- file.path(D,'Proteome_HC_Age_JointModel_Sex_WithBatch')
paths <- c('UDCiN_Age','UDC_Age','Interaction')
if(!dir.exists(analysis_root)) stop('Missing original analysis folder: ',analysis_root)
original <- lapply(paths,function(nm) {
  f <- file.path(analysis_root,nm,'All_Proteins.csv')
  if(!file.exists(f)) stop('Missing original contrast file: ',f)
  z <- readr::read_csv(f,show_col_types=FALSE)
  required <- c('ProteinID','SYMBOL','logFC','P.Value','adj.P.Val')
  if(!all(required %in% names(z))) stop('Missing original result columns in ',f,': ',paste(setdiff(required,names(z)),collapse=', '))
  z <- as.data.frame(z)
  if(anyDuplicated(as.character(z$ProteinID))) stop('Duplicated ProteinID in ',f)
  data.frame(ProteinID=as.character(z$ProteinID),Gene=as.character(z$SYMBOL),
    Contrast=switch(nm,UDCiN_Age='UDC-iN',UDC_Age='UDC',Interaction='Interaction'),
    Slope=as.numeric(z$logFC),P=as.numeric(z$P.Value),FDR=as.numeric(z$adj.P.Val),
    stringsAsFactors=FALSE)
})
results <- dplyr::bind_rows(original)
# The original ProteinID is the 1-based accession-valid row index in the Excel
# after the original filtering. Confirm it still refers to the same protein.
orig_ids <- unique(results$ProteinID)
if(!all(orig_ids %in% rownames(mat))) {
  stop('Original ProteinIDs are not all in the current Excel matrix. Use the exact Excel used for the original analysis.')
}
# Verify ProteinID identity using accession, not display SYMBOL.
# Original SYMBOL may be filled using UniProt mapping, unlike Excel Gene Symbol.
annotation_file <- file.path(analysis_root,'Protein_Annotation.csv')
if(!file.exists(annotation_file)) stop('Missing original Protein_Annotation.csv: ',annotation_file)
original_annotation <- as.data.frame(readr::read_csv(annotation_file,show_col_types=FALSE))
if(!all(c('ProteinID','Accession') %in% names(original_annotation))) stop('Protein_Annotation.csv requires ProteinID and Accession.')
if(anyDuplicated(as.character(original_annotation$ProteinID))) stop('Duplicated original annotation ProteinID.')
original_acc <- as.character(original_annotation$Accession[match(rownames(mat),as.character(original_annotation$ProteinID))])
current_acc <- clean(prot[['Accession']][as.integer(rownames(mat))])
if(anyNA(original_acc) || any(original_acc != current_acc)) stop('Original accession does not match Excel at same ProteinID. Check source version.')
if(length(setdiff(rownames(mat),orig_ids)) || length(setdiff(orig_ids,rownames(mat)))) stop('Protein universe differs from original analysis; use same inputs.')
orig_genes <- original[[1]]$Gene[match(rownames(mat),original[[1]]$ProteinID)]
matched <- !is.na(orig_genes) & nzchar(orig_genes) & nzchar(syms)
symbol_mismatches <- sum(matched & orig_genes != syms)
if(symbol_mismatches) message('Note: ',symbol_mismatches,' SYMBOL annotations differ; ProteinID and Accession identity verified.')
# Select ALL significant cytosolic RPL/RPS genes from the original UDC-iN FDR.
selected_tests <- results %>%
  dplyr::filter(Contrast=='UDC-iN',
                grepl('^(RPL|RPS)[A-Z0-9]+$',Gene),
                is.finite(FDR),FDR<0.05,is.finite(Slope)) %>%
  dplyr::arrange(dplyr::desc(abs(Slope)),FDR,Gene) %>%
  dplyr::distinct(Gene,.keep_all=TRUE)
if(!nrow(selected_tests)) stop('No significant cytosolic RPL/RPS proteins in original UDC-iN age results.')
selected <- selected_tests[,c('Gene','ProteinID','Slope','P','FDR')]
selected$RibosomeClass <- 'Cytosolic'
panel <- selected$Gene
message('Selected ',length(panel),' genes from original analysis (no new FDR calculation).')
readr::write_csv(selected,file.path(output,'01_Selected_All_Significant_Cytosolic_Ribosomal_Proteins.csv'))
results <- dplyr::inner_join(results,selected[,c('ProteinID','Gene')],by=c('ProteinID','Gene'))
if(nrow(results)!=3*nrow(selected)) stop('Some selected proteins are missing a contrast in original results.')
results$Gene <- factor(results$Gene,levels=panel)
results$Contrast <- factor(results$Contrast,levels=c('UDC-iN','UDC','Interaction'))
readr::write_csv(results,file.path(output,'02_Selected_Age_Slopes_and_FDR.csv'))
# Reconstruct original joint design for visualization ONLY.
# This does NOT overwrite any saved P, FDR, or slope.
meta$Age10 <- (meta$Age_map-mean(unique(meta[,c('Donor','Age_map')])$Age_map))/10
meta$CellType <- factor(meta$CellType,levels=c('UDC','UDC-iN'))
meta$Sex <- droplevels(factor(meta$Sex))
meta$Batch <- droplevels(factor(meta$Batch))
design <- model.matrix(~ Age10 * CellType + Sex + Batch,data=meta)
if(qr(design)$rank!=ncol(design)) stop('Rank-deficient design.')
# Use the ORIGINAL donor correlation saved in Analysis_Settings.txt.
settings_file <- file.path(analysis_root,'Analysis_Settings.txt')
if(!file.exists(settings_file)) stop('Missing original Analysis_Settings.txt: ',settings_file)
settings <- readLines(settings_file,warn=FALSE)
corline <- grep('^Consensus correlation:',settings,value=TRUE)
if(length(corline)!=1) stop('Cannot find original consensus correlation in Analysis_Settings.txt.')
original_correlation <- suppressWarnings(as.numeric(trimws(sub('^Consensus correlation:','',corline))))
if(!is.finite(original_correlation) || abs(original_correlation)>=1) stop('Invalid original consensus correlation.')
fit <- limma::lmFit(mat[selected$ProteinID,,drop=FALSE],design,
  block=meta$Donor,correlation=original_correlation)
# Guard against stale/mismatched input: fitted slope must match original CSV.
interaction_col <- grep('^Age10:CellType',colnames(design),value=TRUE)
if(length(interaction_col)!=1) stop('Interaction term not found.')
calc_udc <- fit$coefficients[,'Age10']
calc_in <- calc_udc+fit$coefficients[,interaction_col]
orig_in <- selected$Slope[match(rownames(fit$coefficients),selected$ProteinID)]
if(any(abs(calc_in-orig_in)>1e-4,na.rm=TRUE)) {
  stop('Original slopes do not match source Excel and original donor correlation. Check inputs; no figure saved.')
}
# ---------------------------------------------------------------------------
# 4. Sex/batch adjusted expression; age and cell-type effects are retained.
#    Subtract nuisance fitted contribution and add its common reference mean.
#    This ensures UDC and UDC-iN share the SAME adjustment reference.
# ---------------------------------------------------------------------------
nuisance_cols <- grep('^(Sex|Batch)',colnames(design),value=TRUE)
points <- list()
for(k in seq_len(nrow(selected))) {
  id <- selected$ProteinID[k]; gene <- selected$Gene[k]
  if(length(nuisance_cols)) {
    nu <- as.numeric(design[,nuisance_cols,drop=FALSE] %*% fit$coefficients[id,nuisance_cols])
  } else nu <- rep(0,nrow(meta))
  adjusted <- as.numeric(mat[id,])-nu+mean(nu)
  points[[k]] <- data.frame(Gene=gene,Donor=meta$Donor,CellType=as.character(meta$CellType),
    Age=meta$Age_map,Sex=as.character(meta$Sex),Batch=as.character(meta$Batch),
    RawLog2=as.numeric(mat[id,]),AdjustedLog2=adjusted)
}
points <- dplyr::bind_rows(points)
points$Gene <- factor(points$Gene,levels=panel)
points$CellType <- factor(points$CellType,levels=c('UDC-iN','UDC'))
# Z-score per protein across BOTH cell types (not separately), preserving
# the magnitude of differences between UDC and UDC-iN.
points <- points %>% dplyr::group_by(Gene) %>%
  dplyr::mutate(Z=(AdjustedLog2-mean(AdjustedLog2))/sd(AdjustedLog2)) %>% dplyr::ungroup()
readr::write_csv(points,file.path(output,'03_Adjusted_Donor_Expression.csv'))
# ---------------------------------------------------------------------------
# Figure 1: donor expression versus age, UDC-iN LEFT and UDC RIGHT.
# Red (UDC-iN) and blue (UDC) model-derived slopes; shared y-axis per gene.
# After nuisance adjustment, model slope is unchanged.
# ---------------------------------------------------------------------------
sl <- results[as.character(results$Contrast) %in% c('UDC-iN','UDC'),c('Gene','Contrast','Slope')]
sl$Gene <- as.character(sl$Gene)
sl$CellType <- as.character(sl$Contrast)
line_data <- dplyr::left_join(dplyr::mutate(points,Gene=as.character(Gene),CellType=as.character(CellType)),sl[,c('Gene','CellType','Slope')],by=c('Gene','CellType'))
line_data <- line_data %>% dplyr::group_by(Gene,CellType) %>%
  dplyr::summarise(AgeMin=min(Age),AgeMax=max(Age),AgeMean=mean(Age),
                   MeanExpr=mean(AdjustedLog2),Slope=dplyr::first(Slope),.groups='drop')
line_data <- dplyr::bind_rows(
  data.frame(Gene=line_data$Gene,CellType=line_data$CellType,
             Age=line_data$AgeMin,Fit=line_data$MeanExpr+line_data$Slope*(line_data$AgeMin-line_data$AgeMean)/10),
  data.frame(Gene=line_data$Gene,CellType=line_data$CellType,
             Age=line_data$AgeMax,Fit=line_data$MeanExpr+line_data$Slope*(line_data$AgeMax-line_data$AgeMean)/10))
line_data$Gene <- factor(as.character(line_data$Gene),levels=panel[panel %in% as.character(selected$Gene)])
line_data$CellType <- factor(line_data$CellType,levels=c('UDC-iN','UDC'))
p1 <- ggplot(points,aes(x=Age,y=AdjustedLog2))+
  geom_point(color='#A5A9AD',size=1.6,alpha=.9)+
  geom_line(data=line_data,
            mapping=aes(x=Age,y=Fit,group=interaction(Gene,CellType),color=CellType),
            inherit.aes=FALSE,linewidth=.85)+
  scale_color_manual(values=c('UDC-iN'='#C5524C','UDC'='#3C78A6'),guide='none')+
  facet_grid(Gene~CellType,scales='free_y')+
  labs(x='Donor age (years)',y=expression('Sex- and batch-adjusted log'[2]*' abundance'))+
  theme_classic(base_size=10)+theme(strip.background=element_blank(),
    strip.text.y=element_text(angle=0,face='italic',size=9),
    strip.text.x=element_text(face='bold',size=11),
    panel.spacing.y=grid::unit(.65,'lines'))
ggsave(file.path(output,'Figure_1_Age_Regression_UDCiN_Left_UDC_Right.png'),p1,width=10,height=max(12,1.45*length(panel)),dpi=600,bg='white')
ggsave(file.path(output,'Figure_1_Age_Regression_UDCiN_Left_UDC_Right.pdf'),p1,width=10,height=max(12,1.45*length(panel)))
# ---------------------------------------------------------------------------
# Figure 2: 2-column heatmap, x=actual age (continuous), y=protein.
# Each tile represents ONE donor. Tile widths adapt to local age spacing;
# two panels share age scale and a single z-score color scale.
# ---------------------------------------------------------------------------
ages <- sort(unique(points$Age))
if(anyDuplicated(ages)) {
  # Multiple donors of the same age: position them with a tiny fixed offset
  # solely for display, without changing age in the model or source CSV.
  donor_ages <- unique(points[,c('Donor','Age')]); donor_ages <- donor_ages[order(donor_ages$Age,donor_ages$Donor),]
  donor_ages$PlotAge <- donor_ages$Age
  for(a in unique(donor_ages$Age)) {
    ii <- which(donor_ages$Age==a)
    if(length(ii)>1) donor_ages$PlotAge[ii] <- a+seq(-.2,.2,length.out=length(ii))
  }
  points <- dplyr::left_join(points,donor_ages[,c('Donor','PlotAge')],by='Donor')
} else points$PlotAge <- points$Age
# Use a constant, non-overlapping tile width in years, and label real ages.
distances <- diff(sort(unique(points$PlotAge)))
tile_width <- if(length(distances)) min(5,max(.5,min(distances)*.75)) else 2
points$Gene <- factor(points$Gene,levels=rev(panel[panel %in% selected$Gene]))
lim <- max(1.5,as.numeric(stats::quantile(abs(points$Z),.98,na.rm=TRUE)))
# Independent panel plots with separate headers and a clear central divider.
# A single common color scale is used for both cell types.
make_heat_panel <- function(cell_type, header_fill) {
  dat <- points[as.character(points$CellType)==cell_type,,drop=FALSE]
  ggplot(dat,aes(x=PlotAge,y=Gene,fill=Z))+
    geom_tile(width=tile_width,height=.85,color='white',linewidth=.2)+
    scale_fill_gradient2(low='#3E75A4',mid='white',high='#BD514A',
      midpoint=0,limits=c(-lim,lim),oob=scales::squish,
      name='Adjusted\nZ-score')+
    scale_x_continuous(limits=range(points$Age)+c(-1.2,1.2),
      breaks=pretty(range(points$Age),n=5),expand=expansion(mult=0))+
    labs(x='Donor age (years)',y=NULL,title=cell_type)+
    theme_classic(base_size=12)+
    theme(
      plot.title=element_text(face='bold',size=13,hjust=.5,
        margin=margin(t=9,b=9)),
      plot.title.position='panel',
      plot.background=element_rect(fill='white',colour=NA),
      panel.background=element_rect(fill='white',colour='#BFC4CA',linewidth=.65),
      panel.border=element_rect(fill=NA,colour='#BFC4CA',linewidth=.65),
      axis.text.y=element_text(face='italic',size=10),
      axis.ticks.y=element_blank(),axis.line.y=element_blank(),
      axis.line.x=element_line(linewidth=.4),
      plot.margin=margin(9,12,9,12)
    )+
    theme()
}
# Add distinct colored title bands as independent, non-data plots.
make_header <- function(label, fill) {
  ggplot()+
    annotate('rect',xmin=0,xmax=1,ymin=0,ymax=1,fill=fill)+
    annotate('text',x=.5,y=.5,label=label,fontface='bold',size=5,color='#222222')+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
    theme_void()+theme(plot.margin=margin(0,0,0,0))
}
# Remove duplicate internal title, because each panel has a dedicated header.
heat_left <- make_heat_panel('UDC-iN','#F5DEDB')+labs(title=NULL)+
  theme(legend.position='none')
heat_right <- make_heat_panel('UDC','#DCE9F3')+labs(title=NULL)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position='right')
left_panel <- make_header('UDC-iN','#F5DEDB') / heat_left +
  patchwork::plot_layout(heights=c(.07,1))
right_panel <- make_header('UDC','#DCE9F3') / heat_right +
  patchwork::plot_layout(heights=c(.07,1))
# Blank spacer is intentionally wide to visually separate the cell types.
# Keep each header-and-heatmap pair as one panel to avoid flattening.
p2 <- patchwork::wrap_elements(full=left_panel) +
  patchwork::plot_spacer() +
  patchwork::wrap_elements(full=right_panel) +
  patchwork::plot_layout(ncol=3, widths=c(1,.16,1.18))
ggsave(file.path(output,'Figure_2_Age_Expression_Heatmap.png'),p2,
  width=14.5,height=max(6.5,2.2+.38*length(panel)),dpi=600,bg='white')
ggsave(file.path(output,'Figure_2_Age_Expression_Heatmap.pdf'),p2,
  width=14.5,height=max(6.5,2.2+.38*length(panel)),bg='white')
writeLines(c('Selection uses ORIGINAL UDCiN_Age/All_Proteins.csv FDR < 0.05; original UDC_Age and Interaction CSV used for all slopes/FDR.',
 'Keep ALL significant cytosolic RPL/RPS genes (exclude MRPL/MRPS), ordered by absolute UDC-iN age slope (tie: smaller FDR).',
 'No GO term restriction. This is a significance-selected descriptive panel, not independent validation. The count may differ from 15 when rerun on the current source Excel.',
 'Original 9 paired HC donors; same positive quantified protein universe verified against original results.',
 'Model: ~ Age10 * CellType + Sex + Batch; limma duplicateCorrelation(block=Donor).',
 'Figure 1: donor age on x-axis, adjusted abundance on y-axis; red slope UDC-iN, blue slope UDC. UDC-iN left, UDC right.',
 'Figure 2: sex/batch-adjusted log2 abundances, z-scored per gene across BOTH cell types; separate bordered panels, colored headers, central whitespace.',
 'Actual donor ages on x axis. If identical ages, display positions slightly offset.',
 'BH FDR is read from original full-proteome analysis; not recomputed. Original consensus correlation reused for visualization-only fit.',
 'If multiple significant protein entries share a gene, choose largest absolute original UDC-iN slope, then lowest original FDR.',
 'Age association and between-cell-type slope differences are distinct tests.'),
 file.path(output,'Analysis_Notes.txt'))
message('Completed: ',output)
