# ============================================================================
# PROTEOMICS PCA | Healthy controls only | UDC-iN vs UDC
# Matched to supplied RNA-seq PCA style: paired lines, 95% data ellipses,
# donor age labels, blue UDC / red UDC-iN, no legend, 600-dpi PNG.
# Input: 20230316_iN_batch_All_ys_Proteins_KM edited.xlsx
# Method: log2(TMT Scaled abundance), no DESeq2/VST; all eligible proteins.
# HC selection: Dx 1.HC ONLY; YoungControl excluded.
# Label layout: non-overlapping outside columns; PCA coordinates unchanged.
# Diagnosis and age are verified against Proteome_対応表(2).xlsx.
# ============================================================================


pkgs <- c('readxl','ggplot2','ggrepel','dplyr','tibble')
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE))
  install.packages(p,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(ggrepel)
  library(dplyr); library(tibble)
})
set.seed(123)


# ---- Paths ----
# Set PROTEOMICS_DATA_DIR to the directory containing both input files.
desktop <- path.expand(Sys.getenv('PROTEOMICS_DATA_DIR', unset='~/Desktop'))
input_file <- file.path(desktop,'20230316_iN_batch_All_ys_Proteins_KM edited.xlsx')
if (!file.exists(input_file)) {
  candidates <- list.files(desktop,pattern='^20230316_iN_batch_All_ys_Proteins_KM edited.*\\.xlsx$',full.names=TRUE)
  if (length(candidates)==1) input_file <- candidates
}
if (!file.exists(input_file)) stop('Cannot find the proteomics Excel file in the input directory. Edit input_file.')
mapping_file <- file.path(desktop,'Proteome_対応表(2).xlsx')
if (!file.exists(mapping_file)) {
  candidates <- list.files(desktop,pattern='^Proteome_対応表.*\\.xlsx$',full.names=TRUE)
  if (length(candidates)==1) mapping_file <- candidates
}
if (!file.exists(mapping_file)) stop('Cannot find Proteome_対応表(2).xlsx in the input directory.')
outdir <- file.path(desktop,'Proteome_PCA_HC_UDC_vs_UDCiN_MainFigure')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)


# ---- Authoritative diagnosis mapping from Proteome_対応表(2) ----
# Mapping sheet has NO header: col1 = sample File Name, col2 = Dx,
# col3 = Age. It may include control channels with empty Dx/Age.
map_raw <- as.data.frame(read_excel(mapping_file, col_names=FALSE,
                                    .name_repair='minimal'))
if (ncol(map_raw)<3) stop('Mapping sheet must have >=3 columns.')
map <- data.frame(
  Sample=trimws(as.character(map_raw[[1]])),
  Dx_map=trimws(as.character(map_raw[[2]])),
  Age_map=suppressWarnings(as.numeric(map_raw[[3]])),
  stringsAsFactors=FALSE)
map <- map[!is.na(map$Sample) & nzchar(map$Sample),,drop=FALSE]
if (anyDuplicated(map$Sample)) stop('Duplicated sample names in mapping file.')


meta_raw <- as.data.frame(read_excel(input_file,sheet='metadata',.name_repair='minimal'))
required <- c('File Name','Type','ID2','Patient','Dx','Age','batch')
if (!all(required %in% names(meta_raw)))
  stop('Metadata columns missing: ',paste(setdiff(required,names(meta_raw)),collapse=', '))
meta_all <- meta_raw %>%
  transmute(Sample=trimws(as.character(`File Name`)),
            CellType=trimws(as.character(Type)),
            Donor=trimws(as.character(Patient)),
            ID2=trimws(as.character(ID2)),
            Dx_meta=trimws(as.character(Dx)),
            Age_meta=suppressWarnings(as.numeric(Age)),
            Batch=trimws(as.character(batch)))
if (anyDuplicated(meta_all$Sample)) stop('Duplicated sample names in metadata.')
if (!all(meta_all$Sample %in% map$Sample))
  stop('Some metadata samples are absent from the diagnosis mapping: ',
       paste(setdiff(meta_all$Sample,map$Sample),collapse=' | '))
meta_all <- left_join(meta_all,map,by='Sample')


# Strictly use the diagnosis in the supplied correspondence table.
meta <- meta_all %>%
  filter(CellType %in% c('UDC','iN'),Dx_map=='1.HC') %>%
  mutate(CellType=ifelse(CellType=='iN','UDC-iN','UDC'))
if (!nrow(meta)) stop('No 1.HC UDC/iN samples found in mapping file.')
if (anyNA(meta$Age_map) || anyNA(meta$Donor) || any(!nzchar(meta$Donor)))
  stop('Missing HC age or donor ID.')
if (any(meta$Dx_meta != meta$Dx_map,na.rm=TRUE))
  stop('Dx mismatch between proteome metadata and diagnosis mapping for selected HC.')
if (any(meta$Age_meta != meta$Age_map,na.rm=TRUE))
  stop('Age mismatch between proteome metadata and diagnosis mapping for selected HC.')
meta$Age <- meta$Age_map
meta$Dx <- meta$Dx_map
if (!all(table(meta$Donor,meta$CellType)==1))
  stop('HC donor pairs are incomplete or duplicated; inspect metadata.')
if (any(meta %>% group_by(Donor) %>%
        summarise(n_age=n_distinct(Age),.groups='drop') %>% pull(n_age) != 1))
  stop('Age differs between matched donor samples.')
meta <- meta %>% arrange(Age,Donor,CellType)
meta$CellType <- factor(meta$CellType,levels=c('UDC','UDC-iN'))
meta$Donor <- factor(meta$Donor,levels=unique(meta$Donor))
meta$Label <- paste0('HC. ',format(meta$Age,trim=TRUE,scientific=FALSE),'y')
write.csv(meta,file.path(outdir,'HC_Selected_Sample_Metadata.csv'),row.names=FALSE)
cat('Selected:',n_distinct(meta$Donor),'HC donors and',nrow(meta),'samples\n')
cat('Excluded YoungControl samples:',
    sum(meta_all$Dx_map=='6.YoungControl' & meta_all$CellType %in% c('UDC','iN'),na.rm=TRUE),'\n')


# ---- Read abundance matrix; exact match to metadata File Name ----
all_sheets <- excel_sheets(input_file)
main_sheet <- setdiff(all_sheets,'metadata')[1]
prot <- as.data.frame(read_excel(input_file,sheet=main_sheet,.name_repair='minimal'))
if (!('Accession' %in% names(prot))) stop('Accession column not found.')
if (!all(meta$Sample %in% names(prot))) {
  missing <- setdiff(meta$Sample,names(prot))
  stop('Selected sample columns missing in proteome sheet: ',paste(missing,collapse=' | '))
}
# Original protein table may have repeated accession values; keep unique protein rows
protein_id <- trimws(as.character(prot$Accession))
if (anyNA(protein_id) || any(!nzchar(protein_id))) {
  good_id <- !is.na(protein_id) & nzchar(protein_id)
  prot <- prot[good_id,,drop=FALSE]; protein_id <- protein_id[good_id]
}
if (anyDuplicated(protein_id)) {
  warning('Repeated Accession IDs; appending unique row suffixes (no averaging).')
  protein_id <- make.unique(protein_id,sep='_row')
}
ab <- as.matrix(data.frame(lapply(prot[meta$Sample],function(z)
  suppressWarnings(as.numeric(as.character(z)))),check.names=FALSE))
rownames(ab) <- protein_id
colnames(ab) <- meta$Sample
if (any(ab < 0,na.rm=TRUE)) stop('Negative Scaled abundance values found.')


# ---- Filtering and transformation ----
# TMT scaled abundances may be missing; for a robust primary PCA use proteins
# measured with strictly positive intensity in EVERY selected HC sample.
# This avoids imputation-driven PCA and preserves an all-proteins PCA.
keep <- rowSums(is.finite(ab) & ab>0)==ncol(ab)
ab_complete <- ab[keep,,drop=FALSE]
if (nrow(ab_complete)<100)
  stop('Too few complete positive proteins (',nrow(ab_complete),'). Review missingness.')
log2_mat <- log2(ab_complete)
# Exclude zero-variance proteins; PCA centers proteins, does NOT z-scale them.
variances <- apply(log2_mat,1,var)
log2_mat <- log2_mat[is.finite(variances) & variances>0,,drop=FALSE]
if (nrow(log2_mat)<2) stop('Insufficient variable proteins.')
# Important: no batch centering/correction applied; avoids removing biological
# signal without a batch-effect design review. Input is already TMT Scaled.
pca <- prcomp(t(log2_mat),center=TRUE,scale.=FALSE)
variance <- 100*pca$sdev^2/sum(pca$sdev^2)
pca_df <- as.data.frame(pca$x[,1:2,drop=FALSE]) %>%
  rownames_to_column('Sample') %>% left_join(meta,by='Sample')
stopifnot(!anyNA(pca_df$CellType),nrow(pca_df)==nrow(meta))


# ---- Empirical 95% ellipses (same chi-square method as RNA-seq code) ----
make_ellipse <- function(df,level=0.95,n=300) {
  xy <- as.matrix(df[,c('PC1','PC2')]); mu <- colMeans(xy)
  sigma <- cov(xy); eig <- eigen(sigma,symmetric=TRUE)
  if (any(eig$values<=0)) stop('Ellipse covariance singular for ',as.character(df$CellType[1]))
  theta <- seq(0,2*pi,length.out=n)
  circle <- rbind(cos(theta),sin(theta))
  transform_mat <- eig$vectors %*% diag(sqrt(eig$values),nrow=2)
  coords <- sweep(sqrt(qchisq(level,df=2))*transform_mat %*% circle,1,mu,'+')
  data.frame(PC1=coords[1,],PC2=coords[2,],CellType=as.character(df$CellType[1]))
}
ellipse_df <- bind_rows(lapply(levels(meta$CellType),function(ct)
  make_ellipse(pca_df[pca_df$CellType==ct,,drop=FALSE])))
ellipse_df$CellType <- factor(ellipse_df$CellType,levels=levels(meta$CellType))


# ---- Non-overlapping labels in separate columns (PCA coordinates unchanged) ----
# Labels are placed in left/right columns outside each group's ellipse.
# This is more reliable than ggrepel for tightly clustered paired samples.
cell_colors <- c('UDC'='#3C78A6','UDC-iN'='#CC6859')


spread_y <- function(y, min_gap) {
  if (length(y) <= 1) return(y)
  ord <- order(y)
  yy <- y[ord]
  for (i in 2:length(yy)) yy[i] <- max(yy[i], yy[i-1]+min_gap)
  yy <- yy - mean(yy) + mean(y)
  out <- numeric(length(y)); out[ord] <- yy
  out
}


span_x <- diff(range(c(pca_df$PC1,ellipse_df$PC1)))
span_y <- diff(range(c(pca_df$PC2,ellipse_df$PC2)))
if (span_x <= 0) span_x <- 1
if (span_y <= 0) span_y <- 1


# Assign each donor a distinct label row within its cell type.
# Put UDC labels LEFT of UDC ellipse, UDC-iN labels RIGHT of UDC-iN ellipse.
# Unlike point-repelling labels, fixed columns guarantee no label-label overlap.
label_step <- max(0.075 * span_y, 4.5)
make_labels <- function(ct, side) {
  pts <- pca_df[pca_df$CellType==ct,,drop=FALSE]
  ell <- ellipse_df[ellipse_df$CellType==ct,,drop=FALSE]
  pts <- pts[order(pts$PC2,pts$Age),,drop=FALSE]
  pts$LabelY <- spread_y(pts$PC2,label_step)
  if (side=='left') {
    pts$LabelX <- min(ell$PC1,pts$PC1) - 0.085*span_x
    pts$Hjust <- 1
  } else {
    pts$LabelX <- max(ell$PC1,pts$PC1) + 0.085*span_x
    pts$Hjust <- 0
  }
  pts
}
label_df <- bind_rows(make_labels('UDC','left'),make_labels('UDC-iN','right'))


# Group names ABOVE all sample labels and ellipses; never inside the data cloud.
group_labels <- bind_rows(lapply(levels(meta$CellType),function(ct) {
  pts <- pca_df[pca_df$CellType==ct,,drop=FALSE]
  ell <- ellipse_df[ellipse_df$CellType==ct,,drop=FALSE]
  labs <- label_df[label_df$CellType==ct,,drop=FALSE]
  data.frame(CellType=ct,PC1=mean(pts$PC1),
             PC2=max(c(pts$PC2,ell$PC2,labs$LabelY)) + 0.11*span_y)
}))
group_labels$CellType <- factor(group_labels$CellType,levels=levels(meta$CellType))


# Limits include text anchors, group titles and ellipses; extra horizontal margin
# is reserved for text extending beyond the label anchor points.
all_x <- c(pca_df$PC1,ellipse_df$PC1,label_df$LabelX,group_labels$PC1)
all_y <- c(pca_df$PC2,ellipse_df$PC2,label_df$LabelY,group_labels$PC2)
plot_span_x <- diff(range(all_x)); plot_span_y <- diff(range(all_y))
x_limits <- range(all_x) + c(-0.22,0.22)*plot_span_x
y_limits <- range(all_y) + c(-0.07,0.08)*plot_span_y


# ---- Plot: original RNA-seq palette, ellipses, donor connections ----
p_main <- ggplot(pca_df,aes(x=PC1,y=PC2)) +
  geom_line(aes(group=Donor),color='#C9CDD1',linewidth=0.28,alpha=0.45) +
  geom_polygon(data=ellipse_df,
    aes(x=PC1,y=PC2,group=CellType,fill=CellType),inherit.aes=FALSE,
    alpha=0.045,color=NA) +
  geom_path(data=ellipse_df,
    aes(x=PC1,y=PC2,group=CellType,color=CellType),inherit.aes=FALSE,
    linewidth=0.75) +
  geom_segment(data=label_df,
    aes(x=PC1,y=PC2,xend=LabelX,yend=LabelY),inherit.aes=FALSE,
    linewidth=0.23,color='#C5C9CE') +
  geom_point(aes(fill=CellType),shape=21,size=3.3,stroke=0.55,color='white') +
  geom_text(data=label_df,
    aes(x=LabelX,y=LabelY,label=Label,color=CellType,hjust=Hjust),
    inherit.aes=FALSE,size=3.05,family='sans',show.legend=FALSE) +
  geom_text(data=group_labels,
    aes(x=PC1,y=PC2,label=CellType,color=CellType),inherit.aes=FALSE,
    size=4.6,fontface='bold',family='sans',show.legend=FALSE) +
  scale_color_manual(values=cell_colors) +
  scale_fill_manual(values=cell_colors) +
  labs(x=sprintf('PC1 (%.1f%%)',variance[1]),
       y=sprintf('PC2 (%.1f%%)',variance[2])) +
  coord_cartesian(xlim=x_limits,ylim=y_limits,clip='off') +
  theme_classic(base_size=12) +
  theme(text=element_text(family='sans',color='#202020'),
        axis.title=element_text(size=13,color='#202020'),
        axis.text=element_text(size=10.5,color='#303030'),
        axis.line=element_line(linewidth=0.55,color='#303030'),
        axis.ticks=element_line(linewidth=0.45,color='#303030'),
        legend.position='none',panel.grid=element_blank(),
        plot.margin=margin(30,42,28,42))


# ---- Save ----
png_file <- file.path(outdir,'Proteome_PCA_HC_UDC_vs_UDCiN_600dpi.png')
pdf_file <- file.path(outdir,'Proteome_PCA_HC_UDC_vs_UDCiN.pdf')
ggsave(png_file,p_main,device='png',width=8.5,height=6.5,units='in',dpi=600,bg='white')
ggsave(pdf_file,p_main,device='pdf',width=8.5,height=6.5,units='in',useDingbats=FALSE)
ggsave(file.path(outdir,'Proteome_PCA_HC_Preview_300dpi.png'),p_main,
       device='png',width=8.5,height=6.5,units='in',dpi=300,bg='white')
write.csv(pca_df,file.path(outdir,'PCA_Coordinates.csv'),row.names=FALSE)
write.csv(label_df,file.path(outdir,'PCA_Label_Positions.csv'),row.names=FALSE)
write.csv(data.frame(PC=paste0('PC',seq_along(variance)),VariancePercent=variance),
          file.path(outdir,'PCA_Variance.csv'),row.names=FALSE)
write.csv(data.frame(Accession=rownames(log2_mat)),
          file.path(outdir,'PCA_Proteins_Used.csv'),row.names=FALSE)
write.csv(log2_mat,file.path(outdir,'PCA_Log2ScaledAbundance_Matrix.csv'))
write.csv(data.frame(TotalProteinRows=nrow(ab),CompletePositive=nrow(ab_complete),
                     PCAProteins=nrow(log2_mat),Samples=ncol(log2_mat),
                     Donors=n_distinct(meta$Donor)),
          file.path(outdir,'PCA_QC_Summary.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(outdir,'R_SessionInfo.txt'))
print(p_main)
cat('\n=== PROTEOME HC PCA COMPLETED ===\n')
cat('Donors:',n_distinct(meta$Donor),' Samples:',nrow(meta),
    ' Proteins:',nrow(log2_mat),'\n')
cat('PC1:',round(variance[1],1),'%  PC2:',round(variance[2],1),'%\n')
cat('Output:',outdir,'\n')
cat('Note: PCA is exploratory; TMT batch effects are not additionally corrected.\n')
