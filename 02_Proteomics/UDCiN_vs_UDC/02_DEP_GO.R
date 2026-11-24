# ============================================================================
# PROTEOME | HC ONLY | Paired UDC-iN vs UDC | Volcano + GO BP dot plots
# RNA-seq reference style; HC excludes YoungControl.
# Differential abundance: limma (~ Donor + CellType), BH FDR < 0.05,
#                         abs(log2FC) > 1 (post-hoc effect-size threshold).
# GO BP ORA: UP/DOWN separately, BH FDR < 0.05; top 10 by adjusted P.
# GO dot X = GeneRatio; top-to-bottom order = GeneRatio descending.
# Top 10 terms remain selected by BH-adjusted P; colors UP red, DOWN blue.
# Proteomics input: TMT Scaled abundance; complete positive proteins in all
# 18 HC samples, matching the preceding PCA's protein inclusion policy.
# No extra batch correction: donor fixed effects account for donor-level
# differences, including any batch perfectly confounded with donor.
# ============================================================================


cran <- c('readxl','dplyr','ggplot2','ggrepel','stringr','tibble')
bioc <- c('limma','clusterProfiler','org.Hs.eg.db','AnnotationDbi')
for (p in cran) if (!requireNamespace(p,quietly=TRUE))
  install.packages(p,repos='https://cloud.r-project.org')
if (!requireNamespace('BiocManager',quietly=TRUE))
  install.packages('BiocManager',repos='https://cloud.r-project.org')
for (p in bioc) if (!requireNamespace(p,quietly=TRUE))
  BiocManager::install(p,ask=FALSE,update=FALSE)
suppressPackageStartupMessages({
  library(readxl);library(dplyr);library(ggplot2);library(ggrepel)
  library(stringr);library(tibble);library(limma)
  library(clusterProfiler);library(org.Hs.eg.db);library(AnnotationDbi)
})
set.seed(123)


# ---- 1. Inputs / settings ----
# Set PROTEOMICS_DATA_DIR to the directory containing both input files.
desktop <- path.expand(Sys.getenv('PROTEOMICS_DATA_DIR', unset='~/Desktop'))
find_excel <- function(exact,regex) {
  p <- file.path(desktop,exact)
  if (file.exists(p)) return(p)
  v <- list.files(desktop,pattern=regex,full.names=TRUE)
  if (length(v)!=1) stop('Expected one matching file for ',exact,
                        '; found ',length(v),'. Set path manually.')
  v
}
input_file <- find_excel('20230316_iN_batch_All_ys_Proteins_KM edited.xlsx',
                         '^20230316_iN_batch_All_ys_Proteins_KM edited.*\\.xlsx$')
mapping_file <- find_excel('Proteome_対応表(2).xlsx',
                           '^Proteome_対応表.*\\.xlsx$')
outdir <- file.path(desktop,'Proteome_HC_UDCiN_vs_UDC_Volcano_GO_Top10')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
fdr_cutoff <- 0.05; lfc_cutoff <- 1; go_fdr_cutoff <- 0.05; top_n <- 10
color_up <- '#CC6859';color_down <- '#3C78A6';color_ns <- '#C9CDD1'
write_csv <- function(x,name) write.csv(x,file.path(outdir,name),row.names=FALSE,na='')
save_fig <- function(p,stem,width,height) {
  ggsave(file.path(outdir,paste0(stem,'.png')),p,width=width,height=height,
         units='in',dpi=600,bg='white',limitsize=FALSE)
  ggsave(file.path(outdir,paste0(stem,'.pdf')),p,width=width,height=height,
         units='in',device='pdf',limitsize=FALSE)
}


# ---- 2. Diagnosis mapping and paired HC metadata ----
map_raw <- as.data.frame(read_excel(mapping_file,col_names=FALSE,
                                     .name_repair='minimal'))
if(ncol(map_raw)<3) stop('Mapping must contain sample name, diagnosis, age.')
map <- data.frame(Sample=trimws(as.character(map_raw[[1]])),
                  Dx_map=trimws(as.character(map_raw[[2]])),
                  Age_map=suppressWarnings(as.numeric(map_raw[[3]])))
map <- map[!is.na(map$Sample)&nzchar(map$Sample),,drop=FALSE]
if(anyDuplicated(map$Sample)) stop('Duplicate sample names in mapping.')
meta_raw <- as.data.frame(read_excel(input_file,sheet='metadata',.name_repair='minimal'))
req <- c('File Name','Type','Patient','Dx','Age','batch')
if(!all(req %in% names(meta_raw))) stop('Missing metadata columns: ',
                                       paste(setdiff(req,names(meta_raw)),collapse=', '))
meta_all <- meta_raw %>% transmute(
  Sample=trimws(as.character(`File Name`)),
  CellType=trimws(as.character(Type)),
  Donor=trimws(as.character(Patient)),
  Dx_meta=trimws(as.character(Dx)),
  Age_meta=suppressWarnings(as.numeric(Age)),
  Batch=trimws(as.character(batch))) %>% left_join(map,by='Sample')
if(anyDuplicated(meta_all$Sample)) stop('Duplicate sample names in metadata.')
meta <- meta_all %>% filter(Dx_map=='1.HC',CellType %in% c('UDC','iN')) %>%
  mutate(CellType=ifelse(CellType=='iN','UDC-iN','UDC')) %>%
  arrange(Age_map,Donor,CellType)
if(nrow(meta)!=18 || n_distinct(meta$Donor)!=9)
  stop('Expected 9 HC donor pairs / 18 samples; found ',nrow(meta),' samples.')
if(anyNA(meta$Donor)||anyNA(meta$Age_map)||any(meta$Dx_meta!=meta$Dx_map)||
   any(meta$Age_meta!=meta$Age_map)) stop('Metadata / diagnosis mapping mismatch.')
if(!all(table(meta$Donor,meta$CellType)==1)) stop('Incomplete donor pairing.')
meta$CellType <- factor(meta$CellType,levels=c('UDC','UDC-iN'))
meta$Donor <- factor(meta$Donor,levels=unique(meta$Donor))
write_csv(meta,'01_HC_Selected_Samples.csv')


# ---- 3. Proteome matrix: same filtering as HC PCA ----
main_sheet <- setdiff(excel_sheets(input_file),'metadata')[1]
prot <- as.data.frame(read_excel(input_file,sheet=main_sheet,.name_repair='minimal'))
if(!'Accession' %in% names(prot)) stop('Missing Accession column.')
if(!all(meta$Sample %in% names(prot)))
  stop('Missing proteome columns: ',paste(setdiff(meta$Sample,names(prot)),collapse=' | '))
acc <- trimws(as.character(prot$Accession))
good <- !is.na(acc)&nzchar(acc)
prot <- prot[good,,drop=FALSE];acc <- acc[good]
protein_id <- make.unique(acc,sep='_row')
ab <- as.matrix(data.frame(lapply(prot[meta$Sample],function(z)
  suppressWarnings(as.numeric(as.character(z)))),check.names=FALSE))
rownames(ab) <- protein_id;colnames(ab) <- meta$Sample
if(any(ab<0,na.rm=TRUE)) stop('Negative abundance values.')
keep <- rowSums(is.finite(ab)&ab>0)==ncol(ab)
log_mat <- log2(ab[keep,,drop=FALSE])
variances <- apply(log_mat,1,var)
log_mat <- log_mat[is.finite(variances)&variances>0,,drop=FALSE]
if(nrow(log_mat)<100) stop('Too few complete positive proteins.')


# ---- 4. Paired limma: donor fixed effect, no additional batch correction ----
design <- model.matrix(~Donor+CellType,data=meta)
if(qr(design)$rank!=ncol(design)) stop('Design matrix is not full rank.')
coef_name <- grep('^CellType',colnames(design),value=TRUE)
if(length(coef_name)!=1) stop('Could not identify UDC-iN vs UDC coefficient.')
fit <- limma::lmFit(log_mat,design)
fit <- limma::eBayes(fit)
res <- limma::topTable(fit,coef=coef_name,number=Inf,sort.by='none',
                       adjust.method='BH')
res$ProteinID <- rownames(res)
res$Accession <- sub('_row[0-9]+$','',res$ProteinID)
res$Regulation <- ifelse(!is.na(res$adj.P.Val)&res$adj.P.Val<fdr_cutoff&
                           res$logFC>lfc_cutoff,'UP',
                         ifelse(!is.na(res$adj.P.Val)&res$adj.P.Val<fdr_cutoff&
                                  res$logFC< -lfc_cutoff,'DOWN','NS'))
res$negLog10Padj <- -log10(pmax(res$adj.P.Val,.Machine$double.xmin))


# ---- 5. UniProt accession -> human gene symbols / Entrez IDs ----
# Proteomics uses UniProt accessions, NOT RNA-seq gene symbols.
# Map via org.Hs.eg.db UNIPROT. Keep one Entrez ID per accession (first).
keys <- unique(res$Accession)
keys <- keys[!is.na(keys)&nzchar(keys)]
valid_keys <- intersect(keys,AnnotationDbi::keys(org.Hs.eg.db,keytype='UNIPROT'))
if(length(valid_keys)<10) stop('Very few UniProt accessions match org.Hs.eg.db.')
entrez_map <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=valid_keys,
                column='ENTREZID',keytype='UNIPROT',multiVals='first')
symbol_map <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=valid_keys,
                column='SYMBOL',keytype='UNIPROT',multiVals='first')
res$ENTREZID <- unname(as.character(entrez_map[res$Accession]))
res$SYMBOL <- unname(as.character(symbol_map[res$Accession]))
write_csv(res,'02_Limma_All_Proteins.csv')
write_csv(res[res$Regulation=='UP',,drop=FALSE],'03_Limma_UP.csv')
write_csv(res[res$Regulation=='DOWN',,drop=FALSE],'04_Limma_DOWN.csv')
write_csv(res[res$Regulation!='NS',,drop=FALSE],'05_Limma_All_Significant.csv')


# ---- 6. Volcano ----
vol <- res[is.finite(res$logFC)&is.finite(res$negLog10Padj),,drop=FALSE]
labels <- vol %>% filter(Regulation!='NS') %>%
  group_by(Regulation) %>% arrange(adj.P.Val,desc(abs(logFC)),.by_group=TRUE) %>%
  slice_head(n=8) %>% ungroup() %>%
  mutate(Display=ifelse(!is.na(SYMBOL)&nzchar(SYMBOL),SYMBOL,Accession)) %>%
  distinct(Display,.keep_all=TRUE)
xmax <- max(2,ceiling(max(abs(vol$logFC))))
p_volcano <- ggplot(vol,aes(logFC,negLog10Padj))+
  geom_point(data=filter(vol,Regulation=='NS'),color=color_ns,size=0.85,alpha=0.4)+
  geom_point(data=filter(vol,Regulation!='NS'),aes(color=Regulation),
             size=1.15,alpha=0.8)+
  geom_vline(xintercept=c(-lfc_cutoff,lfc_cutoff),linetype='dashed',
             linewidth=0.35,color='grey55')+
  geom_hline(yintercept=-log10(fdr_cutoff),linetype='dashed',
             linewidth=0.35,color='grey55')+
  ggrepel::geom_text_repel(data=labels,aes(label=Display,color=Regulation),
     size=3,max.overlaps=Inf,max.iter=10000,max.time=5,box.padding=0.4,
     point.padding=0.3,min.segment.length=0,segment.color='grey65',
     segment.size=0.25,seed=123,show.legend=FALSE)+
  scale_color_manual(values=c(UP=color_up,DOWN=color_down,NS=color_ns))+
  scale_x_continuous(limits=c(-xmax,xmax),expand=expansion(mult=0.04))+
  labs(x=expression(log[2]~'fold change (UDC-iN / UDC)'),
       y=expression(-log[10]~'(adjusted P)'))+
  theme_classic(base_size=12)+
  theme(axis.title=element_text(size=13,color='black'),
        axis.text=element_text(size=10.5,color='black'),
        axis.line=element_line(linewidth=0.55),
        axis.ticks=element_line(linewidth=0.45),
        legend.position='none',plot.margin=margin(20,25,20,20))
save_fig(p_volcano,'01_Volcano_HC_UDCiN_vs_UDC',8,6.5)


# ---- 7. GO BP ORA, UP/DOWN separately ----
# GO universe = Entrez IDs of all proteins tested in limma and mapped.
# Multiple proteins can map to one gene; use unique Entrez IDs.
valid_res <- res[!is.na(res$adj.P.Val)&!is.na(res$ENTREZID)&
                 nzchar(res$ENTREZID),,drop=FALSE]
universe_ids <- unique(valid_res$ENTREZID)
up_ids <- unique(valid_res$ENTREZID[valid_res$Regulation=='UP'])
down_ids <- unique(valid_res$ENTREZID[valid_res$Regulation=='DOWN'])
run_go <- function(ids) {
  if(length(ids)<5) return(NULL)
  clusterProfiler::enrichGO(gene=ids,universe=universe_ids,
      OrgDb=org.Hs.eg.db,keyType='ENTREZID',ont='BP',
      pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,
      minGSSize=10,maxGSSize=500,readable=TRUE)
}
extract_go <- function(x) {
  if(is.null(x)) return(data.frame())
  as.data.frame(x)
}
go_up_all <- extract_go(run_go(up_ids))
go_down_all <- extract_go(run_go(down_ids))
go_up_sig <- if(nrow(go_up_all)) filter(go_up_all,p.adjust<go_fdr_cutoff) else go_up_all
go_down_sig <- if(nrow(go_down_all)) filter(go_down_all,p.adjust<go_fdr_cutoff) else go_down_all
write_csv(go_up_all,'06_GO_BP_UP_All.csv')
write_csv(go_down_all,'07_GO_BP_DOWN_All.csv')
write_csv(go_up_sig,'08_GO_BP_UP_Significant.csv')
write_csv(go_down_sig,'09_GO_BP_DOWN_Significant.csv')
ratio_numeric <- function(x) vapply(strsplit(as.character(x),'/',fixed=TRUE),
                                     function(z) as.numeric(z[1])/as.numeric(z[2]),numeric(1))
select_top <- function(df) {
  if(!nrow(df)) return(df)
  df %>% arrange(p.adjust,desc(Count),ID) %>% slice_head(n=top_n) %>%
    mutate(Term=stringr::str_wrap(Description,width=43),
           GeneRatioNumeric=ratio_numeric(GeneRatio),
           Rank=row_number())
}
up_top <- select_top(go_up_sig)
down_top <- select_top(go_down_sig)
write_csv(up_top,'10_GO_BP_UP_Top10.csv')
write_csv(down_top,'11_GO_BP_DOWN_Top10.csv')


# ---- 8. GO dot plot: GeneRatio X, highest GeneRatio at TOP ----
make_dot <- function(df,direction) {
  if(!nrow(df)) {
    message('No significant GO BP terms for ',direction,'; no plot generated.')
    return(NULL)
  }
  # Select top 10 by BH-adjusted P above; display them by GeneRatio descending.
  df <- df %>% arrange(desc(GeneRatioNumeric),p.adjust,ID)
  df$TermID <- factor(df$ID,levels=rev(df$ID))
  term_labels <- stats::setNames(df$Term,df$ID)
  p <- ggplot(df,aes(x=GeneRatioNumeric,y=TermID))+
    geom_point(aes(size=Count,color=p.adjust),alpha=0.95)+
    scale_y_discrete(labels=term_labels)+
    scale_size_continuous(range=c(3.5,9),name='Gene count')+
    scale_color_gradient(low=if(direction=='UP') '#B9473D' else '#235C91',
                         high=if(direction=='UP') '#F0C6B8' else '#BFD8ED',
                         limits=c(0,go_fdr_cutoff),name='Adjusted P')+
    scale_x_continuous(expand=expansion(mult=c(0.05,0.12)))+
    labs(x='Gene ratio',y=NULL)+theme_classic(base_size=12)+
    theme(axis.title.x=element_text(size=12,color='black'),
          axis.text.x=element_text(size=10,color='black'),
          axis.text.y=element_text(size=10.5,color='black'),
          axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          axis.line.x=element_line(linewidth=0.55),
          legend.title=element_text(size=10),legend.text=element_text(size=9),
          legend.position='right',plot.margin=margin(15,20,15,15))
  save_fig(p,paste0('GO_BP_',direction,'_Top10_GeneRatio_Dot'),10,6)
  p
}
p_go_up <- make_dot(up_top,'UP')
p_go_down <- make_dot(down_top,'DOWN')


# ---- 9. Reproducibility / summary ----
summary <- data.frame(Metric=c('HC donors','HC samples','Proteins tested',
  'Proteins UP','Proteins DOWN','Mapped tested Entrez genes','Mapped UP genes',
  'Mapped DOWN genes','GO BP UP FDR<0.05','GO BP DOWN FDR<0.05'),
  Value=c(n_distinct(meta$Donor),nrow(meta),nrow(res),
  sum(res$Regulation=='UP'),sum(res$Regulation=='DOWN'),
  length(universe_ids),length(up_ids),length(down_ids),
  nrow(go_up_sig),nrow(go_down_sig)))
write_csv(summary,'12_Analysis_Summary.csv')
writeLines(c('HC only: Dx_map == 1.HC; YoungControl excluded',
 'Proteins: complete positive Scaled abundance in all 18 HC samples',
 'Transform: log2(Scaled abundance)',
 'Paired limma design: ~ Donor + CellType; UDC-iN / UDC',
 'No additional batch correction; donor effects included',
 'Significant proteins: BH-adjusted P < 0.05 AND abs(log2FC) > 1',
 'GO: Biological Process ORA; UP/DOWN separately',
 'GO background: mapped tested protein genes (unique Entrez IDs)',
 'GO FDR < 0.05; top 10 by adjusted P (ascending)',
 'GO dot X: GeneRatio; top-to-bottom order: GeneRatio descending',
 'GO dot color: UP red gradient, DOWN blue gradient; size: Count',
 'PNG 600 dpi and vector PDF'),file.path(outdir,'Analysis_Settings.txt'))
capture.output(sessionInfo(),file=file.path(outdir,'R_SessionInfo.txt'))
print(summary)
print(p_volcano)
if(!is.null(p_go_up)) print(p_go_up)
if(!is.null(p_go_down)) print(p_go_down)
cat('\nCompleted. Output: ',outdir,'\n',sep='')
