# ============================================================================
# HC-only proteomics | joint paired age x cell-type model
# Nature-style volcano + GO BP UP/DOWN top-10 GeneRatio dot plots
# Age modeled continuously PER 10 YEARS; UP/DOWN = increase/decrease with age.
# ONE sex-adjusted joint model, donor-blocked; all 3 contrasts from SAME fit.
# Batch-adjusted model; no sensitivity analyses.
# Data: 9 adult HC donors per cell type; YoungControl excluded.
# ============================================================================
cran <- c('readxl','dplyr','ggplot2','ggrepel','stringr','statmod')
bioc <- c('limma','clusterProfiler','org.Hs.eg.db','AnnotationDbi')
for(p in cran) if(!requireNamespace(p,quietly=TRUE)) install.packages(p,repos='https://cloud.r-project.org')
if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org')
for(p in bioc) if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,ask=FALSE,update=FALSE)
suppressPackageStartupMessages({library(readxl);library(dplyr);library(ggplot2);library(ggrepel);library(limma);library(clusterProfiler);library(org.Hs.eg.db);library(AnnotationDbi)})
set.seed(123)
desktop <- path.expand(Sys.getenv("PROTEOMICS_DATA_DIR", unset = "~/Desktop"))
find_excel <- function(exact,pattern) {
 p <- file.path(desktop,exact); if(file.exists(p)) return(p)
 v <- list.files(desktop,pattern=pattern,full.names=TRUE)
 if(length(v)!=1) stop('Could not uniquely identify ',exact,' on Desktop. Found: ',length(v))
 v
}
input_file <- find_excel('20230316_iN_batch_All_ys_Proteins_KM edited.xlsx','^20230316_iN_batch_All_ys_Proteins_KM edited.*\\.xlsx$')
mapping_file <- find_excel('Proteome_対応表(2).xlsx','^Proteome_対応表.*\\.xlsx$')
root <- file.path(desktop,'Proteome_HC_Age_JointModel_Sex_WithBatch')
dir.create(root,recursive=TRUE,showWarnings=FALSE)
# Age effect: BH FDR only; no effect-size cutoff.
fdr_cutoff <- 0.05
age_lfc_cutoff <- 0  # Prior criterion: BH FDR < 0.05 only; no slope cutoff
go_fdr_cutoff <- 0.05
top_n <- 10
# Sex- and batch-adjusted joint model.
red <- '#CC6859'; blue <- '#3C78A6'; gray <- '#C9CDD1'
write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,na='')
save_fig <- function(p,folder,stem,w,h) {
 ggsave(file.path(folder,paste0(stem,'.png')),p,width=w,height=h,units='in',dpi=600,bg='white')
 ggsave(file.path(folder,paste0(stem,'.pdf')),p,width=w,height=h,units='in',device='pdf')
}
# 1. Read mapping (three columns, no header) and validate diagnoses.
map_raw <- as.data.frame(read_excel(mapping_file,col_names=FALSE,.name_repair='minimal'))
if(ncol(map_raw)<3) stop('Mapping must have sample, diagnosis, age columns.')
map <- data.frame(Sample=trimws(as.character(map_raw[[1]])),
 Dx_map=trimws(as.character(map_raw[[2]])),Age_map=suppressWarnings(as.numeric(map_raw[[3]])))
map <- map[!is.na(map$Sample)&nzchar(map$Sample),,drop=FALSE]
if(anyDuplicated(map$Sample)) stop('Duplicated sample names in mapping.')
meta_raw <- as.data.frame(read_excel(input_file,sheet='metadata',.name_repair='minimal'))
req <- c('File Name','Type','ID2','Dx','Age','Sex','batch')
if(!all(req %in% names(meta_raw))) stop('Missing metadata columns: ',paste(setdiff(req,names(meta_raw)),collapse=', '))
meta <- meta_raw %>% transmute(Sample=trimws(as.character(`File Name`)),
 CellType=trimws(as.character(Type)),Donor=trimws(as.character(ID2)),
 Dx_meta=trimws(as.character(Dx)),Age_meta=suppressWarnings(as.numeric(Age)),
  Sex=trimws(as.character(Sex)),
  Batch=trimws(as.character(batch))) %>% left_join(map,by='Sample') %>%
 filter(Dx_map=='1.HC',CellType %in% c('UDC','iN')) %>%
 mutate(CellType=ifelse(CellType=='iN','UDC-iN','UDC')) %>% arrange(Age_map,Donor,CellType)
if(nrow(meta)!=18 || n_distinct(meta$Donor)!=9) stop('Expected 9 adult HC pairs (18 samples); found ',nrow(meta))
if(anyNA(meta$Sex)||any(!nzchar(meta$Sex))||anyNA(meta$Age_map)||anyNA(meta$Donor)) stop('Missing sex, age or donor information.')
if(!all(table(meta$Donor,meta$CellType)==1)) stop('Incomplete HC donor pairs.')
if(any(tapply(meta$Age_map,meta$Donor,function(z)length(unique(z)))!=1)) stop('Age mismatch within donor.')
if(any(tapply(meta$Sex,meta$Donor,function(z)length(unique(z)))!=1)) stop('Sex mismatch within donor.')
write_csv(meta,file.path(root,'HC_Included_Samples.csv'))
# 2. Proteome matrix: complete positive proteins in all 18 samples, matching previous analysis.
main_sheet <- setdiff(excel_sheets(input_file),'metadata')[1]
prot <- as.data.frame(read_excel(input_file,sheet=main_sheet,.name_repair='minimal'))
if(!all(c('Accession','Gene Symbol') %in% names(prot))) stop('Accession or Gene Symbol column missing.')
if(!all(meta$Sample %in% names(prot))) stop('Proteome sample columns missing: ',paste(setdiff(meta$Sample,names(prot)),collapse='; '))
acc <- trimws(as.character(prot$Accession)); valid <- !is.na(acc)&nzchar(acc)
prot <- prot[valid,,drop=FALSE];acc <- acc[valid]
protein_id <- as.character(seq_len(nrow(prot)))
ab <- as.matrix(data.frame(lapply(prot[meta$Sample],function(z) suppressWarnings(as.numeric(as.character(z)))),check.names=FALSE))
rownames(ab) <- protein_id;colnames(ab) <- meta$Sample
# Negative and zero abundance are treated as missing, as in the earlier validation script.
keep <- rowSums(is.finite(ab)&ab>0)==ncol(ab)
log_mat <- log2(ab[keep,,drop=FALSE]); vv <- apply(log_mat,1,var)
log_mat <- log_mat[is.finite(vv)&vv>0,,drop=FALSE]
if(nrow(log_mat)<100) stop('Too few proteins after filtering.')
# UniProt -> Entrez/SYMBOL mapping (one annotation per protein; GO de-duplicated by gene)
accessions <- acc[as.integer(rownames(log_mat))]
keys <- intersect(unique(accessions),AnnotationDbi::keys(org.Hs.eg.db,keytype='UNIPROT'))
if(length(keys)<10) stop('Very few UniProt IDs mapped.')
entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=keys,column='ENTREZID',keytype='UNIPROT',multiVals='first')
symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=keys,column='SYMBOL',keytype='UNIPROT',multiVals='first')
ratio_numeric <- function(x) vapply(strsplit(as.character(x),'/',fixed=TRUE),function(z) as.numeric(z[1])/as.numeric(z[2]),numeric(1))
make_dot <- function(df,direction,folder,cell) {
 if(!nrow(df)) {message(cell,' ',direction,': no significant GO terms');return(NULL)}
 df <- df %>% arrange(desc(GeneRatioNumeric),p.adjust,ID)
 df$TermID <- factor(df$ID,levels=rev(df$ID))
 term_labels <- stats::setNames(df$Term,df$ID)
 p <- ggplot(df,aes(GeneRatioNumeric,TermID))+
   geom_point(aes(size=Count,color=p.adjust),alpha=.95)+
   scale_y_discrete(labels=term_labels)+scale_size_continuous(range=c(3.5,9),name='Gene count')+
   scale_color_gradient(low=if(direction=='UP') '#B9473D' else '#235C91',
     high=if(direction=='UP') '#F0C6B8' else '#BFD8ED',
     limits=c(0,go_fdr_cutoff),name='Adjusted P')+
   scale_x_continuous(expand=expansion(mult=c(.05,.12)))+
   labs(title=paste0(switch(cell, UDC_Age='UDC', UDCiN_Age='UDC-iN', Interaction='Interaction'), ' — ', direction),x='Gene ratio',y=NULL)+theme_classic(base_size=12)+
   theme(plot.title=element_text(size=16,face='plain',hjust=0),axis.title.x=element_text(size=12),axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10.5),axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),axis.line.x=element_line(linewidth=.55),
    legend.position='right',plot.margin=margin(15,20,15,15))
 plot_height <- max(6, 1.8 + 0.31 * nrow(df))
 save_fig(p,folder,paste0('GO_BP_',direction,'_AllSignificant_GeneRatio'),10,plot_height)
 p
}
# ============================================================================
# ONE joint model: ~ Age10 * CellType + Sex, block = Donor
# 3 contrasts from same limma fit; batch-adjusted; no sensitivity analyses.
# ============================================================================
meta$Age10 <- (meta$Age_map - mean(unique(meta[,c('Donor','Age_map')])$Age_map))/10
meta$CellType <- factor(meta$CellType,levels=c('UDC','UDC-iN'))
meta$Sex <- droplevels(factor(meta$Sex))
if(anyNA(meta$Batch) || any(!nzchar(meta$Batch))) stop('Missing batch in selected HC samples.')
meta$Batch <- droplevels(factor(meta$Batch))
design <- model.matrix(~ Age10 * CellType + Sex + Batch,data=meta)
if(qr(design)$rank != ncol(design)) stop('Rank deficient design: ',paste(colnames(design),collapse=', '))
if(nrow(design)-ncol(design)<2) stop('Insufficient residual degrees of freedom.')
coef_inter <- grep('^Age10:CellType',colnames(design),value=TRUE)
if(length(coef_inter)!=1L) stop('Interaction coefficient not found.')
# donor correlation must be estimable
corfit <- limma::duplicateCorrelation(log_mat,design,block=meta$Donor)
if(!is.finite(corfit$consensus.correlation)) stop('Could not estimate within-donor correlation.')
fit <- limma::lmFit(log_mat,design,block=meta$Donor,
                    correlation=corfit$consensus.correlation)
ct <- matrix(0,nrow=ncol(design),ncol=3,
             dimnames=list(colnames(design),c('UDC_Age','UDCiN_Age','Interaction')))
ct['Age10','UDC_Age'] <- 1
ct['Age10','UDCiN_Age'] <- 1
ct[coef_inter,'UDCiN_Age'] <- 1
ct[coef_inter,'Interaction'] <- 1
fit <- limma::contrasts.fit(fit,contrasts=ct)
fit <- limma::eBayes(fit,trend=TRUE,robust=FALSE)
write_csv(data.frame(Sample=meta$Sample,design,check.names=FALSE),
          file.path(root,'Design_Matrix.csv'))
writeLines(c('Model: ~ Age10 * CellType + Sex + Batch',
             'Block: Donor (limma duplicateCorrelation)',
             paste('Consensus correlation:',corfit$consensus.correlation),
             'Batch: adjusted as a fixed effect',
             'Age10: centered age / 10 years',
             'UDC_Age = Age10',
             paste0('UDCiN_Age = Age10 + ',coef_inter),
             paste0('Interaction = ',coef_inter,' (UDC-iN minus UDC)'),
             'Primary protein universe: complete positive finite values in all 18 HC samples',
             'No imputation; exclude zero-variance proteins',
             'eBayes(trend=TRUE,robust=FALSE); BH FDR < 0.05; no FC cutoff',
             'GO BP ORA on mapped tested Entrez genes, up and down separately',
             'GO dot plots include ALL significant GO terms (BH FDR < 0.05), displayed in descending GeneRatio; Top10 CSV retained'),
           file.path(root,'Analysis_Settings.txt'))
# Original Gene Symbol is preferred; UniProt mapping used as fallback.
annotation <- data.frame(ProteinID=rownames(log_mat),
                         Accession=acc[as.integer(rownames(log_mat))],
                         OriginalSymbol=trimws(as.character(prot[['Gene Symbol']][as.integer(rownames(log_mat))])),
                         stringsAsFactors=FALSE)
annotation$SYMBOL <- annotation$OriginalSymbol
annotation$SYMBOL[is.na(annotation$SYMBOL)|!nzchar(annotation$SYMBOL)|
                    grepl('[,;/|]',annotation$SYMBOL)] <- NA_character_
# map SYMBOL -> ENTREZID, avoid treating multi-gene annotations as single genes
all_symbols <- unique(na.omit(annotation$SYMBOL))
symbol_map <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=all_symbols,
                                  column='ENTREZID',keytype='SYMBOL',multiVals='first')
annotation$ENTREZID <- unname(as.character(symbol_map[annotation$SYMBOL]))
# fallback from accession when the original gene symbol is missing/unmappable
missing_entrez <- is.na(annotation$ENTREZID)|!nzchar(annotation$ENTREZID)
annotation$ENTREZID[missing_entrez] <- unname(as.character(entrez[annotation$Accession[missing_entrez]]))
missing_symbol <- is.na(annotation$SYMBOL)|!nzchar(annotation$SYMBOL)
annotation$SYMBOL[missing_symbol] <- unname(as.character(symbol[annotation$Accession[missing_symbol]]))
write_csv(annotation,file.path(root,'Protein_Annotation.csv'))
summaries <- list()
for(contrast in colnames(ct)) {
  folder <- file.path(root,contrast)
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  res <- limma::topTable(fit,coef=contrast,number=Inf,sort.by='none',adjust.method='BH')
  res$ProteinID <- rownames(res)
  res <- merge(annotation,res,by='ProteinID',sort=FALSE)
  res <- res[match(rownames(log_mat),res$ProteinID),,drop=FALSE]
  res$Regulation <- ifelse(is.finite(res$adj.P.Val)&res$adj.P.Val<fdr_cutoff & res$logFC>0,'UP',
                    ifelse(is.finite(res$adj.P.Val)&res$adj.P.Val<fdr_cutoff & res$logFC<0,'DOWN','NS'))
  res$negLog10Padj <- -log10(pmax(res$adj.P.Val,.Machine$double.xmin))
  write_csv(res,file.path(folder,'All_Proteins.csv'))
  write_csv(res[res$Regulation=='UP',,drop=FALSE],file.path(folder,'UP_Proteins.csv'))
  write_csv(res[res$Regulation=='DOWN',,drop=FALSE],file.path(folder,'DOWN_Proteins.csv'))
  # Volcano: log2 slope per decade; interaction is difference in slopes
  vol <- res[is.finite(res$logFC)&is.finite(res$negLog10Padj),,drop=FALSE]
  lab <- vol %>% filter(Regulation!='NS') %>% group_by(Regulation) %>%
    arrange(adj.P.Val,desc(abs(logFC)),.by_group=TRUE) %>% slice_head(n=8) %>%
    ungroup() %>% mutate(Display=ifelse(!is.na(SYMBOL)&nzchar(SYMBOL),SYMBOL,Accession)) %>%
    distinct(Display,.keep_all=TRUE)
  xmax <- max(.5,ceiling(max(abs(vol$logFC))*2)/2)
  x_label <- if(contrast=='Interaction')
    expression(Delta~log[2]~'age slope per 10 years (UDC-iN - UDC)') else
    expression(log[2]~'abundance change per 10 years')
  p <- ggplot(vol,aes(logFC,negLog10Padj))+
    geom_point(data=filter(vol,Regulation=='NS'),color=gray,size=.85,alpha=.4)+
    geom_point(data=filter(vol,Regulation!='NS'),aes(color=Regulation),size=1.15,alpha=.8)+
    geom_vline(xintercept=0,linetype='dashed',linewidth=.35,color='grey55')+
    geom_hline(yintercept=-log10(fdr_cutoff),linetype='dashed',linewidth=.35,color='grey55')+
    ggrepel::geom_text_repel(data=lab,aes(label=Display,color=Regulation),size=3,
      max.overlaps=Inf,max.iter=10000,max.time=5,box.padding=.4,point.padding=.3,
      min.segment.length=0,segment.color='grey65',segment.size=.25,seed=123,show.legend=FALSE)+
    scale_color_manual(values=c(UP=red,DOWN=blue,NS=gray))+
    scale_x_continuous(limits=c(-xmax,xmax),expand=expansion(mult=.04))+
    labs(title=switch(contrast, UDC_Age='UDC', UDCiN_Age='UDC-iN', Interaction='Interaction'),x=x_label,
         y=expression(-log[10]~'(BH-adjusted P)'))+
    theme_classic(base_size=12)+
    theme(axis.title=element_text(size=13),axis.text=element_text(size=10.5),
          axis.line=element_line(linewidth=.55),axis.ticks=element_line(linewidth=.45),
          legend.position='none',plot.margin=margin(20,25,20,20))
  save_fig(p,folder,paste0('Volcano_',contrast),8,6.5)
  valid_res <- res[is.finite(res$adj.P.Val)&!is.na(res$ENTREZID)&nzchar(res$ENTREZID),,drop=FALSE]
  universe <- unique(valid_res$ENTREZID)
  for(direction in c('UP','DOWN')) {
    ids <- unique(valid_res$ENTREZID[valid_res$Regulation==direction])
    all <- data.frame()
    if(length(ids)>=5L) {
      eg <- clusterProfiler::enrichGO(gene=ids,universe=universe,
          OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType='ENTREZID',ont='BP',
          pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,
          minGSSize=10,maxGSSize=500,readable=TRUE)
      if(!is.null(eg)) all <- as.data.frame(eg)
    }
    sig <- if(nrow(all)) all[is.finite(all$p.adjust)&all$p.adjust<go_fdr_cutoff,,drop=FALSE] else all
    top <- if(nrow(sig)) sig %>% arrange(p.adjust,desc(Count),ID) %>%
      slice_head(n=top_n) %>% mutate(Term=stringr::str_wrap(Description,width=43),
                                    GeneRatioNumeric=ratio_numeric(GeneRatio)) else sig
    write_csv(all,file.path(folder,paste0('GO_BP_',direction,'_All.csv')))
    write_csv(sig,file.path(folder,paste0('GO_BP_',direction,'_Significant.csv')))
    write_csv(top,file.path(folder,paste0('GO_BP_',direction,'_Top10.csv')))
    plot_all <- if(nrow(sig)) sig %>% arrange(p.adjust,desc(Count),ID) %>%
      mutate(Term=stringr::str_wrap(Description,width=43),
             GeneRatioNumeric=ratio_numeric(GeneRatio)) else sig
    make_dot(plot_all,direction,folder,contrast)
  }
  summaries[[contrast]] <- data.frame(Contrast=contrast,N_Donors=9,N_Samples=18,
    Proteins=nrow(res),UP=sum(res$Regulation=='UP'),DOWN=sum(res$Regulation=='DOWN'),
    FDR005=sum(res$Regulation!='NS'),Correlation=corfit$consensus.correlation)
}
summary <- dplyr::bind_rows(summaries)
write_csv(summary,file.path(root,'Analysis_Summary.csv'))
capture.output(sessionInfo(),file=file.path(root,'R_SessionInfo.txt'))
print(summary)
cat('\nCompleted. Output:',root,'\n')
if(Sys.info()[['sysname']]=='Darwin') system2('open',root,wait=FALSE)
