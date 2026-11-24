# Donor exclusion: The 72-year-old donor originally classified as AD was excluded from the analyses because the updated clinical diagnosis indicated that the donor did not have AD. Both matched UDC and UDC-iN samples from this donor were excluded, leaving 11 AD UDC-iN donors.

# ============================================================================

# UDC-iN: AD vs HC leading-edge RNA splicing Top 20 -> adjusted 3-group expression

# No rerun of GSEA or limma; MCI does not influence selection.

# ============================================================================

cran <- c('readxl','readr','dplyr','ggplot2')

for (p in cran) if (!requireNamespace(p,quietly=TRUE)) install.packages(p,repos='https://cloud.r-project.org')

if (!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org')

for (p in c('clusterProfiler','org.Hs.eg.db')) if (!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,ask=FALSE,update=FALSE)

library(dplyr)

D <- path.expand(Sys.getenv("PROTEOME_DATA_DIR", unset = "~/Desktop"))

if (!dir.exists(D)) stop('Desktop not found')

# Optional: specify exact input paths here; leave NA to find newest matching Desktop file.

proteome_file <- file.path(D, "20230316_iN_batch_All_ys_Proteins_KM edited.xlsx")

mapping_file <- file.path(D, "Proteome_対応表(2).xlsx")

ranking_file <- NA_character_

gsea_file <- NA_character_

find_input <- function(explicit,pattern,label) {

  if (!is.na(explicit)) {if (!file.exists(explicit)) stop(label,' missing: ',explicit); return(explicit)}

  folders <- c(D,list.dirs(D,recursive=FALSE,full.names=TRUE))

  hits <- unique(unlist(lapply(folders,function(z) list.files(z,pattern=pattern,full.names=TRUE,ignore.case=TRUE))))

  if (!length(hits)) stop('Cannot locate ',label,' on Desktop. Set its path above.')

  hits <- hits[order(file.info(hits)$mtime,decreasing=TRUE)]

  message(label,': ',hits[1]); hits[1]

}

P <- find_input(proteome_file,'^20230316.*Proteins.*edited.*\\.xlsx$','Proteome')

C <- find_input(mapping_file,'^Proteome_対応表.*\\.xlsx$','Diagnosis mapping')

ad11_dir <- file.path(D,'UDCiN_GO_BP_GSEA_MeanT_WithBatch_GlobalBH3Contrasts_nPerm10000_AD11')

if (!dir.exists(ad11_dir)) stop('AD11 GSEA output folder missing: ',ad11_dir)

R <- if (!is.na(ranking_file)) ranking_file else file.path(ad11_dir,'05_Ranked_Genes_MeanT_AD_vs_HC.csv')

if (!file.exists(R)) stop('AD11 ranking CSV missing: ',R)

G <- if (!is.na(gsea_file)) gsea_file else file.path(ad11_dir,'07_GSEA_GO_BP_AD_vs_HC_FDR005.csv')

if (!file.exists(G)) stop('AD11 GSEA CSV missing: ',G)

O <- file.path(D,'UDCiN_ADvsHC_Splicing_7GO_LeadingEdge_Top20_AD11')

dir.create(O,recursive=TRUE,showWarnings=FALSE)

W <- function(x,f) readr::write_csv(as.data.frame(x),file.path(O,f),na='')

clean <- function(x) trimws(as.character(x))

# 1. AD vs HC leading-edge UNION for the seven significant negative GO terms

terms <- c(

 'RNA splicing, via transesterification reactions with bulged adenosine as nucleophile',

 'mRNA splicing, via spliceosome',

 'RNA splicing, via transesterification reactions',

 'mRNA processing',

 'RNA splicing',

 'regulation of mRNA splicing, via spliceosome',

 'regulation of mRNA processing')

go <- readr::read_csv(G,show_col_types=FALSE)

req <- c('ID','Description','NES','core_enrichment','p.adjust_global_3contrasts')

if (!all(req %in% names(go))) stop('GSEA CSV lacks: ',paste(setdiff(req,names(go)),collapse=', '))

g <- go[go$Description %in% terms,,drop=FALSE]

if (nrow(g)!=7 || anyDuplicated(g$Description) || !setequal(g$Description,terms)) stop('All seven specified GO terms must be present')

if (any(!is.finite(g$NES) | g$NES>=0 | !is.finite(g$p.adjust_global_3contrasts) | g$p.adjust_global_3contrasts>=.05)) stop('Expected seven negative GO terms with global FDR < 0.05')

W(g,'01_Seven_Splicing_GO_ADvsHC.csv')

edges <- dplyr::bind_rows(lapply(seq_len(nrow(g)),function(i) {

 ids <- unique(trimws(strsplit(as.character(g$core_enrichment[i]),'/',fixed=TRUE)[[1]]))

 ids <- ids[nzchar(ids)]

 data.frame(ENTREZID=ids,GO_ID=g$ID[i],GO_Term=g$Description[i])

}))

if (!nrow(edges)) stop('core_enrichment is empty')

W(edges,'02_LeadingEdge_Long.csv')

edge_union <- edges %>% dplyr::group_by(ENTREZID) %>%

 dplyr::summarise(N_GO_Terms=dplyr::n_distinct(GO_ID),GO_IDs=paste(sort(unique(GO_ID)),collapse=';'),.groups='drop')

W(edge_union,'03_LeadingEdge_Union.csv')

# 2. Top 20 = most NEGATIVE AD vs HC GSEA RankScore, not top 20 by FDR.

r <- readr::read_csv(R,show_col_types=FALSE)

req <- c('ENTREZID','RankScore','ProteinRows','GeneSymbols')

if (!all(req %in% names(r))) stop('Ranking CSV lacks: ',paste(setdiff(req,names(r)),collapse=', '))

r$ENTREZID <- as.character(r$ENTREZID)

candidates <- dplyr::inner_join(edge_union,r,by='ENTREZID') %>%

 dplyr::filter(is.finite(RankScore),RankScore<0) %>%

 dplyr::arrange(RankScore,ENTREZID)

W(candidates,'04_All_LeadingEdge_Candidates_Ranked.csv')

if (nrow(candidates)<20) stop('Only ',nrow(candidates),' negative leading-edge genes')

top <- head(candidates,20)

top$Top20_Rank <- seq_len(nrow(top))

top <- top[,c('Top20_Rank','ENTREZID','GeneSymbols','RankScore','ProteinRows','N_GO_Terms','GO_IDs')]

W(top,'05_Selected_Top20.csv')

# 3. Match iN donor metadata and diagnosis exactly as original GSEA code.

cl <- as.data.frame(readxl::read_excel(C,col_names=FALSE,.name_repair='minimal'))

if (ncol(cl)<2) stop('Mapping needs two columns')

cl <- data.frame(Sample=clean(cl[[1]]),Label=clean(cl[[2]]))

cl$Dx <- dplyr::case_when(cl$Label=='1.HC'~'HC',

 cl$Label %in% c('1.MCIduetoAD','2.MCIduetoAD')~'MCI',

 cl$Label %in% c('2.AD','3.AD')~'AD',TRUE~NA_character_)

cl <- cl[!is.na(cl$Dx) & !is.na(cl$Sample) & nzchar(cl$Sample),]

if (anyDuplicated(cl$Sample)) stop('Duplicate mapping sample')

md <- as.data.frame(readxl::read_excel(P,sheet='metadata',.name_repair='minimal'))

req <- c('File Name','Type','ID2','Age','Sex','batch')

if (!all(req %in% names(md))) stop('Missing metadata columns')

meta <- data.frame(Sample=clean(md[['File Name']]),CellType=clean(md[['Type']]),

 Donor=clean(md[['ID2']]),Age=suppressWarnings(as.numeric(md[['Age']])),

 Sex=clean(md[['Sex']]),Batch=clean(md[['batch']]))

if (anyDuplicated(meta$Sample)) stop('Duplicate metadata sample')

meta <- merge(meta,cl,by='Sample',sort=FALSE)

# Exclude BOTH samples of the same 72-year-old AD donor before iN selection.

exclude_names <- c(

  UDCiN='Abundances (Scaled): F3: 131C, Sample, 03_AD__iN, 2_iN, 3_AD',

  UDC='Abundances (Scaled): F3: 127C, Sample, 3_AD_UDC, 1_UDC, 3_AD'

)

excluded_records <- meta[match(unname(exclude_names),meta$Sample),,drop=FALSE]

if (nrow(excluded_records)!=2L || anyNA(excluded_records$Sample) ||

    !identical(as.character(excluded_records$Sample),unname(exclude_names)))

  stop('The two specified AD samples were not uniquely found in metadata.')

if (any(excluded_records$Dx!='AD') ||

    any(!is.finite(excluded_records$Age) | excluded_records$Age!=72) ||

    length(unique(excluded_records$Donor))!=1L ||

    !identical(as.character(excluded_records$CellType),c('iN','UDC')))

  stop('Excluded samples are not a matched AD 72y UDC-iN/UDC donor pair.')

W(excluded_records,'00_Excluded_AD_72y_UDC_and_UDCiN.csv')

meta <- meta[!(meta$Sample %in% unname(exclude_names)),,drop=FALSE]

meta <- meta[meta$CellType=='iN' & is.finite(meta$Age) &

 !is.na(meta$Sex) & nzchar(meta$Sex) & !is.na(meta$Donor) & nzchar(meta$Donor) &

 !is.na(meta$Batch) & nzchar(meta$Batch),,drop=FALSE]

meta$Dx <- factor(meta$Dx,levels=c('HC','MCI','AD'))

if (anyDuplicated(meta$Donor)) stop('Repeated donor')

if (any(table(meta$Dx)<2)) stop('Insufficient donors')

if (sum(meta$Dx=='AD')!=11L) stop('Expected 11 AD iN donors after exclusion; got ',sum(meta$Dx=='AD'))

W(meta,'06_Selected_iN_Metadata.csv')

print(table(meta$Dx))

# 4. Original measured TMT intensities; same complete-case filter as GSEA.

sheets <- readxl::excel_sheets(P)

dat_sheet <- sheets[which(tolower(sheets)!='metadata')[1]]

raw <- as.data.frame(readxl::read_excel(P,sheet=dat_sheet,.name_repair='minimal'),check.names=FALSE)

if (!all(c('Accession','Gene Symbol') %in% names(raw))) stop('Missing protein annotations')

if (!all(meta$Sample %in% names(raw))) stop('Missing sample columns')

X0 <- as.matrix(raw[,meta$Sample,drop=FALSE]); suppressWarnings(storage.mode(X0)<-'numeric')

rownames(X0) <- as.character(seq_len(nrow(raw)))

keep <- rowSums(is.finite(X0)&X0>0)==ncol(X0)

X <- log2(X0[keep,,drop=FALSE])

sd0 <- apply(X,1,stats::sd)

X <- X[is.finite(sd0)&sd0>0,,drop=FALSE]

row_id <- as.integer(rownames(X))

anno <- data.frame(RowID=row_id,Accession=clean(raw[['Accession']][row_id]),

 GeneSymbol=clean(raw[['Gene Symbol']][row_id]))

symbols <- unique(anno$GeneSymbol[!is.na(anno$GeneSymbol)&nzchar(anno$GeneSymbol)&

 !grepl('[,;/|]',anno$GeneSymbol)])

mapped <- suppressMessages(clusterProfiler::bitr(symbols,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db::org.Hs.eg.db))

if (is.null(mapped) || !nrow(mapped)) stop('SYMBOL to ENTREZID mapping failed')

entries <- dplyr::inner_join(anno,mapped,by=c('GeneSymbol'='SYMBOL')) %>%

 dplyr::distinct(ENTREZID,RowID,.keep_all=TRUE)

entries$ENTREZID <- as.character(entries$ENTREZID)

entries <- dplyr::inner_join(entries,top[,c('ENTREZID','Top20_Rank')],by='ENTREZID')

if (!all(top$ENTREZID %in% entries$ENTREZID)) stop('Top20 not all mapped to original protein rows')

W(entries,'07_Top20_Original_Protein_Entries.csv')

# Mean of measured log2 protein entries per gene/donor, consistent with

# multiple-protein-per-gene aggregation in the ranking (but NOT a t-test).

long <- data.frame(RowID=rep(row_id,times=ncol(X)),

 Sample=rep(colnames(X),each=nrow(X)),Log2Abundance=as.vector(X))

donor <- dplyr::inner_join(entries,long,by='RowID') %>%

 dplyr::group_by(ENTREZID,Top20_Rank,Sample) %>%

 dplyr::summarise(Log2Abundance=mean(Log2Abundance),N_ProteinEntries=dplyr::n(),.groups='drop') %>%

 dplyr::left_join(meta[,c('Sample','Donor','Dx','Age','Sex','Batch')],by='Sample') %>%

 dplyr::left_join(top[,c('ENTREZID','GeneSymbols','RankScore')],by='ENTREZID') %>%

 dplyr::arrange(Top20_Rank,Dx,Donor)

W(donor,'08_Top20_Donor_Measured_Log2.csv')

groups <- donor %>% dplyr::group_by(ENTREZID,GeneSymbols,Top20_Rank,Dx) %>%

 dplyr::summarise(n=dplyr::n(),Mean=mean(Log2Abundance),.groups='drop')

W(groups,'09_Top20_Raw_Group_MeansOnly.csv')

centered <- groups %>% dplyr::group_by(ENTREZID) %>%

 dplyr::mutate(CenteredMean=Mean-mean(Mean)) %>% dplyr::ungroup()

W(centered,'10_Top20_Raw_Group_Means_Centered.csv')

# 5. Covariate adjustment: Dx retained; Age, Sex, Batch adjusted.

# Reference = overall mean age and empirical pooled sex/batch distribution.

# This makes estimated marginal means independent of an arbitrary reference sex/batch.

# Per-donor adjusted value = observed - fitted nuisance contribution

#                           + mean fitted nuisance contribution.

# The mean of adjusted donor values within each Dx is the adjusted Dx mean.

# Confidence intervals are from the fitted linear model, not the SD of adjusted dots.

if (!requireNamespace('limma',quietly=TRUE))

  BiocManager::install('limma',ask=FALSE,update=FALSE)

meta <- meta[match(unique(donor$Sample),meta$Sample),,drop=FALSE]

meta$Dx <- factor(meta$Dx,levels=c('HC','MCI','AD'))

meta$Age10 <- (meta$Age-mean(meta$Age))/10

meta$Sex <- factor(meta$Sex)

meta$Batch <- factor(meta$Batch)

# Full 3-group model; include sex/batch only if estimable.

form <- stats::as.formula(paste('~ 0 + Dx + Age10',

  if (nlevels(meta$Sex)>1) '+ Sex' else '',

  if (nlevels(meta$Batch)>1) '+ Batch' else ''))

mm <- stats::model.matrix(form,meta)

if (qr(mm)$rank<ncol(mm))

  stop('Dx, age, sex and batch model is rank-deficient (confounded). Check metadata; no adjusted figure generated.')

if (nrow(mm)<=ncol(mm)) stop('No residual degrees of freedom for adjusted model.')

# Matrix of gene-level donor values in EXACT metadata sample order.

mat <- matrix(NA_real_,nrow=nrow(top),ncol=nrow(meta),

  dimnames=list(top$ENTREZID,meta$Sample))

for (k in seq_len(nrow(donor)))

  mat[donor$ENTREZID[k],donor$Sample[k]] <- donor$Log2Abundance[k]

if (any(!is.finite(mat))) stop('Incomplete gene-by-sample abundance matrix.')

# Each row: gene. Columns: donors. Same design as original GSEA limma.

fit <- limma::lmFit(mat,mm)

# Important: adjusted dots and EMMs are derived from fitted coefficients;

# the original AD-vs-HC ranking is unchanged. CIs use unmoderated OLS

# residual variance (no empirical-Bayes moderation) for interpretability.

# Model matrix evaluated with Dx set to each diagnosis, while the

# observed pooled age/sex/batch composition is kept unchanged.

refX <- lapply(levels(meta$Dx),function(group) {

  m <- meta

  m$Dx <- factor(rep(group,nrow(m)),levels=levels(meta$Dx))

  colMeans(stats::model.matrix(form,m))

})

names(refX) <- levels(meta$Dx)

refX <- do.call(rbind,refX)

if (!identical(colnames(refX),colnames(mm))) stop('Model matrix mismatch.')

# Covariate contribution excludes the 3 diagnosis indicator columns.

dx_cols <- grep('^Dx',colnames(mm))

nuisance_cols <- setdiff(seq_len(ncol(mm)),dx_cols)

B <- fit$coefficients

nuisance <- if (length(nuisance_cols)) {

  B[,nuisance_cols,drop=FALSE] %*% t(mm[,nuisance_cols,drop=FALSE])

} else {

  matrix(0,nrow=nrow(mat),ncol=ncol(mat))

}

mean_nuisance <- rowMeans(nuisance)

adj_mat <- mat - nuisance + mean_nuisance

rownames(adj_mat) <- rownames(mat); colnames(adj_mat) <- colnames(mat)

# Estimated marginal means (EMMs) only; no confidence intervals.

stat_rows <- lapply(seq_len(nrow(top)), function(i) {

  do.call(rbind, lapply(seq_len(nrow(refX)), function(j) {

    estimate <- as.numeric(sum(refX[j, ] * B[i, ]))

    data.frame(ENTREZID=top$ENTREZID[i],

               GeneSymbols=top$GeneSymbols[i],

               Top20_Rank=top$Top20_Rank[i],

               Dx=rownames(refX)[j],

               AdjustedMean=estimate)

  }))

})

adj_groups <- dplyr::bind_rows(stat_rows)

adj_donor <- donor

adj_donor$AdjustedLog2Abundance <- mapply(function(id,sample)

  adj_mat[id,sample],adj_donor$ENTREZID,adj_donor$Sample)

adj_donor$Dx <- factor(adj_donor$Dx,levels=c('HC','MCI','AD'))

adj_groups$Dx <- factor(adj_groups$Dx,levels=c('HC','MCI','AD'))

# Sanity check: mean of adjusted dots equals fitted EMM for each group.

check <- adj_donor %>% dplyr::group_by(ENTREZID,Dx) %>%

  dplyr::summarise(MeanAdjustedDots=mean(AdjustedLog2Abundance),.groups='drop') %>%

  dplyr::left_join(adj_groups[,c('ENTREZID','Dx','AdjustedMean')],

                   by=c('ENTREZID','Dx'))

if (max(abs(check$MeanAdjustedDots-check$AdjustedMean))>1e-7)

  stop('Adjusted donor means do not match model EMMs. Aborting.')

W(adj_donor,'11_Top20_Donor_AgeSexBatch_Adjusted.csv')

W(adj_groups,'12_Top20_EMM_AgeSexBatch_MeansOnly.csv')

W(check,'13_Adjustment_Verification.csv')

# 6. Figure A: adjusted donor dots and RED HORIZONTAL MEAN LINES (NO CI).

lev <- rev(top$GeneSymbols)

adj_donor$GenePlot <- factor(adj_donor$GeneSymbols,levels=lev)

adj_groups$GenePlot <- factor(adj_groups$GeneSymbols,levels=lev)

pA <- ggplot2::ggplot(adj_donor,ggplot2::aes(x=Dx,y=AdjustedLog2Abundance)) +

 ggplot2::geom_point(position=ggplot2::position_jitter(width=.12,height=0),

                     color='#69727E',size=1.2,alpha=.7) +

 ggplot2::geom_segment(data=adj_groups,

   ggplot2::aes(x=as.numeric(Dx)-0.27, xend=as.numeric(Dx)+0.27,

                y=AdjustedMean, yend=AdjustedMean),

   inherit.aes=FALSE,color='#D43F3A',linewidth=1.0) +

 ggplot2::facet_wrap(~GenePlot,scales='free_y',ncol=4) +

 ggplot2::labs(x=NULL,y='Covariate-adjusted log2 TMT abundance',

               title='AD vs HC leading-edge RNA splicing Top 20') +

 ggplot2::theme_classic(base_size=11) +

 ggplot2::theme(strip.background=ggplot2::element_blank(),

               strip.text=ggplot2::element_text(face='bold'),

               plot.title=ggplot2::element_text(face='bold'))

ggplot2::ggsave(file.path(O,'FigureA_Top20_AdjustedDonors_EMM_RedMeanLines.png'),

               pA,width=12.5,height=13,dpi=600,bg='white')

ggplot2::ggsave(file.path(O,'FigureA_Top20_AdjustedDonors_EMM_RedMeanLines.pdf'),

               pA,width=12.5,height=13,bg='white')

# 7. Figure B: covariate-adjusted group EMMs, centered per gene.

adj_centered <- adj_groups %>% dplyr::group_by(ENTREZID) %>%

 dplyr::mutate(CenteredAdjustedMean=AdjustedMean-mean(AdjustedMean)) %>%

 dplyr::ungroup()

adj_centered$GenePlot <- factor(adj_centered$GeneSymbols,levels=lev)

W(adj_centered,'14_Top20_Adjusted_EMM_Centered.csv')

pB <- ggplot2::ggplot(adj_centered,ggplot2::aes(x=Dx,y=GenePlot)) +

 ggplot2::geom_point(ggplot2::aes(fill=CenteredAdjustedMean),shape=21,

                     size=7,color='white',stroke=.4) +

 ggplot2::scale_fill_gradient2(low='#3266A8',mid='#F4F4F4',high='#B84945',

                               midpoint=0,

                               name='Centered adjusted mean\n(log2 abundance)') +

 ggplot2::labs(x=NULL,y=NULL,

               title='AD vs HC leading-edge Top 20 | adjusted group means') +

 ggplot2::theme_classic(base_size=12) +

 ggplot2::theme(axis.line=ggplot2::element_blank(),

               axis.ticks=ggplot2::element_blank(),

               axis.text.x=ggplot2::element_text(face='bold'),

               plot.title=ggplot2::element_text(face='bold'))

ggplot2::ggsave(file.path(O,'FigureB_Top20_AdjustedGroupMean_DotMatrix_MeanOnly.png'),

               pB,width=8.5,height=9,dpi=600,bg='white')

ggplot2::ggsave(file.path(O,'FigureB_Top20_AdjustedGroupMean_DotMatrix_MeanOnly.pdf'),

               pB,width=8.5,height=9,bg='white')

writeLines(c(

 'Top20 selection: AD vs HC leading-edge union of seven RNA-splicing GO terms;',

 '20 most negative AD vs HC RankScore. MCI NOT used in selection.',

 'For gene-level measured abundance: mean log2 protein entries per donor.',

 paste('Adjustment model:',deparse(form)),

 'Covariate-adjusted donor dots remove Age/Sex/Batch fitted effects only;',

 'diagnosis-group effects are preserved. Age centered at pooled sample mean.',

 'EMMs average model predictions over the pooled observed Age/Sex/Batch distribution.',

 'Figure A shows adjusted donor dots and RED HORIZONTAL EMM MEAN LINES, without confidence intervals.',

 'Figure B is centered across the 3 adjusted group EMMs within each gene.',

 'No significance tests were added. Adjusted plots are not independent validation.',

 'Compare with raw donor/mean CSVs 08-10 for sensitivity.',

 paste('Proteome:',P),paste('Mapping:',C),paste('Ranking:',R),paste('GSEA:',G)),

 file.path(O,'Analysis_Settings_Adjusted_RedMeanLines.txt'))

capture.output(sessionInfo(),file=file.path(O,'SessionInfo_Adjusted_RedMeanLines.txt'))

print(top)

print(pA);print(pB)

message('DONE: ',O)
