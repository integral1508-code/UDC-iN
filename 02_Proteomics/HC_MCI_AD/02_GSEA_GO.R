# Donor exclusion: The 72-year-old donor originally classified as AD was excluded from the analyses because the updated clinical diagnosis indicated that the donor did not have AD. Both matched UDC and UDC-iN samples from this donor were excluded, leaving 11 AD UDC-iN donors.

# GSEA REANALYSIS: nPermSimple=10000; original AD/MCI design and global BH preserved.

# ================================================================

# UDC-iN | HC, MCI, AD | GO BP preranked GSEA + shared GO analysis

# Contrasts: MCI-HC, AD-HC, AD-MCI

# Model: Diagnosis + Age + Sex + TMT batch (when estimable)

# Ranking: limma moderated t; MEAN t for all protein entries per ENTREZID

# No imputation; complete-case across selected UDC-iN samples

# GO plots: UP x=NES (red); DOWN x=-NES (blue)

# Dot size=setSize, dot color=global BH-adjusted P (3 contrasts)

# Common significant GO IDs: pairwise and all 3; UP/DOWN separately

# ================================================================

options(stringsAsFactors=FALSE, repos=c(CRAN="https://cloud.r-project.org"))

cran <- c("readxl","readr","dplyr","ggplot2","stringr","scales")

bioc <- c("limma","clusterProfiler","org.Hs.eg.db","AnnotationDbi")

for (p in cran) if (!requireNamespace(p,quietly=TRUE))

  install.packages(p,repos="https://cloud.r-project.org")

if (!requireNamespace("BiocManager",quietly=TRUE))

  install.packages("BiocManager",repos="https://cloud.r-project.org")

for (p in bioc) if (!requireNamespace(p,quietly=TRUE))

  BiocManager::install(p,ask=FALSE,update=FALSE)

suppressPackageStartupMessages({

  library(readxl);library(readr);library(dplyr);library(ggplot2)

  library(limma);library(clusterProfiler);library(org.Hs.eg.db)

})

set.seed(20261016)

D <- path.expand(Sys.getenv("PROTEOME_DATA_DIR", unset = "~/Desktop"))

if (!dir.exists(D)) D <- path.expand("~")

choose_input <- function(label,pattern) {

  f <- list.files(D,pattern="\\.xlsx$",full.names=TRUE,ignore.case=TRUE)

  hits <- f[grepl(pattern,basename(f),ignore.case=TRUE)]

  if(length(hits)==1L) return(hits)

  if(length(hits)>1L) {

    # Prefer most recently modified if multiple versions exist.

    hits <- hits[order(file.info(hits)$mtime,decreasing=TRUE)]

    message("Multiple ",label," files; using newest: ",basename(hits[1]))

    return(hits[1])

  }

  if(interactive()) return(file.choose())

  stop("Cannot find ",label," in ",D)

}

P <- file.path(D, "20230316_iN_batch_All_ys_Proteins_KM edited.xlsx")

C <- file.path(D, "Proteome_対応表(2).xlsx")

O <- file.path(D,"UDCiN_GO_BP_GSEA_MeanT_WithBatch_GlobalBH3Contrasts_nPerm10000_AD11")

dir.create(O,showWarnings=FALSE,recursive=TRUE)

W <- function(z,name) readr::write_csv(as.data.frame(z),file.path(O,name),na="")

clean <- function(x) trimws(as.character(x))

message("Proteome: ",P,"\nMapping: ",C,"\nOutput: ",O)

# ---- 1. Metadata and diagnosis correspondence -----------------------------

cl <- as.data.frame(readxl::read_excel(C,col_names=FALSE,.name_repair="minimal"))

if(ncol(cl)<2) stop("Correspondence needs >=2 columns")

cl <- data.frame(Sample=clean(cl[[1]]),Label=clean(cl[[2]]))

cl$Dx <- dplyr::case_when(

  cl$Label=="1.HC" ~ "HC",

  cl$Label %in% c("1.MCIduetoAD","2.MCIduetoAD") ~ "MCI",

  cl$Label %in% c("2.AD","3.AD") ~ "AD",

  TRUE ~ NA_character_

)

cl <- cl[!is.na(cl$Dx)&!is.na(cl$Sample)&nzchar(cl$Sample),]

if(anyDuplicated(cl$Sample)) stop("Duplicate sample in correspondence")

md <- as.data.frame(readxl::read_excel(P,sheet="metadata",.name_repair="minimal"))

needed <- c("File Name","Type","ID2","Sex","Age","batch")

if(!all(needed %in% names(md))) stop("Missing metadata: ",paste(setdiff(needed,names(md)),collapse=", "))

meta <- data.frame(

  Sample=clean(md[["File Name"]]),CellType=clean(md[["Type"]]),

  Donor=clean(md[["ID2"]]),Sex=clean(md[["Sex"]]),

  Batch=clean(md[["batch"]]),

  Age=suppressWarnings(as.numeric(md[["Age"]]))

)

if(anyDuplicated(meta$Sample)) stop("Duplicate sample in metadata")

meta <- merge(meta,cl,by="Sample",sort=FALSE)

# Exclude the two samples from the same 72-year-old AD donor.

# This analysis remains UDC-iN-only; UDC exclusion is recorded for audit.

exclude_names <- c(

  UDCiN="Abundances (Scaled): F3: 131C, Sample, 03_AD__iN, 2_iN, 3_AD",

  UDC="Abundances (Scaled): F3: 127C, Sample, 3_AD_UDC, 1_UDC, 3_AD"

)

excluded <- meta[match(unname(exclude_names), meta$Sample), , drop=FALSE]

if(nrow(excluded)!=2L || anyNA(excluded$Sample) ||

   !identical(as.character(excluded$Sample), unname(exclude_names)))

  stop("Cannot uniquely find both exact AD sample names in merged metadata.")

if(any(excluded$Dx != "AD") ||

   any(!is.finite(excluded$Age) | excluded$Age != 72) ||

   length(unique(excluded$Donor)) != 1L ||

   !identical(as.character(excluded$CellType), c("iN","UDC")))

  stop("Excluded samples are not the expected same 72-year-old AD donor (iN + UDC).")

W(excluded,"01_Excluded_AD_72y_UDC_and_UDCiN.csv")

meta <- meta[!(meta$Sample %in% unname(exclude_names)), , drop=FALSE]

meta <- meta[meta$CellType=="iN" & is.finite(meta$Age) &

               !is.na(meta$Sex)&nzchar(meta$Sex)&

               !is.na(meta$Donor)&nzchar(meta$Donor)&

               !is.na(meta$Batch)&nzchar(meta$Batch),,drop=FALSE]

meta$Dx <- factor(meta$Dx,levels=c("HC","MCI","AD"))

meta$Sex <- droplevels(factor(meta$Sex))

meta$Batch <- droplevels(factor(meta$Batch))

if(anyDuplicated(meta$Sample)) stop("Duplicate selected sample")

if(anyDuplicated(meta$Donor))

  stop("Repeated UDC-iN donors: this code assumes one iN sample per donor.")

if(any(table(meta$Dx)<2)) stop("Need at least 2 donors in each diagnosis")

if(sum(meta$Dx=="AD") != 11L) stop("Expected 11 AD UDC-iN donors after exclusion; found ", sum(meta$Dx=="AD"))

W(meta,"01_Selected_iN_Metadata.csv")

W(as.data.frame(table(meta$Dx)),"01_Diagnosis_Counts.csv")

print(table(meta$Dx));print(table(meta$Dx,meta$Sex))

print(table(meta$Dx,meta$Batch))

W(as.data.frame(table(meta$Dx,meta$Batch)),"01_Diagnosis_By_Batch.csv")

# ---- 2. TMT matrix --------------------------------------------------------

sheets <- readxl::excel_sheets(P)

dat_sheet <- sheets[which(tolower(sheets)!="metadata")[1]]

raw <- as.data.frame(readxl::read_excel(P,sheet=dat_sheet,.name_repair="minimal"),check.names=FALSE)

if(!all(c("Accession","Gene Symbol") %in% names(raw))) stop("Protein annotation columns missing")

if(!all(meta$Sample %in% names(raw)))

  stop("Missing sample columns: ",paste(setdiff(meta$Sample,names(raw)),collapse=", "))

X0 <- as.matrix(raw[,meta$Sample,drop=FALSE])

suppressWarnings(storage.mode(X0)<-"numeric")

rownames(X0) <- as.character(seq_len(nrow(raw)))

keep <- rowSums(is.finite(X0)&X0>0)==ncol(X0)

X <- log2(X0[keep,,drop=FALSE])

sdv <- apply(X,1,sd)

X <- X[is.finite(sdv)&sdv>0,,drop=FALSE]

if(nrow(X)<100) stop("Too few complete-case proteins")

anno <- data.frame(

  RowID=as.integer(rownames(X)),

  Accession=clean(raw[["Accession"]][as.integer(rownames(X))]),

  GeneSymbol=clean(raw[["Gene Symbol"]][as.integer(rownames(X))])

)

W(data.frame(Metric=c("iN_samples","total_protein_rows","complete_case_proteins"),

             Value=c(ncol(X),nrow(raw),nrow(X))),"02_QC.csv")

# ---- 3. limma age + sex + batch (when estimable) --------------------------

meta$Age10 <- (meta$Age-mean(meta$Age))/10

covars <- c("Dx","Age10",

            if(nlevels(meta$Sex)>1) "Sex",

            if(nlevels(meta$Batch)>1) "Batch")

formula <- as.formula(paste("~ 0 +",paste(covars,collapse=" + ")))

design <- model.matrix(formula,data=meta)

if(qr(design)$rank<ncol(design)) {

  bad <- colnames(design)[qr(design)$pivot[seq.int(qr(design)$rank+1,ncol(design))]]

  stop("Rank-deficient design; aliased columns: ",paste(bad,collapse=", "),

       ". Check diagnosis versus batch. Batch is NOT silently removed.")

}

message("Model: ",deparse(formula))

if(nrow(design)-ncol(design)<3)

  stop("Insufficient residual degrees of freedom for age/sex-adjusted model")

W(data.frame(Sample=meta$Sample,design,check.names=FALSE),"03_Design_Matrix.csv")

ct <- limma::makeContrasts(

  MCI_vs_HC=DxMCI-DxHC,

  AD_vs_HC=DxAD-DxHC,

  AD_vs_MCI=DxAD-DxMCI,

  levels=design

)

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(X,design),ct),trend=TRUE)

# ---- 4. GSEA, mean moderated t per gene ---------------------------------

# Map SYMBOL -> ENTREZID first, then calculate mean t per ENTREZID.

# All valid mapped protein entries contribute equally to the mean.

# One rank score per ENTREZID; no DEP significance pre-filter.

gsea_results <- list()

summary_rows <- list()

go_fdr <- 0.05

top_n <- 10

save_plot <- function(g,stem,height=6) {

  ggplot2::ggsave(file.path(O,paste0(stem,".png")),g,

                  width=9.5,height=height,dpi=600,bg="white")

  ggplot2::ggsave(file.path(O,paste0(stem,".pdf")),g,

                  width=9.5,height=height,bg="white")

}

all_gsea_tabs <- list()

ranked_gene_counts <- list()

for (cn in colnames(ct)) {

  tt <- limma::topTable(fit,coef=cn,number=Inf,sort.by="none",adjust.method="BH")

  tt$RowID <- as.integer(rownames(tt))

  tt <- dplyr::left_join(tt,anno,by="RowID")

  W(tt,paste0("04_limma_",cn,"_ALL.csv"))

  gs <- tt %>%

    dplyr::filter(!is.na(GeneSymbol),nzchar(GeneSymbol),

                  !grepl("[,;/|]",GeneSymbol),is.finite(t))

  mapped <- suppressMessages(clusterProfiler::bitr(

    unique(gs$GeneSymbol),fromType="SYMBOL",toType="ENTREZID",

    OrgDb=org.Hs.eg.db::org.Hs.eg.db))

  if(is.null(mapped)||nrow(mapped)<100) stop("Gene mapping failed: ",cn)

  mapped_proteins <- dplyr::inner_join(gs,mapped,by=c("GeneSymbol"="SYMBOL")) %>%

    dplyr::distinct(ENTREZID,RowID,.keep_all=TRUE)

  # Average the moderated t-statistics for every protein row mapped to a gene.

  # This is a ranking score, not a newly estimated gene-level t-test.

  ranked <- mapped_proteins %>%

    dplyr::group_by(ENTREZID) %>%

    dplyr::summarise(

      RankScore=mean(t),

      ProteinRows=dplyr::n(),

      GeneSymbols=paste(sort(unique(GeneSymbol)),collapse=";"),

      .groups="drop"

    ) %>%

    dplyr::arrange(dplyr::desc(RankScore),ENTREZID)

  ranks <- ranked$RankScore

  names(ranks) <- ranked$ENTREZID

  if(anyDuplicated(names(ranks))) stop("Nonunique ENTREZ IDs")

  ranks <- sort(ranks,decreasing=TRUE)

  W(mapped_proteins,paste0("05a_Mapped_Protein_Entries_",cn,".csv"))

  W(ranked,paste0("05_Ranked_Genes_MeanT_",cn,".csv"))

  set.seed(20261016)

  go <- clusterProfiler::gseGO(

    geneList=ranks,OrgDb=org.Hs.eg.db::org.Hs.eg.db,

    keyType="ENTREZID",ont="BP",minGSSize=10,maxGSSize=500,

    pvalueCutoff=1,pAdjustMethod="BH",eps=0,verbose=FALSE,

    by="fgsea",seed=TRUE,nPermSimple=10000

  )

  tab <- as.data.frame(go)

  if(nrow(tab)) {

    tab$Contrast <- cn

    tab$Direction <- ifelse(tab$NES>0,"UP","DOWN")

  }

  all_gsea_tabs[[cn]] <- tab

  ranked_gene_counts[[cn]] <- length(ranks)

}

combined_gsea <- dplyr::bind_rows(all_gsea_tabs)

if (!nrow(combined_gsea)) stop("No GSEA terms returned across contrasts")

combined_gsea$p.adjust_within_contrast <- combined_gsea$p.adjust

valid_p <- is.finite(combined_gsea$pvalue)

combined_gsea$p.adjust_global_3contrasts <- NA_real_

combined_gsea$p.adjust_global_3contrasts[valid_p] <-

  stats::p.adjust(combined_gsea$pvalue[valid_p], method="BH")

combined_gsea$p.adjust <- combined_gsea$p.adjust_global_3contrasts

W(combined_gsea,"06_GSEA_GO_BP_ALL_3Contrasts_GlobalBH.csv")

W(data.frame(Contrast=names(all_gsea_tabs),

             TestedTerms=vapply(all_gsea_tabs,nrow,integer(1))),

  "06_GSEA_Number_Of_Tests_By_Contrast.csv")

for (cn in colnames(ct)) {

  tab <- combined_gsea[combined_gsea$Contrast==cn,,drop=FALSE]

  W(tab,paste0("06_GSEA_GO_BP_",cn,"_ALL.csv"))

  sig <- if(nrow(tab)) tab[is.finite(tab$p.adjust)&tab$p.adjust<go_fdr,,drop=FALSE] else tab

  W(sig,paste0("07_GSEA_GO_BP_",cn,"_FDR005.csv"))

  gsea_results[[cn]] <- sig

  summary_rows[[cn]] <- data.frame(

    Contrast=cn,RankedGenes=ranked_gene_counts[[cn]],TestedTerms=nrow(tab),

    Significant=nrow(sig),UP=sum(sig$NES>0),DOWN=sum(sig$NES<0))

  for(dir in c("UP","DOWN")) {

    z <- sig[sig$Direction==dir,,drop=FALSE]

    z <- z[order(z$p.adjust,-abs(z$NES),z$ID),,drop=FALSE]

    z <- head(z,top_n)

    W(z,paste0("08_GO_BP_",cn,"_",dir,"_Top10.csv"))

    if(!nrow(z)) {

      message(cn," ",dir,": no FDR<0.05 terms")

      next

    }

    z$PlotNES <- if(dir=="UP") z$NES else -z$NES

    z$Term <- stringr::str_wrap(z$Description,width=43)

    z <- z[order(-z$PlotNES,z$p.adjust,z$ID),,drop=FALSE]

    z$TermID <- factor(z$ID,levels=rev(z$ID))

    labs_term <- stats::setNames(z$Term,z$ID)

    cols <- if(dir=="UP") c("#B9473D","#F0C6B8") else c("#235C91","#BFD8ED")

    g <- ggplot(z,aes(x=PlotNES,y=TermID))+

      geom_point(aes(size=setSize,color=p.adjust),alpha=.95)+

      scale_y_discrete(labels=labs_term)+

      scale_size_continuous(range=c(3.5,9),name="Gene-set size")+

      scale_color_gradient(low=cols[1],high=cols[2],

        limits=c(0,go_fdr),oob=scales::squish,name="Global BH FDR")+

      scale_x_continuous(expand=expansion(mult=c(.04,.12)))+

      labs(title=paste(cn,"| GO BP |",dir),

           x=if(dir=="UP") "NES" else expression(-NES),y=NULL)+

      theme_classic(base_size=12)+

      theme(plot.title=element_text(face="bold",size=14),

            axis.text.y=element_text(size=9),

            axis.line.y=element_blank(),axis.ticks.y=element_blank())

    save_plot(g,paste0("09_GO_BP_",cn,"_",dir,"_NES_Dot"),

              height=max(5.5,2.7+.36*nrow(z)))

  }

}

W(dplyr::bind_rows(summary_rows),"10_GSEA_Summary.csv")

# ---- 5. Significant GO overlaps: exact GO ID and same direction ----------

# A GO is 'shared' only if it has BH FDR<0.05 in EACH selected comparison,

# AND its NES has the SAME SIGN in all comparisons.

cnames <- colnames(ct)

pairs <- combn(cnames,2,simplify=FALSE)

overlap_summary <- list()

overlap_tables <- list()

for(dir in c("UP","DOWN")) {

  sets <- lapply(gsea_results,function(z) {

    if(!nrow(z)) return(character(0))

    unique(as.character(z$ID[z$Direction==dir]))

  })

  triple <- Reduce(intersect,sets)

  groupings <- c(pairs,list(cnames))

  for(group in groupings) {

    ids <- Reduce(intersect,sets[group])

    ids <- sort(ids)

    label <- paste(group,collapse="__AND__")

    result <- data.frame(ID=ids)

    for(cn in group) {

      z <- gsea_results[[cn]]

      if(nrow(z)) {

        zz <- z[match(ids,z$ID),c("ID","Description","NES","p.adjust","setSize"),drop=FALSE]

        names(zz)[-1] <- paste0(cn,"_",names(zz)[-1])

        result <- dplyr::left_join(result,zz,by="ID")

      }

    }

    # retain one human-readable description

    if(nrow(result)) {

      description_col <- paste0(group[1],"_Description")

      result <- result[,c("ID",description_col,setdiff(names(result),c("ID",description_col))),drop=FALSE]

      names(result)[2] <- "Description"

    }

    W(result,paste0("11_Shared_",dir,"_",label,".csv"))

    overlap_tables[[paste(dir,label,sep="_")]] <- result

    overlap_summary[[paste(dir,label,sep="_")]] <-

      data.frame(Direction=dir,Comparisons=paste(group,collapse=" & "),

                 NumberOfComparisons=length(group),SharedGO=nrow(result))

  }

  # GO significant in exactly two or three comparisons (not necessarily pair-specific)

  union_ids <- sort(unique(unlist(sets,use.names=FALSE)))

  freq <- vapply(union_ids,function(id)sum(vapply(sets,function(z)id %in% z,logical(1))),integer(1))

  frequency <- data.frame(ID=union_ids,N_Comparisons=freq,

                          MCI_vs_HC=union_ids %in% sets[["MCI_vs_HC"]],

                          AD_vs_HC=union_ids %in% sets[["AD_vs_HC"]],

                          AD_vs_MCI=union_ids %in% sets[["AD_vs_MCI"]])

  if(nrow(frequency)) {

    for(cn in cnames) {

      zz <- gsea_results[[cn]]

      frequency[[paste0("NES_",cn)]] <- zz$NES[match(frequency$ID,zz$ID)]

      frequency[[paste0("FDR_",cn)]] <- zz$p.adjust[match(frequency$ID,zz$ID)]

    }

    allterms <- dplyr::bind_rows(gsea_results)

    frequency$Description <- allterms$Description[match(frequency$ID,allterms$ID)]

    frequency <- frequency[,c("ID","Description",setdiff(names(frequency),c("ID","Description")))]

  }

  W(frequency,paste0("12_GO_",dir,"_Presence_All_Comparisons.csv"))

  W(frequency[frequency$N_Comparisons>=2,,drop=FALSE],

    paste0("13_GO_",dir,"_Shared_AtLeast2.csv"))

  W(frequency[frequency$N_Comparisons==3,,drop=FALSE],

    paste0("14_GO_",dir,"_Shared_All3.csv"))

}

W(dplyr::bind_rows(overlap_summary),"15_Shared_GO_Counts.csv")

# ---- 6. Documentation -----------------------------------------------------

writeLines(c(

 "UDC-iN only; HC/MCI/AD three contrasts: MCI-HC, AD-HC, AD-MCI",

 "Excluded one 72y AD donor (both UDC and UDC-iN sample names); AD iN donors=11.",

 "Protein matrix: TMT scaled abundance -> log2, complete cases across all selected iN samples.",

 paste("Model:",deparse(formula),". TMT batch included when >1 batch; no silent removal of confounded batch terms."),

 "One sample per donor; no duplicateCorrelation.",

 "limma eBayes(trend=TRUE), moderated t ranking, no DEP significance prefilter.",

 "Gene ID mapping SYMBOL->ENTREZID; gene ranking score = arithmetic mean of protein-entry moderated t per ENTREZID. This score is NOT a gene-level t-test.",

 "GO BP gseGO by fgsea with nPermSimple=10000; set size 10-500. All GO p-values from all 3 contrasts pooled for one global BH FDR <0.05; within-contrast BH retained as p.adjust_within_contrast.",

 "GO UP NES>0; GO DOWN NES<0; DOWN plot uses -NES only for x-axis.",

 "Dot color adjusted P (red UP, blue DOWN); dot size = gene-set size.",

 "Overlap uses EXACT GO ID, global 3-contrast BH FDR<0.05, SAME NES DIRECTION, across comparisons.",

 "For comparisons, sign direction: MCI-HC, AD-HC, AD-MCI.",

 paste("Proteome file:",P),paste("Mapping file:",C)),

 file.path(O,"Analysis_Settings.txt"))

capture.output(sessionInfo(),file=file.path(O,"SessionInfo.txt"))

message("DONE: ",O)

print(dplyr::bind_rows(summary_rows))

print(dplyr::bind_rows(overlap_summary))
