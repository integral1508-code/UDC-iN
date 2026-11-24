# Donor exclusion: The 72-year-old donor originally classified as AD was excluded from the analyses because the updated clinical diagnosis indicated that the donor did not have AD. Both matched UDC and UDC-iN samples from this donor were excluded, leaving 11 AD UDC-iN donors.

# ============================================================================

# HC / MCI / AD TMT proteomics: Volcano + GO BP UP/DOWN (UDC-iN ONLY)

# Contrasts: MCI-HC, AD-HC, AD-MCI (positive = first-named group higher)

# limma empirical Bayes, independent UDC-iN donors, age+sex+batch adjusted

# No imputation; complete-case over all selected UDC-iN HC/MCI/AD samples

# BH FDR < .05, no fold-change cutoff; GO top 10 by FDR, shown GeneRatio desc

# ============================================================================

options(stringsAsFactors=FALSE)

cran <- c("readxl","ggplot2","ggrepel","stringr","statmod")

bioc <- c("limma","clusterProfiler","org.Hs.eg.db","AnnotationDbi")

for (pkg in cran) if (!requireNamespace(pkg,quietly=TRUE))

  install.packages(pkg,repos="https://cloud.r-project.org")

if (!requireNamespace("BiocManager",quietly=TRUE))

  install.packages("BiocManager",repos="https://cloud.r-project.org")

for (pkg in bioc) if (!requireNamespace(pkg,quietly=TRUE))

  BiocManager::install(pkg,ask=FALSE,update=FALSE)

suppressPackageStartupMessages({

 library(readxl);library(ggplot2);library(ggrepel)

 library(limma);library(clusterProfiler);library(org.Hs.eg.db)

 library(AnnotationDbi)

})

set.seed(123)

desktop <- path.expand(Sys.getenv("PROTEOME_DATA_DIR", unset = "~/Desktop"))

find_excel <- function(exact,regex) {

  f <- file.path(desktop,exact)

  if(file.exists(f)) return(f)

  matches <- list.files(desktop,pattern=regex,full.names=TRUE)

  if(length(matches)!=1L) stop("Input ambiguous/missing: ",exact,

                               "\nCandidates: ",paste(basename(matches),collapse=", "))

  matches

}

P <- file.path(desktop, "20230316_iN_batch_All_ys_Proteins_KM edited.xlsx")
C <- file.path(desktop, "Proteome_対応表(2).xlsx")
stopifnot(file.exists(P), file.exists(C))

out <- file.path(desktop,"Proteome_UDCiN_HC_MCI_AD_11AD_AgeSexBatchAdjusted")

dir.create(out,recursive=TRUE,showWarnings=FALSE)

csv <- function(x,f) write.csv(as.data.frame(x),f,row.names=FALSE,na="")

save_fig <- function(g,stem,w=9,h=6.5) {

  ggsave(paste0(stem,".png"),g,width=w,height=h,dpi=600,bg="white")

  ggsave(paste0(stem,".pdf"),g,width=w,height=h,device="pdf")

}

red <- "#CC6859"; blue <- "#3C78A6"; gray <- "#C9CDD1"

fdr <- .05; go_fdr <- .05; top_n <- 10

trim <- function(z) trimws(as.character(z))

# --------------------------------------------------------------------------

# 1. Diagnosis correspondence: classify by suffix, not numeric prefix.

# The mapping contains both 1.MCIduetoAD/2.MCIduetoAD and 2.AD/3.AD.

# Avoid accidentally excluding those donors due to inconsistent numbering.

# --------------------------------------------------------------------------

mp <- as.data.frame(read_excel(C,col_names=FALSE,.name_repair="minimal"))

stopifnot(ncol(mp)>=3)

map <- data.frame(Sample=trim(mp[[1]]),DxRaw=trim(mp[[2]]),

                  Age=suppressWarnings(as.numeric(mp[[3]])))

map$Dx <- ifelse(grepl("^[0-9]+\\.HC$",map$DxRaw),"HC",

          ifelse(grepl("^[0-9]+\\.MCIduetoAD$",map$DxRaw),"MCI",

          ifelse(grepl("^[0-9]+\\.AD$",map$DxRaw),"AD",NA_character_)))

map <- map[!is.na(map$Dx)&is.finite(map$Age)&!is.na(map$Sample),]

if(anyDuplicated(map$Sample)) stop("Duplicated mapping sample names.")

md <- as.data.frame(read_excel(P,sheet="metadata",.name_repair="minimal"))

need <- c("File Name","Type","ID2","Sex","batch")

if(!all(need %in% names(md))) stop("Missing metadata: ",

                                    paste(setdiff(need,names(md)),collapse=", "))

meta <- data.frame(Sample=trim(md[["File Name"]]),

                   Cell=trim(md[["Type"]]),

                   Donor=trim(md[["ID2"]]),Sex=trim(md[["Sex"]]),

                   Batch=trim(md[["batch"]]))

# Identify BOTH supplied abundance columns in the original metadata before

# restricting the analysis to UDC-iNs.  The analysis itself remains iN-only.

exclude_in <- "Abundances (Scaled): F3: 131C, Sample, 03_AD__iN, 2_iN, 3_AD"

exclude_udc <- "Abundances (Scaled): F3: 127C, Sample, 3_AD_UDC, 1_UDC, 3_AD"

requested_exclusions <- c(exclude_in, exclude_udc)

if(anyDuplicated(meta$Sample)) stop("Duplicate File Name values in metadata.")

excl_idx <- match(requested_exclusions, meta$Sample)

if(anyNA(excl_idx)) stop("Cannot find exact excluded abundance headers in metadata: ",

                         paste(requested_exclusions[is.na(excl_idx)],collapse=" | "))

excluded_both <- meta[excl_idx,,drop=FALSE]

if(length(unique(excluded_both$Donor)) != 1L ||

   anyNA(excluded_both$Donor) || !nzchar(excluded_both$Donor[1]))

  stop("Excluded UDC and UDC-iN do not map to the same donor ID.")

if(!setequal(excluded_both$Cell,c("UDC","iN")))

  stop("Expected exactly one UDC and one iN in the two excluded samples.")

excluded_donor <- excluded_both$Donor[1]

meta <- merge(meta,map,by="Sample",sort=FALSE)

# Verify the UDC-iN sample is AD and the mapped donor age is 72 years.

excl_in_map <- meta[meta$Sample==exclude_in,,drop=FALSE]

if(nrow(excl_in_map)!=1L || excl_in_map$Dx!="AD" ||

   excl_in_map$Donor!=excluded_donor ||

   !is.finite(excl_in_map$Age) || excl_in_map$Age!=72)

  stop("Excluded UDC-iN must map to AD, age 72, and the same donor.")

# If the UDC sample has a diagnosis/age mapping, validate it as well.

excl_udc_map <- meta[meta$Sample==exclude_udc,,drop=FALSE]

if(nrow(excl_udc_map)>1L ||

   (nrow(excl_udc_map)==1L &&

    (excl_udc_map$Dx!="AD" || excl_udc_map$Donor!=excluded_donor ||

     !is.finite(excl_udc_map$Age) || excl_udc_map$Age!=72)))

  stop("Excluded UDC metadata conflicts with AD, age 72, same donor.")

# Explicitly exclude both matched samples from the metadata before analysis.

meta <- meta[!(meta$Sample %in% requested_exclusions),,drop=FALSE]

if(any(meta$Donor==excluded_donor))

  stop("Other samples from the excluded donor remain in the input metadata.")

meta <- meta[meta$Cell == "iN" &

               !is.na(meta$Donor)&nzchar(meta$Donor)&

               !is.na(meta$Sex)&nzchar(meta$Sex)&

               !is.na(meta$Batch)&nzchar(meta$Batch),]

if(anyDuplicated(meta$Sample)) stop("Duplicated sample names.")

if(anyDuplicated(meta$Donor)) stop("Multiple UDC-iN samples per donor: independent-donor model requires one sample per donor. Inspect metadata.")

if(sum(meta$Dx=="AD")!=11L)

  stop("Expected 11 AD UDC-iN donors after excluding the matched 72y donor; found ",sum(meta$Dx=="AD"))

excluded_both$Dx <- "AD"

excluded_both$Age <- 72

csv(excluded_both,file.path(out,"Excluded_AD_UDC_and_UDCiN.csv"))

message("Excluded 72y AD donor: ",excluded_donor,

        " (both UDC and UDC-iN; 11 AD UDC-iN donors remain)")

meta <- meta[order(meta$Donor),]

meta$Dx <- factor(meta$Dx,levels=c("HC","MCI","AD"))

meta$Sex <- droplevels(factor(meta$Sex))

meta$Batch <- droplevels(factor(meta$Batch))

meta$Age10 <- (meta$Age-mean(meta$Age))/10

donor_counts <- table(meta$Dx[!duplicated(meta$Donor)])

print(donor_counts)

if(any(donor_counts==0)) stop("Missing at least one diagnosis group.")

csv(meta,file.path(out,"Included_Samples.csv"))

csv(data.frame(Diagnosis=names(donor_counts),Donors=as.integer(donor_counts)),

    file.path(out,"Group_Counts.csv"))

message("NOTE: diagnosis groups are recognized by text suffix; inspect Included_Samples.csv ",

        "to confirm the updated diagnosis classification.")

# --------------------------------------------------------------------------

# 2. Complete-case proteins across ALL selected samples, no imputation.

# --------------------------------------------------------------------------

sheet <- setdiff(excel_sheets(P),"metadata")[1]

prot <- as.data.frame(read_excel(P,sheet=sheet,.name_repair="minimal"))

if(!all(c("Accession","Gene Symbol") %in% names(prot)))

  stop("Missing protein annotations.")

if(!all(meta$Sample %in% names(prot)))

  stop("Missing abundance columns: ",paste(setdiff(meta$Sample,names(prot)),collapse="; "))

raw <- as.matrix(data.frame(lapply(prot[,meta$Sample,drop=FALSE],

                                    function(z)suppressWarnings(as.numeric(as.character(z)))),

                            check.names=FALSE))

colnames(raw) <- meta$Sample

rownames(raw) <- as.character(seq_len(nrow(raw)))

keep <- rowSums(is.finite(raw)&raw>0)==ncol(raw)

x <- log2(raw[keep,,drop=FALSE])

v <- apply(x,1,var)

x <- x[is.finite(v)&v>0,,drop=FALSE]

if(nrow(x)<100) stop("Too few complete-case proteins; review missingness.")

message("Complete-case UDC-iN proteins: ",nrow(x)," / ",nrow(raw))

anno <- data.frame(ProteinRowID=rownames(x),

  Accession=trim(prot[["Accession"]][as.integer(rownames(x))]),

  SYMBOL=trim(prot[["Gene Symbol"]][as.integer(rownames(x))]))

anno$SYMBOL[is.na(anno$SYMBOL)|!nzchar(anno$SYMBOL)|

               grepl("[,;/|]",anno$SYMBOL)] <- NA_character_

syms <- unique(na.omit(anno$SYMBOL))

sym_map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,

  keys=syms,keytype="SYMBOL",column="ENTREZID",multiVals="first")

anno$ENTREZID <- unname(as.character(sym_map[anno$SYMBOL]))

# UniProt fallback for missing SYMBOL mapping

missing <- is.na(anno$ENTREZID)|!nzchar(anno$ENTREZID)

acc <- unique(anno$Accession[missing])

acc <- intersect(acc,AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype="UNIPROT"))

if(length(acc)) {

  amap <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=acc,

       keytype="UNIPROT",column="ENTREZID",multiVals="first")

  anno$ENTREZID[missing] <- unname(as.character(amap[anno$Accession[missing]]))

}

csv(anno,file.path(out,"Protein_Annotation.csv"))

# --------------------------------------------------------------------------

# 3. UDC-iN only, independent donors, age + sex + batch adjusted.

# Contrasts: MCI-HC, AD-HC, AD-MCI.

# --------------------------------------------------------------------------

design <- model.matrix(~0+Dx+Age10+Sex+Batch,data=meta)

colnames(design) <- sub("^Dx","",colnames(design))

if(qr(design)$rank!=ncol(design)) stop("Design rank deficient: check Age/Sex/Dx/Batch confounding.")

if(nrow(design)-ncol(design)<3) stop("Insufficient residual df.")

fit <- lmFit(x,design)

contrast_text <- c(

 MCI_vs_HC="MCI-HC",

 AD_vs_HC="AD-HC",

 AD_vs_MCI="AD-MCI")

cm <- makeContrasts(contrasts=unname(contrast_text),levels=design)

colnames(cm) <- names(contrast_text)

fit <- eBayes(contrasts.fit(fit,cm),trend=TRUE,robust=FALSE)

csv(data.frame(Sample=meta$Sample,design,check.names=FALSE),

    file.path(out,"Design_Matrix.csv"))

writeLines(c("Proteomics HC/MCI/AD, UDC-iN ONLY",

 "Matched 72y AD donor excluded (UDC F3:127C and UDC-iN F3:131C); remaining AD UDC-iN donors = 11.",

 "Diagnosis labels matched by text suffix (HC, MCIduetoAD, AD); check Included_Samples.csv.",

 "Model: ~0 + Dx + Age10 + Sex + Batch; Age10 centered per decade.",

 "One UDC-iN sample per donor; no duplicateCorrelation (independent donors).",

 "TMT batch included as a fixed-effect covariate. Age and sex adjusted; common age and sex effects across groups.",

 "No imputation; complete-case across selected UDC-iN samples; remove zero variance.",

 "eBayes(trend=TRUE, robust=FALSE), BH FDR < 0.05; no fold-change cutoff.",

 "Direction is first-named diagnosis minus reference diagnosis: MCI-HC, AD-HC, AD-MCI.",

 "GO BP ORA: UP and DOWN separate; universe = all mapped tested proteins.",

 "GO top 10 chosen by BH adjusted P, ordered by descending GeneRatio."),

 file.path(out,"Analysis_Settings.txt"))

# --------------------------------------------------------------------------

# 4. Volcano + GO ORA, Top10 dot

# --------------------------------------------------------------------------

ratio <- function(s) vapply(strsplit(as.character(s),"/",fixed=TRUE),

  function(z)as.numeric(z[1])/as.numeric(z[2]),numeric(1))

summ <- list()

for (contrast in colnames(cm)) {

  dest <- file.path(out,contrast)

  dir.create(dest,recursive=TRUE,showWarnings=FALSE)

  z <- topTable(fit,coef=contrast,number=Inf,sort.by="none",adjust.method="BH")

  z$ProteinRowID <- rownames(z)

  z <- merge(anno,z,by="ProteinRowID",sort=FALSE)

  z <- z[match(rownames(x),z$ProteinRowID),,drop=FALSE]

  z$Direction <- ifelse(is.finite(z$adj.P.Val)&z$adj.P.Val<fdr&z$logFC>0,"UP",

                  ifelse(is.finite(z$adj.P.Val)&z$adj.P.Val<fdr&z$logFC<0,"DOWN","NS"))

  z$negLog10FDR <- -log10(pmax(z$adj.P.Val,.Machine$double.xmin))

  csv(z,file.path(dest,"All_Proteins.csv"))

  csv(z[z$Direction=="UP",],file.path(dest,"UP_Proteins.csv"))

  csv(z[z$Direction=="DOWN",],file.path(dest,"DOWN_Proteins.csv"))

  zz <- z[is.finite(z$logFC)&is.finite(z$negLog10FDR),]

  labels <- zz[zz$Direction!="NS" & !is.na(zz$SYMBOL)&nzchar(zz$SYMBOL),]

  labels <- labels[order(labels$adj.P.Val,-abs(labels$logFC)),]

  labels <- do.call(rbind,lapply(c("UP","DOWN"),function(d)head(labels[labels$Direction==d,],8)))

  if(is.null(labels)) labels <- zz[FALSE,]

  labels <- labels[!duplicated(labels$SYMBOL),]

  xmax <- max(.5,ceiling(max(abs(zz$logFC))*2)/2)

  title <- paste0("UDC-iN: ",contrast)

  title <- gsub("_vs_"," vs ",title,fixed=TRUE)

  g <- ggplot(zz,aes(logFC,negLog10FDR))+

    geom_point(data=zz[zz$Direction=="NS",],color=gray,size=.85,alpha=.42)+

    geom_point(data=zz[zz$Direction!="NS",],aes(color=Direction),size=1.15,alpha=.8)+

    geom_vline(xintercept=0,linetype="dashed",color="grey55",linewidth=.35)+

    geom_hline(yintercept=-log10(fdr),linetype="dashed",color="grey55",linewidth=.35)+

    ggrepel::geom_text_repel(data=labels,aes(label=SYMBOL,color=Direction),

      size=3,max.overlaps=Inf,box.padding=.4,point.padding=.3,

      segment.color="grey65",min.segment.length=0,seed=123,show.legend=FALSE)+

    scale_color_manual(values=c(UP=red,DOWN=blue))+

    coord_cartesian(xlim=c(-xmax,xmax))+

    labs(title=title,x=expression(log[2]~"abundance fold change"),

         y=expression(-log[10]~"(BH-adjusted P)"))+

    theme_classic(base_size=12)+theme(legend.position="none",

         plot.title=element_text(size=16),plot.margin=margin(15,20,15,15))

  save_fig(g,file.path(dest,paste0("Volcano_",contrast)),8,6.5)

  universe <- unique(na.omit(z$ENTREZID[is.finite(z$adj.P.Val)]))

  go_counts <- c(UP=0L,DOWN=0L)

  for (direction in c("UP","DOWN")) {

    genes <- unique(na.omit(z$ENTREZID[z$Direction==direction]))

    all <- data.frame()

    if(length(genes)>=5) {

      eg <- clusterProfiler::enrichGO(gene=genes,universe=universe,

        OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType="ENTREZID",ont="BP",

        pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,

        minGSSize=10,maxGSSize=500,readable=TRUE)

      if(!is.null(eg)) all <- as.data.frame(eg)

    }

    sig <- if(nrow(all)) all[is.finite(all$p.adjust)&all$p.adjust<go_fdr,] else all

    go_counts[direction] <- nrow(sig)

    top <- if(nrow(sig)) head(sig[order(sig$p.adjust,-sig$Count,sig$ID),],top_n) else sig

    csv(all,file.path(dest,paste0("GO_BP_",direction,"_All.csv")))

    csv(sig,file.path(dest,paste0("GO_BP_",direction,"_Significant.csv")))

    if(nrow(top)) {

      top$GeneRatioNumeric <- ratio(top$GeneRatio)

      top$Term <- stringr::str_wrap(top$Description,width=43)

      top <- top[order(-top$GeneRatioNumeric,top$p.adjust,top$ID),]

      top$TermID <- factor(top$ID,levels=rev(top$ID))

      lab <- stats::setNames(top$Term,top$ID)

      gg <- ggplot(top,aes(GeneRatioNumeric,TermID))+

        geom_point(aes(size=Count,color=p.adjust),alpha=.95)+

        scale_y_discrete(labels=lab)+

        scale_size_continuous(range=c(3.5,9),name="Gene count")+

        scale_color_gradient(low=if(direction=="UP")"#B9473D" else "#235C91",

          high=if(direction=="UP")"#F0C6B8" else "#BFD8ED",

          limits=c(0,go_fdr),name="Adjusted P")+

        scale_x_continuous(expand=expansion(mult=c(.05,.12)))+

        labs(title=paste0(title," — ",direction),x="Gene ratio",y=NULL)+

        theme_classic(base_size=12)+theme(plot.title=element_text(size=16),

          axis.line.y=element_blank(),axis.ticks.y=element_blank(),

          plot.margin=margin(15,20,15,15))

      save_fig(gg,file.path(dest,paste0("GO_BP_",direction,"_Top10_GeneRatio")),10,6)

    } else message(contrast," ",direction,": no significant GO BP terms.")

    csv(top,file.path(dest,paste0("GO_BP_",direction,"_Top10.csv")))

  }

  summ[[contrast]] <- data.frame(Contrast=contrast,Tested=nrow(z),

    UP=sum(z$Direction=="UP"),DOWN=sum(z$Direction=="DOWN"),

    GO_UP=go_counts["UP"],GO_DOWN=go_counts["DOWN"],row.names=NULL)

}

csv(do.call(rbind,summ),file.path(out,"Summary.csv"))

writeLines(c(capture.output(sessionInfo()),paste("Proteomics file:",P),

             paste("Diagnosis mapping:",C)),file.path(out,"SessionInfo.txt"))

message("DONE: ",out)

print(do.call(rbind,summ))
