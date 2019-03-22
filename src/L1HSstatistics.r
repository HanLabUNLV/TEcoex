library(stringr)
library(ggplot2)
library(DESeq2)

setwd("../")
args = commandArgs(trailingOnly=TRUE)
genedatadir = "data/rsem/" # default gene data dir
TEdatadir = "data/TET/" # default TE data dir
print (args)
if (length(args)==2) {
  genedatadir = args[1]
  if (substr(genedatadir, nchar(genedatadir), nchar(genedatadir)) != "/") {
    genedatadir = paste0(genedatadir, "/")
  }
  TEdatadir = args[2]
  if (substr(TEdatadir, nchar(TEdatadir), nchar(TEdatadir)) != "/") {
    TEdatadir = paste0(TEdatadir, "/")
  }
}
resultdir = paste0("result.", unlist(str_split(genedatadir, "/"))[2], ".", unlist(str_split(TEdatadir, "/"))[2], "/")
dir.create(resultdir) 
print(resultdir)


if (file.exists(paste0(resultdir,"ddsGenes.RData"))) {
  load(paste0(resultdir,"ddsGenes.RData"))
} else {
  
  if (grepl("TET", genedatadir)) {
    file.names <- dir(genedatadir, pattern ="ele.cntTable", recursive=TRUE)
    file.names <- file.names[!grepl ("discount", file.names)]
    file.names <- file.names[!grepl ("v1", file.names)]
    
    tissueID <- data.frame(matrix(unlist(str_split(file.names, "/")), ncol=2, byrow=TRUE))
    tissues <- factor(tissueID[,1])
    patientIDs <- factor(data.frame(matrix(unlist(str_split(tissueID[,2], "[.]", n=2)), ncol=2, byrow=TRUE))[,1])
  } else if (grepl("rsem", genedatadir)) {
    file.names <- dir(genedatadir, pattern ="genes.results", recursive=TRUE)
    names <- data.frame(matrix(unlist(strsplit(file.names, "/")), ncol=2, byrow=TRUE))
    tissues <- factor(names[,1])
    patientIDs <- factor(substr(names[,2], 1, 12))
  }
  file.names = paste0(genedatadir,file.names)
  patient_files <- cbind.data.frame(patientIDs, file.names)
  colnames(patient_files) <- c("patient", "filename")
  
  genenames <- rownames(read.table(as.character(patient_files[1,2]), header=TRUE, row.names=1))  
  TEidx <- grepl (":", genenames)
  geneidx <-!TEidx
  coldata <- cbind.data.frame(patientIDs, tissues)
  colnames(coldata) = c("patient", "tissue")
  batchids <- read.table("data/patient.info", header=TRUE, sep="\t")
  rownames(batchids) <- batchids$patient
  batchids <- batchids[coldata$patient,]
  coldata <- cbind.data.frame(coldata, batchids[-c(1,2)])
  rownames(coldata) <- paste(coldata$patient, coldata$tissue, sep=".")
  rm(patientIDs, tissues)
  
  coldata$BatchId <- factor(coldata$BatchId)
  #need to make batch nested within tissue 
#  coldata <- cbind.data.frame(coldata, rep(0, nrow(coldata)))
#  colnames(coldata)[ncol(coldata)] <- "nested.batch"
#  tissuetypes = levels(coldata$tissue)
#  for (i in 1:length(tissuetypes)) {
#    coltissue <- coldata[coldata$tissue==tissuetypes[i],]
#    coltissue$BatchId <- droplevels(coltissue$BatchId)
#    coldata[coldata$tissue==tissuetypes[i],"nested.batch"]=as.numeric(coltissue$BatchId)
#  }
  coldata$nested.batch = factor(coldata$nested.batch)
  dim(coldata)
  
  
  cntmatrix = matrix(rep(0, nrow(patient_files)*sum(geneidx)),  nrow=sum(geneidx))
  for (i in 1:nrow(patient_files)) {
    cnts <- read.table(as.character(patient_files[i,2]), header=TRUE, row.names=1)
    cntmatrix[,i]=as.matrix(cnts[geneidx,1])
  }
  dim(cntmatrix)
  rownames(cntmatrix) <- genenames[geneidx]
  colnames(cntmatrix) <- rownames(coldata)
  write.table(cntmatrix, file=paste0(resultdir,"cntmatrix.gene.txt"), quote=FALSE, row.names=TRUE, sep="\t")
 

  if (file.exists(paste0(resultdir, "nonDEGgenes.txt"))) {
    nonDEGgenes <- read.table(file=paste0(resultdir,"nonDEGgenes.txt"))
    ctrlgenesidx <-  match(nonDEGgenes[,1], rownames(cntmatrix))

  } else {
 
  ddsGenes <- DESeqDataSetFromMatrix(
    countData = round(cntmatrix),
    colData = coldata,
    design =  ~ nested.batch + tissue
  )
  

  # step 1.1
  FDR <- 0.1
  floorPDEG <- 0.05
  filter <- apply(counts(ddsGenes), 1, function(x) { all(x > 0)})
  ddsGenes <- ddsGenes[filter,]
  
  
  if (file.exists(paste0(resultdir,"sizefactors.0.txt"))) {
    sf <- read.table(paste0(resultdir,"sizefactors.0.txt"), row.names=1, sep="\t", header=TRUE) 
    sizeFactors(ddsGenes) <- c(t(sf))
  } else {
    ddsGenes <- estimateSizeFactors(ddsGenes)
    sizeFactors(ddsGenes) <- sizeFactors(ddsGenes) / mean(sizeFactors(ddsGenes))
    sf <- sizeFactors(ddsGenes)
    write.table(sf, file=paste0(resultdir,"sizefactors.0.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  }
  
  # step 1.2
  if (file.exists(paste0(resultdir,"ddsGenes.LRT.RData"))) {
    load(paste0(resultdir,"ddsGenes.LRT.RData"))
  } else {
    ddsGenes <- estimateDispersions(ddsGenes)
    ddsGenes <- nbinomLRT(ddsGenes, reduced = ~ nested.batch, maxit=1000)
    if (any(mcols(ddsGenes)$betaConv)) {
      ddsGenes <- ddsGenes[which(mcols(ddsGenes)$betaConv),]
    }
    save(ddsGenes, file=paste0(resultdir,"ddsGenes.LRT.RData"))
  }
  
  if (file.exists(paste0(resultdir,"res.LRT.RData"))) {
    load(paste0(resultdir,"res.LRT.RData"))
  } else {
    res <- results(ddsGenes)
    save(res, file=paste0(resultdir,"res.LRT.RData"))
  }
  
  FDR = 1.0e-50
  floorPDEG = 1.0e-50
  pval <- res$pvalue
  pval[is.na(pval)] <- 1
  qval <- p.adjust(pval, method = "BH")
  if (sum(qval < FDR) > (floorPDEG * nrow(ddsGenes))) {
    is.DEG <- as.logical(qval < FDR)
  } else {
    is.DEG <- as.logical(rank(pval, ties.method = "min") <= nrow(ddsGenes) * floorPDEG)
  }
  
  # step 1.3
  nonDEGgenes <- rownames(counts(ddsGenes, normalized=FALSE)[!is.DEG,])
  write.table(nonDEGgenes, file=paste0(resultdir,"nonDEGgenes.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  ctrlgenesidx <-  match(nonDEGgenes, rownames(cntmatrix))

  } 

  ddsGenes <- DESeqDataSetFromMatrix(
    countData = round(cntmatrix),
    colData = coldata,
    design = ~ nested.batch + tissue )
  ddsGenes <- estimateSizeFactors(ddsGenes, controlGenes = ctrlgenesidx)

  sizefactors.1 <- sizeFactors(ddsGenes)
  write.table(sizefactors.1, file=paste0(resultdir,"sizefactors.1.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  normcnts.1 = counts(ddsGenes, normalized=TRUE)
 


  # step 2 normalize within tissue 
  sizefactors.2 <- rep(1, length(sizefactors.1))
  normcnts.2 = matrix(rep(0, nrow(cntmatrix)*ncol(cntmatrix)),  nrow=nrow(cntmatrix))
  tissuelist <- levels(coldata$tissue)
  for (i in 1:length(tissuelist)) {
    keep <- coldata$tissue == tissuelist[i]
    dds <- DESeqDataSetFromMatrix(
      countData = round(normcnts.1[,keep]),
      colData = coldata[keep,],
      design = ~ 1 )
    dds <- estimateSizeFactors(dds)
    sizefactors.2[keep] = sizeFactors(dds)
    normcnts.2[,keep] = counts(dds, normalized=TRUE)
  }
  write.table(sizefactors.2, file=paste0(resultdir,"sizefactors.2.txt"), quote=FALSE, row.names=TRUE, sep="\t")

  stopifnot (normcnts.2 == t(round(t(round(cntmatrix))/sizefactors.1)/sizefactors.2))
  sf <- sizefactors.1*sizefactors.2
  sizeFactors(ddsGenes) <- sf
  write.table(sf, file=paste0(resultdir,"sizefactors.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  save(ddsGenes, file=paste0(resultdir,"ddsGenes.RData"))
  
}



if (file.exists(paste0(resultdir,"ddsNew.RData"))) {
  load(file = paste0(resultdir,"ddsNew.RData"))
} else {
  
  sf = sizeFactors(ddsGenes)
  coldata = colData(ddsGenes)
  cntmatrix = read.table(file=paste0(resultdir,"cntmatrix.gene.txt"), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
  
  # construct counts from discounted elements 
  file.names <- dir(TEdatadir, pattern ="ele.cntTable", recursive=TRUE)
  file.names <- file.names[grepl ("discount", file.names)]
  file.names <- file.names[!grepl ("v1", file.names)]
  
  tissueID <- data.frame(matrix(unlist(str_split(file.names, "/")), ncol=2, byrow=TRUE))
  tissues <- factor(tissueID[,1])
  patientIDs <- factor(data.frame(matrix(unlist(str_split(tissueID[,2], "[.]", n=2)), ncol=2, byrow=TRUE))[,1])
  file.names = paste0(TEdatadir,file.names)
  patient_files <- cbind.data.frame(patientIDs, file.names)
  colnames(patient_files) <- c("patient", "filename")
  TEnames <- rownames(read.table(as.character(file.names[1]), header=TRUE, row.names=1))
  
  TEcntmatrix = matrix(rep(0, length(file.names)*length(TEnames)),  nrow=length(TEnames))
  for (i in 1:length(file.names)) {
    cnts <- read.table(as.character(file.names[i]), header=TRUE, row.names=1)
    TEcntmatrix[,i]=as.matrix(cnts)
  }
  dim(TEcntmatrix)
  rownames(TEcntmatrix) <- TEnames
  colnames(TEcntmatrix) <- paste(patientIDs, tissues, sep=".")
  
  # merge with existing gene cntmatrix
  t_TEcntmatrix = t(TEcntmatrix)
  t_cntmatrix = t(cntmatrix)
  
  m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
  new_patients = m[,1]
  new_genenames = colnames(m[,-1])
  new_cntmatrix = t(m[,-1])
  colnames(new_cntmatrix) = new_patients
  write.table(new_cntmatrix, file=paste0(resultdir,"cntmatrix.geneTE.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  m_sf <- as.matrix(sf)
  new_sf <- m_sf[new_patients,]  # need patient ids in same order to use previously estimated sizeFactors
  
  stopifnot( names(new_sf) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
  new_coldata <- coldata[names(new_sf),]
  colnames(new_cntmatrix) = NULL
  ddsnew <- DESeqDataSetFromMatrix(
    countData = round(new_cntmatrix),
    colData = new_coldata,
    design = ~ nested.batch + tissue )
  # set sizefactors 
  sizeFactors(ddsnew) <- new_sf
  save(ddsnew, file=paste0(resultdir,"ddsNew.RData"))
  ddsnewcnts <- counts (ddsnew, normalized=TRUE)
  write.table(ddsnewcnts, file=paste0(resultdir,"normalizedcnt.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  
  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
  write.table(L1HSnormalized, file=paste0(resultdir,"L1HS.normalized.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  
  normalizedcolsums <- colSums(ddsnewcnts)
  log2colsums <- log2(normalizedcolsums)
  write.table(log2colsums, file=paste0(resultdir,"log2colsums.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  
}



 
 
 if (file.exists(paste0(resultdir,"VSTcnts.DEGES.RData"))) {
  load(file = paste0(resultdir,"VSTcnts.DEGES.RData"))
 } else {
  vsd <- vst(ddsnew, blind = FALSE)
  VSTcnts <- assay(vsd)
  write.table(VSTcnts["L1HS:L1:LINE",], paste0(resultdir,"L1HS.VST.cnts.txt"), quote=FALSE, row.names = TRUE)
  save(VSTcnts, file = paste0(resultdir,"VSTcnts.DEGES.RData"))
  write.table(VSTcnts, file=paste0(resultdir,"VSTcnt.txt"), quote=FALSE, row.names=TRUE, sep="\t")
 
  dir.create(paste0(resultdir, "VSTcnts")) 
  genenames = rownames(ddsnew)
  genenames <- gsub("ALR/Alpha", "ALR_Alpha", genenames)
  genenames <- gsub("BSR/Beta", "BSR_Beta", genenames)
  genenames <- gsub("THRA1/BTR", "THRA1_BTR", genenames)
  rownames(ddsnew) <- genenames
  coldataNew = colData(ddsnew)
  for (i in 1:nrow(VSTcnts)) {
    fname = paste0(resultdir,"VSTcnts/",genenames[i], ".txt")
    cnts <- matrix(VSTcnts[i,], ncol=1)
    rownames(cnts) <- coldataNew$patient
    colnames(cnts) <- c("patient\tgene")
    write.table(cnts, file=fname, quote=FALSE, sep="\t", row.names=TRUE)
  }
  L1HSnormalized = read.table(paste0(resultdir, "L1HS.normalized.txt"), row.names=1, sep="\t", check.names = FALSE)
  L1HStable <- cbind.data.frame(coldataNew, L1HSnormalized )
  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
  colnames(L1HStable)[ncol(coldataNew)+1:2] = c("normalizedcnt", "VSTcnts")
  write.table(L1HStable, file=paste0(resultdir,"L1HS.VST.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  
  lineVSTcnts <- VSTcnts[grep("LINE", rownames(VSTcnts)), ]
  LINEcnt <- 2^(lineVSTcnts)
  LINEtotal <- colSums(LINEcnt)
  L1HScnts <- LINEcnt["L1HS:L1:LINE",]
  L1PA2cnts <- LINEcnt["L1PA2:L1:LINE",]
  L1PA3cnts <- LINEcnt["L1PA3:L1:LINE",]
  L1PA4cnts <- LINEcnt["L1PA4:L1:LINE",]
  L1PA5cnts <- LINEcnt["L1PA5:L1:LINE",]
  LINEyoung <- L1HScnts+L1PA2cnts+L1PA3cnts+L1PA4cnts+L1PA5cnts
  LINEold <- LINEtotal - LINEyoung
  LINEVSTold <- cbind.data.frame(coldataNew$patient, log2(LINEold))
  write.table(LINEVSTold, file=paste0(resultdir,"oldLINE.VST.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  
 }


 
# violin plot L1HS
L1HS_bytype = read.table(paste0(resultdir,"L1HS.VST.txt"), header=TRUE, row.names=1, sep="\t")
oldLINE = read.table(paste0(resultdir,"oldLINE.VST.txt"), header=TRUE, row.names=1, sep="\t")
colnames(oldLINE) = c("patient", "oldLINEVSTcnts")

stopifnot(! any(oldLINE$patient!=L1HS_bytype$patient))
oldLINE <- cbind(oldLINE[,2], L1HS_bytype[,1:2])
colnames(oldLINE) <- colnames(L1HS_bytype[,c(36, 1, 2)])

L1HS_5prime = read.table(paste0(TEdatadir, "L1HS5prime.300.txt"), header=FALSE, sep="\t", stringsAsFactor=FALSE)
rownames(L1HS_5prime)= L1HS_5prime[,1]
L1HS_5prime = L1HS_5prime[as.character(L1HS_bytype$patient),]

stopifnot(! any(rownames(L1HS_5prime)!=L1HS_bytype$patient))
L1HS_5prime <- cbind(log2(L1HS_5prime[,2]), L1HS_bytype[,c(1,2,4,5)])
L1HS_5prime[L1HS_5prime[,1]<0,1] = -1
colnames(L1HS_5prime) <- colnames(L1HS_bytype[,c(36, 1,2,4,5)])
write.table(L1HS_5prime, file=paste0(resultdir,"L1HS5prime.VST.txt"), quote=FALSE, row.names=TRUE, sep="\t")
L1HS_5prime <- L1HS_5prime[,1:3]

housekeepinggenes <- read.table("../mobiledna.housekeeping/data/housekeeping.txt")
gnames <- unlist(strsplit(rownames(VSTcnts), "[|]"))[c(TRUE, FALSE)]
ctrlgenesidx <-  match(housekeepinggenes[,1], gnames)
ctrlgenesidx <- ctrlgenesidx[!is.na(ctrlgenesidx)]
housemeans <- cbind(colMeans(VSTcnts[ctrlgenesidx,]), L1HS_bytype[,1:2])
stopifnot(! any(rownames(housemeans) != rownames(L1HS_bytype)))
colnames(housemeans) <- colnames(L1HS_bytype[,c(36, 1, 2)])

L1_house <- rbind(housemeans, L1HS_5prime)
L1_house$L1 = c(rep(0, nrow(housemeans)), rep(2, nrow(L1HS_5prime)))
L1_house$L1 <- as.factor(L1_house$L1)

dodge <- position_dodge(width = 0.6)
ggplot(L1_house, aes(x=tissue, y=VSTcnts, fill=L1)) + geom_violin(position = dodge, draw_quantiles = 0.5) + 
    scale_fill_manual(values=c("#A4A4A4", "#E69F00", "#FF6666", "#D55E00")) +
    theme_bw()
ggsave(paste0(resultdir, 'violin.bytype.L1HS.pdf'), width = 20, height = 7, dpi=600)


# scatter plot L1HS5prime vs oldLINE
L1HS_oldLINE <- cbind.data.frame(L1HS_5prime[,1], oldLINE)
colnames(L1HS_oldLINE)[1] <- "L1HS"
bp <- ggplot(L1HS_oldLINE, aes(x=L1HS, y=VSTcnts)) + geom_point() 
pdf(paste0(resultdir, "L1HS5prime2oldLINE.pdf"), height=28, width=7)
bp + facet_grid(tissue ~ .)
dev.off()
L1exp = L1HS_oldLINE
L1exp[,1:2] = 2^L1HS_oldLINE[,1:2]
bp <- ggplot(L1exp, aes(x=L1HS, y=VSTcnts)) + geom_point() 
pdf(paste0(resultdir, "L1HS5prime2oldLINE.exp.pdf"), height=28, width=7)
bp + facet_grid(tissue ~ .)
dev.off()


# violin plot family means
TEmeans <- function(familyname)
{
  LINEsidx <-  grepl(familyname, gnames)
  LINEsidx <- LINEsidx[!is.na(LINEsidx)]
  LINEmeans <- cbind(colMeans(VSTcnts[LINEsidx,]), L1HS_bytype[,1:2])
  stopifnot(! any(rownames(LINEmeans) != rownames(L1HS_bytype)))
  colnames(LINEmeans) <- colnames(L1HS_bytype[,c(36, 1, 2)])
  LINEmeans
}

LINEmeans <- TEmeans(":LINE")
DNAmeans <- TEmeans(":DNA")
LTRmeans <- TEmeans(":LTR")
SINEmeans <- TEmeans(":SINE")

family_house <- rbind(housemeans, LINEmeans, DNAmeans, LTRmeans, SINEmeans)
family_house$fam = c(rep(0, nrow(housemeans)), rep(1, nrow(LINEmeans)), rep(2, nrow(DNAmeans)), rep(3, nrow(LTRmeans)), rep(4, nrow(SINEmeans)) )
family_house$fam <- as.factor(family_house$fam)

dodge <- position_dodge(width = 0.6)
ggplot(family_house, aes(x=tissue, y=VSTcnts, fill=fam)) + geom_violin(position = dodge, draw_quantiles = 0.5) +
      scale_fill_manual(values=c("#A4A4A4", "#FDBF6F", "#FB9A99", "#B2DF8A", "#CAB2D6")) +
      #stat_summary(fun.y=mean, geom="point", shape=20, size=7, color="red", fill="red") +
      theme_bw()
ggsave(paste0(resultdir, 'violin.bytype.family.pdf'), width = 20, height = 7, dpi=600)

exit()



 

  if (file.exists(paste0(resultdir,"ddsNew.LRT.RData"))) {
    load(paste0(resultdir,"ddsNew.LRT.RData"))
  } else {
    ddsnew <- estimateDispersions(ddsnew)
    ddsnew <- nbinomLRT(ddsnew, reduced = ~ nested.batch, maxit=2000)
    if (any(mcols(ddsnew)$betaConv)) {
      ddsnew <- ddsnew[which(mcols(ddsnew)$betaConv),]
    }
    save(ddsnew, file=paste0(resultdir,"ddsNew.LRT.RData"))
  }

  #contrasts <- resultsNames(ddsnew)
  coldataNew = colData(ddsnew)
  tissuelist <- levels(coldataNew$tissue)
  dir.create(paste0(resultdir, "resNEW.LRT")) 
  for (i in 1:length(tissuelist)) {
    for (j in (i+1):length(tissuelist)) {
      filename = paste0(resultdir,"resNEW.LRT/resNew.LRT.", i, ".", j, ".RData")
      if (file.exists(filename)) {
        load(filename)
      } else {
        res <- results(ddsnew, contrast=c("tissue", as.character(tissuelist[i]), as.character(tissuelist[j])))
        save(res, file=filename)
      }
    }
  }

 
  for (i in 1:length(tissuelist)) {
    for (j in (i+1):length(tissuelist)) {
      resfilename = paste0(resultdir,"resNEW.LRT/resNew.LRT.", i, ".", j, ".RData")
      pdffilename = paste0(resultdir,"resNEW.LRT/MAplot.", tissuelist[i], ".", tissuelist[j], ".pdf")
      if (file.exists(resfilename)) {
        load(resfilename)
        pdf(pdffilename)
        plotMA(res)
        dev.off()
      }
    }
  }




