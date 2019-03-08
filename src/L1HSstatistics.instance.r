library(stringr)
library(ggplot2)
library(DESeq2)
library(pheatmap)

#library("RColorBrewer")
# library("DESeq2")
# library("geneplotter")
# library("EDASeq")
# library("genefilter")


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
resultdir = paste0("result.", unlist(str_split(genedatadir, "/"))[2], ".", unlist(str_split(TEdatadir, "/"))[2], ".instance",  "/")
generesultdir = str_replace(resultdir, ".instance", "")
dir.create(resultdir) 
print(resultdir)
print(generesultdir)


if (file.exists(paste0(generesultdir,"ddsGenes.RData"))) {
  load(paste0(generesultdir,"ddsGenes.RData"))
} else  {
  print("run L1HSstatistics.r before running L1HSstatistics.instance.r")
  quit()
}


if (file.exists(paste0(resultdir,"ddsNew.RData"))) {
  load(file = paste0(resultdir,"ddsNew.RData"))
} else {
  
  sf = sizeFactors(ddsGenes)
  coldata = colData(ddsGenes)
  cntmatrix = read.table(file=paste0(generesultdir,"cntmatrix.gene.txt"), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
  
  # construct counts from discounted instances 
  TEinstances <- read.table(paste0(TEdatadir, "rawcnts.new.maxgt5.txt"), header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  t_TEinstances = t(TEinstances)
  rownames(t_TEinstances) <- colnames(TEinstances)
  t_cntmatrix = t(cntmatrix)
  rownames(t_cntmatrix) <- coldata$patient

  # merge with existing gene cntmatrix
  m <- merge(t_cntmatrix, t_TEinstances, by.x="row.names", by.y="row.names")
  new_patients = m[,1]
  new_genenames = colnames(m[,-1])
  new_cntmatrix = t(m[,-1])
  colnames(new_cntmatrix) = new_patients
  write.table(new_cntmatrix, file=paste0(resultdir,"cntmatrix.geneTE.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  m_sf <- as.matrix(sf)
  sf_patient_tissue <- data.frame(matrix(unlist(strsplit(rownames(m_sf), "[.]")), ncol=2, byrow=TRUE))
  m_sf <- cbind.data.frame(m_sf, sf_patient_tissue)
  rownames(m_sf) <- sf_patient_tissue[,1]
  m_sf <- m_sf[new_patients,]  # need patient ids in same order to use previously estimated sizeFactors
  new_sf <- m_sf[,1]  # need patient ids in same order to use previously estimated sizeFactors
  names(new_sf) <- paste0(m_sf[,2], ".", m_sf[,3])
  
  
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
  

}



 
 
 if (file.exists(paste0(resultdir,"VSTcnts.DEGES.RData"))) {
  load(file = paste0(resultdir,"VSTcnts.DEGES.RData"))
 } else {
  vsd <- vst(ddsnew, blind = FALSE)
  VSTcnts <- assay(vsd)
  save(VSTcnts, file = paste0(resultdir,"VSTcnts.DEGES.RData"))
  write.table(VSTcnts, file=paste0(resultdir,"VSTcnt.txt"), quote=FALSE, row.names=TRUE, sep="\t")
  dir.create(paste0(resultdir, "VSTcnts")) 
  
  genenames = rownames(ddsnew)
  genenames <- gsub("ALR/Alpha", "ALR_Alpha", genenames)
  genenames <- gsub("BSR/Beta", "BSR_Beta", genenames)
  genenames <- gsub("THRA1/BTR", "THRA1_BTR", genenames)
  genenames <- substr(genenames, 1, nchar(genenames)-2)
  rownames(ddsnew) <- genenames
  coldataNew = colData(ddsnew)
  for (i in 1:nrow(VSTcnts)) {
    fname = paste0(resultdir,"VSTcnts/",genenames[i], ".txt")
    cnts <- matrix(VSTcnts[i,], ncol=1)
    rownames(cnts) <- coldataNew$patient
    colnames(cnts) <- c("patient\tgene")
    write.table(cnts, file=fname, quote=FALSE, sep="\t", row.names=TRUE)
  }
 }




