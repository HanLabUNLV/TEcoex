#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)

setwd("../../")
datadir = "data/"
resultdir = "result.rsem.TET.instance/"
generesultdir = "result.rsem.TET/"
args = commandArgs(trailingOnly = TRUE)

#Enter the L1family, ZFPgene, and cancertype arguments exactly how they appear on the L1 filename/ZFPgene filename/cancertype column.
#The directories are hard coded, so they may have to be modified if you're running different tests.

#ZNF662|389114_L2a:L2:LINE_KIRP_coexpressed_minus.txt
#BRCA/ZNF641|121274_MER97a:hAT-Tip100:DNA_results_positive_lm.txt
L1family = "MER97a:hAT-Tip100:DNA"
ZFPgene = "ZNF641|121274"
cancertype = "BRCA"
if (length(args) > 0) {
  L1family = args[1]
  ZFPgene = args[2]
  cancertype = args[3]
}
print (L1family)
print (ZFPgene)
print (cancertype)


cancerdir = cancertype
cancername = cancertype



coexpressed_plus = c()
coexpressed_minus = c()
fname = paste(resultdir, "ZNF/", cancerdir, "/", ZFPgene, "_", L1family, "_results_positive_lm.txt", sep="")
if (file.exists(fname)) {
  eval <- read.table(fname, sep="\t", header=TRUE)
  eval <- eval[order(-eval$partialeta2),]
  minuscoef <- eval[eval$generatio_coef<0,]
  pluscoef <- eval[eval$generatio_coef>0,]
  dim(minuscoef)
  dim(pluscoef)
  sigplus <- pluscoef[pluscoef$generatio_qval<0.00001,]
  sigminus <- minuscoef[minuscoef$generatio_qval<0.00001,]
  coexpressed_plus = rbind.data.frame(coexpressed_plus, cbind.data.frame(rep(cancerdir, nrow(sigplus)), sigplus))
  coexpressed_minus = rbind.data.frame(coexpressed_minus, cbind.data.frame(rep(cancerdir, nrow(sigminus)), sigminus))
}

if (nrow(coexpressed_plus)>0) {
  filename = paste(resultdir, "ZNF/", ZFPgene, "_", L1family, "_", cancername, "_coexpressed_plus.txt", sep="")
  write.table(coexpressed_plus, file=filename, sep="\t", quote = FALSE, row.names = FALSE )
}
if (nrow(coexpressed_minus)>0) {
  filename = paste(resultdir, "ZNF/", ZFPgene, "_", L1family, "_", cancername, "_coexpressed_minus.txt", sep="")
  write.table(coexpressed_minus, file=filename, sep="\t", quote = FALSE, row.names = FALSE )
}

