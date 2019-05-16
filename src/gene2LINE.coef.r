#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)

setwd("../")
args = commandArgs(trailingOnly=TRUE)
resultdir = "result.rsem.TET/gene2L1HS5prime"
cutoff = 0.0005
print (args)
if (length(args)>0) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
  cutoff = args[2]
}
print(resultdir)


cancertypes = list("BLCA", "BRCA", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), "LIHC", c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC")
cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))



coexpressed_plus = c()
coexpressed_minus = c()
for (j in 1:length(cancertypes)) {
  type <- cancertypes[j]
  cancerdir <- cancerdirs[j]
  fname = paste(resultdir, cancerdir, "results_positive_lm.txt", sep="/")
  if (file.exists(fname)) {
    eval <- read.table(fname, sep="\t", header=TRUE)
    eval <- eval[order(-eval$partialeta2),]
    minuscoef <- eval[eval$generatio_coef<0,]
    pluscoef <- eval[eval$generatio_coef>0,]
    dim(minuscoef)
    dim(pluscoef)
    sigplus <- pluscoef[pluscoef$generatio_qval<cutoff,]
    sigminus <- minuscoef[minuscoef$generatio_qval<cutoff,]
    coexpressed_plus = rbind.data.frame(coexpressed_plus, cbind.data.frame(rep(cancerdir, nrow(sigplus)), sigplus))
    coexpressed_minus = rbind.data.frame(coexpressed_minus, cbind.data.frame(rep(cancerdir, nrow(sigminus)), sigminus))
  }
}
write.table(coexpressed_plus, file=paste(resultdir, "coexpressed_plus.txt", sep="/"), sep="\t", quote = FALSE, row.names = FALSE )
write.table(coexpressed_minus, file=paste(resultdir, "coexpressed_minus.txt", sep="/"), sep="\t", quote = FALSE, row.names = FALSE )



