#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)

setwd("../")
args = commandArgs(trailingOnly=TRUE)
resultdir = "result.rsem.TET/"
print (args)
if (length(args)>0) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
}
print(resultdir)
subdir = "gene2L1HS5prime"





generatio <- function(genename)
{
  normalfname <- paste(resultdir, "VSTcnts/", genename, ".txt", sep="")
  
  if (file.info(normalfname)$size == 0) {
    print ("file size zero")
    return (0)
  }
  gene <- read.table(normalfname, header=TRUE)
  colnames(gene) <- c("patient", "gene_n")
  return (gene)
}



L1HS <- read.table(paste0(resultdir, "L1HS5prime.VST.txt"), sep="\t", header=TRUE, row.names=1)
#L1HS <- L1HS[L1HS$condition=="normal",]
#colnames(L1HS)[7:8] <-c("L1HSnRPM", "VSTcnts")

log2TPMsum <- read.table(file=paste0(resultdir, "log2colsums.txt"), sep="\t", header = TRUE, row.names=1)
log2TPMsum <- cbind(substr(rownames(log2TPMsum), 1, 12), log2TPMsum)
colnames(log2TPMsum) <- c("patient", "log2TPMsum")

#oldLINEVST <- read.table(file=paste0(resultdir, "oldLINE.VST.txt"), sep="\t", header = TRUE, row.names=1)
#colnames(oldLINEVST) <- c("patient", "oldLINEVST")
oldLINEVST <- read.table(file=paste0(resultdir, "/WGCNA/consensus/eigengeneME1.txt"), sep="\t", header = FALSE)
colnames(oldLINEVST) <- c("patient", "oldLINEVST")

file.names <- dir(paste0(resultdir, "VSTcnts/"), pattern =".txt")
file.names.nchar <- sapply(file.names, nchar)
genenames <- substr(file.names, 1, file.names.nchar-4)
genenames <- unique(genenames)

#cancertypes = list("BLCA", "BRCA", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), "LIHC", c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC")
cancertypes = list("BLCA", "BRCA", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), "LIHC", c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC", "LAML", "ACC")
cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))
dir.create(paste0(resultdir, subdir))
newdirs = paste(resultdir, subdir, cancerdirs, sep="/")
sapply(newdirs, dir.create)

THCAradiation = read.table(file="data/THCA.radiation.txt", sep="\t", header=FALSE, row.names=1)
THCAradiation <- cbind.data.frame(substr(rownames(THCAradiation), 1, 12), THCAradiation)
colnames(THCAradiation) <- c("patient", "THCAradiation")


for (j in 1:length(cancertypes)) {
  cancertype <- cancertypes[[j]]
  cancerdir <- cancerdirs[j]
  filename = paste(resultdir, subdir, cancerdir, "results_positive_lm.txt", sep="/")
  print(cancerdir)

  subdata <- L1HS[L1HS$tissue %in% cancertype,]
  if (nrow(subdata)<8) {
    print(paste0("skipping ", cancerdir))
    next;
  }

  pearson_r <- 0
  pearson_pval <- 1
  spearman_rho <- 0
  spearman_pval <- 1
  bestmodelnames <- "mod"
  fixed_generatio <- 0 
  generatio_coef <- 0 
  generatio_pval <- 1 
  partialeta2 <- 0 
  tested <- 0 
  
#  for(i in 1:length(genenames)){
#    print (i)
#    print (genenames[i])
    XBP1 <- read.table("./data/XBP1.isoforms.TCGA.rsem.rawcnts.txt", header=TRUE)
    XBP1 <- merge(L1HS[,c("patient", "VSTcnts", "sizeFactor")], XBP1, by = "patient")
    XBP1 = cbind.data.frame(XBP1, XBP1u=rowSums(XBP1[,4:6])/XBP1$sizeFactor)
    XBP1 = cbind.data.frame(XBP1, XBP1s=XBP1$uc011akl.1/XBP1$sizeFactor)
    XBP1 = cbind.data.frame(XBP1, gene_n=(XBP1$XBP1u/(XBP1$XBP1u+XBP1$XBP1s)))
    XBP1psi <- read.table("./data/XBP1.totalpsi.MDanderson.txt", header=TRUE)
    XBP1 = merge(XBP1, XBP1psi, by="patient")
    XBP1$totalpsi = as.numeric(as.character(XBP1$totalpsi))
    XBP1 = XBP1[!is.na(XBP1$totalpsi),]
    #XBP1$gene_n = XBP1$totalpsi
    plot( totalpsi ~ gene_n, data=XBP1)
 
#    gene<- generatio(genenames[i])
#    if (is.null(dim(gene))) next;

    data = merge(L1HS, XBP1[,c(1,4:11)], by="patient")
    data <- merge(data, log2TPMsum, by = "patient")
    data <- merge(data, oldLINEVST, by = "patient")
    if (cancertype == "THCA") {
      data <- merge(data, THCAradiation, by = "patient")
    }
    
    subdata <- data[data$tissue %in% cancertype,]
    if (nrow(subdata)<8) next;
    subdata$tissue <- droplevels(subdata$tissue)
    subdata$BatchId <- factor(subdata$BatchId)
    if (length(cancertype)>1) {
      subdata$BatchId <- droplevels(subdata$BatchId)
      subdata$nested.batch=as.numeric(subdata$BatchId)
    }
    subdata$nested.batch = factor(subdata$nested.batch)
    dataratio <- subdata    
    
    test<-cor.test(dataratio$gene_n,dataratio$VSTcnts, method="pearson")
    pearson_r <- test$estimate
    pearson_pval <- test$p.value
    test<-cor.test(dataratio$gene_n,dataratio$VSTcnts, method="spearman")
    spearman_rho <- test$estimate
    spearman_pval <- test$p.value

    AICc = NULL    
    set.seed(0)
    garbage <- rnorm(length(dataratio$VSTcnts))  
    if (cancertype == "THCA") {
      if (length(summary(dataratio$nested.batch))>1) {
        mod0 <- lm(VSTcnts ~ garbage + oldLINEVST, data=dataratio)
        mod1 <- lm(VSTcnts ~ log2TPMsum + oldLINEVST, data=dataratio)
        mod2 <- lm(VSTcnts ~ gene_n + oldLINEVST, data=dataratio)
        mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
        mod4 <- lm(VSTcnts ~ nested.batch + oldLINEVST, data=dataratio)
        mod5 <- lm(VSTcnts ~ log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
        mod6 <- lm(VSTcnts ~ gene_n + nested.batch + oldLINEVST, data=dataratio)
        mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
        mod8 <- lm(VSTcnts ~ THCAradiation + oldLINEVST, data=dataratio)
        mod9 <- lm(VSTcnts ~ log2TPMsum + THCAradiation + oldLINEVST, data=dataratio)
        mod10 <- lm(VSTcnts ~ gene_n + THCAradiation + oldLINEVST, data=dataratio)
        mod11 <- lm(VSTcnts ~ gene_n + log2TPMsum + THCAradiation + oldLINEVST, data=dataratio)
        mod12 <- lm(VSTcnts ~ nested.batch + THCAradiation + oldLINEVST, data=dataratio)
        mod13 <- lm(VSTcnts ~ log2TPMsum + nested.batch + THCAradiation + oldLINEVST, data=dataratio)
        mod14 <- lm(VSTcnts ~ gene_n + nested.batch + THCAradiation + oldLINEVST, data=dataratio)
        mod15 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch + THCAradiation + oldLINEVST, data=dataratio)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
      } else {
        mod0 <- lm(VSTcnts ~ garbage + oldLINEVST, data=dataratio)
        mod1 <- lm(VSTcnts ~ log2TPMsum + oldLINEVST, data=dataratio)
        mod2 <- lm(VSTcnts ~ gene_n + oldLINEVST, data=dataratio)
        mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
        mod4 <- lm(VSTcnts ~ THCAradiation + oldLINEVST, data=dataratio)
        mod5 <- lm(VSTcnts ~ log2TPMsum + THCAradiation + oldLINEVST, data=dataratio)
        mod6 <- lm(VSTcnts ~ gene_n + THCAradiation + oldLINEVST, data=dataratio)
        mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + THCAradiation + oldLINEVST, data=dataratio)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      }      
    } else {
    if (length(summary(dataratio$nested.batch))>1) {
      mod0 <- lm(VSTcnts ~ garbage + oldLINEVST, data=dataratio)
      mod1 <- lm(VSTcnts ~ log2TPMsum + oldLINEVST, data=dataratio)
      mod2 <- lm(VSTcnts ~ gene_n + oldLINEVST, data=dataratio)
      mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
      mod4 <- lm(VSTcnts ~ nested.batch + oldLINEVST, data=dataratio)
      mod5 <- lm(VSTcnts ~ log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
      mod6 <- lm(VSTcnts ~ gene_n + nested.batch + oldLINEVST, data=dataratio)
      mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
    } else {
      mod0 <- lm(VSTcnts ~ garbage + oldLINEVST, data=dataratio)
      mod1 <- lm(VSTcnts ~ log2TPMsum + oldLINEVST, data=dataratio)
      mod2 <- lm(VSTcnts ~ gene_n + oldLINEVST, data=dataratio)
      mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3)
    }
    }
    
    bestmodel <- eval(getCall(AICc, 1))
    bestmodelnames <- rownames(AICc)[1]
    sm <- summary(bestmodel)
    print(sm)
    fixed_generatio <- any(rownames(sm$coefficients)=="gene_n")
    if (fixed_generatio) {
      plot(bestmodel$model[,"gene_n"], resid(bestmodel), ylab="residuals", xlab="gene_n")
      generatio_coef <- sm$coefficients["gene_n","Estimate"]
      generatio_pval <- sm$coefficients["gene_n","Pr(>|t|)"]
      if (genenames != "L1HS:L1:LINE" && !any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2 <- effectsize$Effects["gene_n","pEta-sqr"]
      } else {
        partialeta2 <- NA
      }
    }
    tested <- 1
  
  #}

  testedidx = which(tested==1)
  filename = paste(resultdir, subdir, cancerdir, "XBP1_cor.txt", sep="/")
  results <- cbind.data.frame(cancerdir, pearson_r, pearson_pval, spearman_rho, spearman_pval)
  colnames(results) <- c("tissue", "pearson_r" , "pearson_pval","spearman_rho", "spearman_pval")
  results <- results[order(-results[,"spearman_rho"]),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  filename = paste(resultdir, subdir, cancerdir, "XBP1_lm.txt", sep="/")
  results <- cbind.data.frame(cancerdir, fixed_generatio, bestmodelnames, generatio_coef, generatio_pval, partialeta2)
  colnames(results) <- c("tissue", "fixed_generatio", "bestmodelnames", "generatio_coef", "generatio_pval", "partialeta2")
  results <- results[order(-results$partialeta2),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
}




