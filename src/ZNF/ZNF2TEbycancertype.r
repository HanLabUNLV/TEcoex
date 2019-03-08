#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)

setwd("../../")
datadir = "data/"
resultdir = "result.rsem.TET/"
args = commandArgs(trailingOnly=TRUE)
TEfile = "result.rsem.TET/FullInfo/LINElist/L1HS:L1:LINE.txt"
print (args)
print (args[1])
if (length(args)==2) {
  resultdir = args[1]
  TEfile = args[2]
}
TEfamilyname <- basename(TEfile)
TEfamilyname <- substr(TEfamilyname,1,nchar(TEfamilyname)-4) 


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


TEcnts <- read.table(TEfile, sep="\t", header=TRUE)
colnames(TEcnts)[which(colnames(TEcnts)=="gene")] <-"VSTcnts"

log2TPMsum <- read.table(file=paste0(resultdir, "log2colsums.txt"), sep="\t", header = TRUE, row.names=1)
log2TPMsum <- cbind(substr(rownames(log2TPMsum), 1, 12), log2TPMsum)
colnames(log2TPMsum) <- c("patient", "log2TPMsum")

file.names <- read.table(file=paste(datadir, "ZNF", "ZNFlist.txt", sep="/"), header = FALSE, stringsAsFactors=FALSE)
file.names <- file.names[,1]
#genenames <- data.frame(matrix(unlist(strsplit(file.names, "[|]")), ncol=2, byrow=TRUE))[,1]
genenames <- substr(file.names, 1, nchar(file.names)-4)
stopifnot( genenames == unique(genenames))

cancertypes = list("BLCA", "BRCA", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), "LIHC", c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC")
cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))
dir.create(paste0(resultdir, "ZNF"))
newdirs = paste(resultdir, "ZNF", cancerdirs, sep="/")
sapply(newdirs, dir.create)

THCAradiation = read.table(file="data/THCA.radiation.txt", sep="\t", header=FALSE, row.names=1)
THCAradiation <- cbind.data.frame(substr(rownames(THCAradiation), 1, 12), THCAradiation)
colnames(THCAradiation) <- c("patient", "THCAradiation")


for (j in 1:length(cancertypes)) {
  cancertype <- cancertypes[[j]]
  cancerdir <- cancerdirs[j]
  filename = paste(resultdir, "ZNF", cancerdir, "results_positive_lm.txt", sep="/")
  if (file.exists(filename)) {
    print(paste0("skipping ", cancerdir))
    next
  }

  pearson_r <- rep(0, length(genenames))
  pearson_pval <- rep(1, length(genenames))
  spearman_rho <- rep(0, length(genenames))
  spearman_pval <- rep(1, length(genenames))
  bestmodelnames <- rep("mod", length(genenames))
  fixed_generatio <- rep(0, length(genenames))
  generatio_coef <- rep(0,length(genenames))
  generatio_pval <- rep(1,length(genenames))
  #R2marginal <- rep(0, length(genenames))
  #R2conditional <- rep(0, length(genenames))
  partialeta2 <- rep(0, length(genenames))
  tested <- rep(0, length(genenames))
  
  for(i in 1:length(genenames)){
    print (i)
    print (genenames[i])
    gene<- generatio(genenames[i])
    if (is.null(dim(gene))) next;

    data <- merge(TEcnts, gene, by = "patient")
    data <- merge(data, log2TPMsum, by = "patient")
    if (cancertype == "THCA") {
      data <- merge(data, THCAradiation, by = "patient")
    }
    
    subdata <- data[data$tissue %in% cancertype,]
    if (nrow(subdata)<8) next;
    if (all(colQuantiles(as.matrix(subdata$gene_n))["25%"]<log2(2))) next;
    subdata$tissue <- droplevels(subdata$tissue)
    subdata$BatchId <- factor(subdata$BatchId)
    if (length(cancertype)>1) {
      subdata$BatchId <- droplevels(subdata$BatchId)
      subdata$nested.batch=as.numeric(subdata$BatchId)
    }
    subdata$nested.batch = factor(subdata$nested.batch)
    dataratio <- subdata    
    
    test<-cor.test(dataratio$gene_n,dataratio$VSTcnts, method="pearson")
    pearson_r[i] <- test$estimate
    pearson_pval[i] <- test$p.value
    test<-cor.test(dataratio$gene_n,dataratio$VSTcnts, method="spearman")
    spearman_rho[i] <- test$estimate
    spearman_pval[i] <- test$p.value

    AICc = NULL    
    set.seed(0)
    garbage <- rnorm(length(dataratio$VSTcnts))  
    if (cancertype == "THCA") {
      if (length(summary(dataratio$nested.batch))>1) {
        mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
        mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
        mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
        mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
        mod4 <- lm(VSTcnts ~ nested.batch, data=dataratio)
        mod5 <- lm(VSTcnts ~ log2TPMsum + nested.batch, data=dataratio)
        mod6 <- lm(VSTcnts ~ gene_n + nested.batch, data=dataratio)
        mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch, data=dataratio)
        mod8 <- lm(VSTcnts ~ THCAradiation, data=dataratio)
        mod9 <- lm(VSTcnts ~ log2TPMsum + THCAradiation, data=dataratio)
        mod10 <- lm(VSTcnts ~ gene_n + THCAradiation, data=dataratio)
        mod11 <- lm(VSTcnts ~ gene_n + log2TPMsum + THCAradiation, data=dataratio)
        mod12 <- lm(VSTcnts ~ nested.batch + THCAradiation, data=dataratio)
        mod13 <- lm(VSTcnts ~ log2TPMsum + nested.batch + THCAradiation, data=dataratio)
        mod14 <- lm(VSTcnts ~ gene_n + nested.batch + THCAradiation, data=dataratio)
        mod15 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch + THCAradiation, data=dataratio)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
      } else {
        mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
        mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
        mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
        mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
        mod4 <- lm(VSTcnts ~ THCAradiation, data=dataratio)
        mod5 <- lm(VSTcnts ~ log2TPMsum + THCAradiation, data=dataratio)
        mod6 <- lm(VSTcnts ~ gene_n + THCAradiation, data=dataratio)
        mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + THCAradiation, data=dataratio)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      }      
    } else {
    if (length(summary(dataratio$nested.batch))>1) {
      mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
      mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
      mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
      mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
      mod4 <- lm(VSTcnts ~ nested.batch, data=dataratio)
      mod5 <- lm(VSTcnts ~ log2TPMsum + nested.batch, data=dataratio)
      mod6 <- lm(VSTcnts ~ gene_n + nested.batch, data=dataratio)
      mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
    } else {
      mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
      mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
      mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
      mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3)
    }
    }
    
    bestmodel <- eval(getCall(AICc, 1))
    bestmodelnames[i] <- rownames(AICc)[1]
    sm <- summary(bestmodel)
    fixed_generatio[i] <- any(rownames(sm$coefficients)=="gene_n")
    if (fixed_generatio[i]) {
      plot(bestmodel$model[,"gene_n"], resid(bestmodel), ylab="residuals", xlab="gene_n")
      generatio_coef[i] <- sm$coefficients["gene_n","Estimate"]
      generatio_pval[i] <- sm$coefficients["gene_n","Pr(>|t|)"]
      if (!any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2[i] <- effectsize$Effects["gene_n","pEta-sqr"]
      } else {
        partialeta2[i] <- NA
      }
    }
    tested[i] <- 1
    
  }

  testedidx = which(tested==1)
  if (length(testedidx) > 0) {
  pearson_pval[is.na(pearson_pval)] = 1
  pearson_qval = qvalue(pearson_pval[testedidx])$qvalues
  spearman_pval[is.na(spearman_pval)] = 1
  spearman_qval = qvalue(spearman_pval[testedidx])$qvalues
  filename = paste(resultdir, "ZNF/", cancerdir, "/", TEfamilyname,"_results_cor.txt", sep="")
  results <- cbind.data.frame(genenames[testedidx], pearson_r[testedidx], pearson_pval[testedidx], pearson_qval, spearman_rho[testedidx], spearman_pval[testedidx], spearman_qval)
  colnames(results) <- c("genenames", "pearson_r" , "pearson_pval", "pearson_qval" ,"spearman_rho", "spearman_pval", "spearman_qval" )
  results <- results[order(-results[,"spearman_rho"]),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  generatio_qval = qvalue(generatio_pval[testedidx])$qvalues 
  filename = paste(resultdir, "ZNF/", cancerdir, "/", TEfamilyname, "_results_lm.txt", sep="")
  results <- cbind.data.frame(genenames[testedidx], fixed_generatio[testedidx], bestmodelnames[testedidx], generatio_coef[testedidx], generatio_pval[testedidx], generatio_qval, partialeta2[testedidx])
  colnames(results) <- c("genenames", "fixed_generatio", "bestmodelnames", "generatio_coef", "generatio_pval", "generatio_qval", "partialeta2")
  results <- results[order(-results$partialeta2),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
  filename = paste(resultdir, "ZNF/",  cancerdir, "/", TEfamilyname, "_results_positive_lm.txt", sep="")
  positive<-as.data.frame(results[results$fixed_generatio==1,])
  positive[,4] <- as.numeric(as.character(positive[,4]))
  positive[,5] <- as.numeric(as.character(positive[,5]))
  positive[,6] <- as.numeric(as.character(positive[,6]))
  positive <- positive[order(-positive$partialeta2),]
  write.table(positive, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  }
}




coexpressed_plus = c()
coexpressed_minus = c()
for (j in 1:length(cancertypes)) {
  type <- cancertypes[j]
  cancerdir <- cancerdirs[j]
  fname = paste(resultdir, "ZNF/", cancerdir, "/", TEfamilyname, "_results_positive_lm.txt", sep="")
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
}
write.table(coexpressed_plus, file=paste0(resultdir, "ZNF/", TEfamilyname, "_coexpressed_plus.txt"), sep="\t", quote = FALSE, row.names = FALSE )
write.table(coexpressed_minus, file=paste0(resultdir, "ZNF/", TEfamilyname, "_coexpressed_minus.txt"), sep="\t", quote = FALSE, row.names = FALSE )



