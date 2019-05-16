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
L1family = "ERV3-16A3_LTR:ERVL:LTR"
ZFPgene = "ZNF133|7692"
cancertype = "ESCASTAD"
if (length(args) > 0) {
  resultdir = args[1]
  L1family = args[2]
  ZFPgene = args[3]
  cancertype = args[4]
}
print (L1family)
print (ZFPgene)
print (cancertype)


cancerdir = cancertype
cancername = cancertype

if (cancertype == "CHOLLIHC") {
  cancertype <- c("CHOL", "LIHC")
} else if (cancertype == "COADREAD") {
  cancertype <- c("COAD", "READ")
} else if (cancertype == "ESCASTAD") {
  cancertype <- c("ESCA", "STAD")
} else if (cancertype == "KICHKIRCKIRP") {
  cancertype <- c("KICH", "KIRC", "KIRP")
} else if (cancertype == "LUADLUSC") {
  cancertype <- c("LUAD", "LUSC")
}





generatio <- function(genename)
{
  normalfname <- paste(generesultdir, "VSTcnts/", genename, ".txt", sep="")
  
  if (file.info(normalfname)$size == 0) {
    print ("file size zero")
    return (0)
  }
  gene <- read.table(normalfname, header=TRUE)
  colnames(gene) <- c("patient", "gene_n")
  return (gene)
}


readL1 <- function(L1dup, patientinfo)
{
  L1HS <- read.table(paste(resultdir, "VSTcnts/", L1dup, ":", L1family, ".txt", sep=""), header=TRUE)
  colnames(L1HS)[2] <-"VSTcnts"
  L1HS <- merge(L1HS, patientinfo, by="patient")
  return (L1HS)
}

patientinfo <- read.table(file=paste0(datadir, "patient.info"), sep="\t", header = TRUE)

log2TPMsum <- read.table(file=paste0(generesultdir, "log2colsums.txt"), sep="\t", header = TRUE, row.names=1)
log2TPMsum <- cbind(substr(rownames(log2TPMsum), 1, 12), log2TPMsum)
colnames(log2TPMsum) <- c("patient", "log2TPMsum")

#instead of using a list of genes, we will use a list of L1 dups and one gene
file.names <- dir(paste0(resultdir, "VSTcnts/"), pattern = paste0(L1family, ".txt"))
L1names <- substr(file.names, 1, nchar(file.names)-4)
L1names = unlist(strsplit(L1names, ":"))[c(TRUE, FALSE, FALSE, FALSE)]
stopifnot( L1names == unique(L1names))
names(file.names) <- L1names

#filter TEs uniquely mappable and no overlap with genes 1K
L1famname <- strsplit(L1family, ":")[[1]][1]
noovp <- read.table(file=paste0(resultdir,"gene.noovp.1000.uniq.TEs.bed"), sep="\t", header=TRUE)
noovpTEnames <- as.character(noovp[,4])
noovpTEnames <- cbind.data.frame(noovpTEnames, unlist(lapply(strsplit(noovpTEnames, "_"), "[[", 1)))
file.names.noovp <- as.character(noovpTEnames[noovpTEnames[,2]==L1famname,1])
file.names <- file.names[file.names.noovp]
L1names = names(file.names[!is.na(file.names)])
if (length(L1names) == 0) {
  exit
}

cancertypes = list("BLCA", "BRCA", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), "LIHC", c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC")
cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))
dir.create(paste0(resultdir, "ZNF"))
newdirs = paste(resultdir, "ZNF", cancerdirs, sep="/")
sapply(newdirs, dir.create)


THCAradiation = read.table(file="data/THCA.radiation.txt", sep="\t", header=FALSE, row.names=1)
THCAradiation <- cbind.data.frame(substr(rownames(THCAradiation), 1, 12), THCAradiation)
colnames(THCAradiation) <- c("patient", "THCAradiation")


pearson_r <- rep(0, length(L1names))
pearson_pval <- rep(1, length(L1names))
bestmodelnames <- rep("mod", length(L1names))
fixed_generatio <- rep(0, length(L1names))
generatio_coef <- rep(0,length(L1names))
generatio_pval <- rep(1,length(L1names))
#R2marginal <- rep(0, length(L1names))
#R2conditional <- rep(0, length(L1names))
partialeta2 <- rep(0, length(L1names))
tested <- rep(0, length(L1names))

for(i in 1:length(L1names)){
  #print (i)
  #print (L1names[i])
  #print (cancertype)

  filename = paste(resultdir, "ZNF/", cancerdir, "/", ZFPgene, "_", L1family, "_results_positive_lm.txt", sep="")
  if (file.exists(filename)) {
    print(paste0("skipping ", ZFPgene, "_", L1family))
    next
  }
  
  gene<- generatio(ZFPgene)
  L1HS<-readL1(L1names[i], patientinfo)
  if (is.null(dim(gene))) {
    print("null dimensions for gene")
    next;
  }

  data <- merge(L1HS, gene, by = "patient")
  data <- merge(data, log2TPMsum, by = "patient")
  if (cancertype == "THCA") {
    data <- merge(data, THCAradiation, by = "patient")
  }

  subdata <- data[data$tissue %in% cancertype,]
  if (all(sapply(subdata$VSTcnts, identical, subdata$VSTcnts[1]))) {
    #print("identical L1HS count")
    next;
  }

  if (nrow(subdata)<8) {
    #print("less than 8 rows")
    next;
  }
  if (all(colQuantiles(as.matrix(subdata$gene_n))["25%"]<log2(2))) {
    #print("col quantiles")
    next;
  }

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
pearson_qval = p.adjust(pearson_pval[testedidx])
filename = paste(resultdir, "ZNF/", cancerdir, "/", ZFPgene, "_", L1family, "_results_cor.txt", sep="")
results <- cbind.data.frame(L1names[testedidx], pearson_r[testedidx], pearson_pval[testedidx], pearson_qval)
colnames(results) <- c("L1names", "pearson_r" , "pearson_pval", "pearson_qval" )
results <- results[order(-results[,"pearson_r"]),]
write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)


generatio_qval = p.adjust(generatio_pval[testedidx])
filename = paste(resultdir, "ZNF/", cancerdir, "/", ZFPgene, "_", L1family, "_results_lm.txt", sep="")
results <- cbind.data.frame(L1names[testedidx], fixed_generatio[testedidx], bestmodelnames[testedidx], generatio_coef[testedidx], generatio_pval[testedidx], generatio_qval, partialeta2[testedidx])
colnames(results) <- c("L1names", "fixed_generatio", "bestmodelnames", "generatio_coef", "generatio_pval", "generatio_qval", "partialeta2")
results <- results[order(-results$partialeta2),]
write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)


filename = paste(resultdir, "ZNF/",  cancerdir, "/", ZFPgene, "_", L1family, "_results_positive_lm.txt", sep="")
positive<-as.data.frame(results[results$fixed_generatio==1,])
positive[,4] <- as.numeric(as.character(positive[,4]))
positive[,5] <- as.numeric(as.character(positive[,5]))
positive[,6] <- as.numeric(as.character(positive[,6]))
positive <- positive[order(-positive$partialeta2),]
write.table(positive, file=filename, sep="\t", row.names = FALSE, quote=FALSE)

}





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

