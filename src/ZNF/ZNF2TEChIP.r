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
#L1family = "THE1B-int:ERVL-MaLR:LTR"
#ZFPgene = "ZNF133|7692"
ChIPlistfile = "/data2/han_lab/adrian/ZFP_TE/results.maxgt5.txt"
#ChIPlistfile = "/data2/han_lab/mobiledna/data/TET/random.maxgt5.unbound.txt"
cancertype = "BLCA"
if (length(args) > 0) {
  resultdir = args[1]
#  L1family = args[2]
#  ZFPgene = args[3]
#  cancertype = args[4]
  ChIPlistfile = args[2]
  cancertype = args[3]
}
#print (L1family)
#print (ZFPgene)
print (resultdir)
print (ChIPlistfile)
ChIPbasename = basename(ChIPlistfile)
print (ChIPbasename)
print (cancertype)


cancerdir = cancertype
cancername = cancertype

if (cancertype == "COADREAD") {
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
  genename = sub("-rep2", "", genename)
  genename = sub("-rep3", "", genename)

  normalfname <- list.files(paste0(generesultdir, "VSTcnts/"), pattern=paste0("^",genename, "\\|[0-9]+\\.txt$"))
  if (length(normalfname) < 1) {
    print ("gene file not found")
    return (0)
  }
  if (length(normalfname) > 1) {
    print ("multiple files found")
    return (0)
  }
  normalfname <- paste0(generesultdir, "VSTcnts/", normalfname)
  if (file.info(normalfname[1])$size == 0) {
    print ("file size zero")
    return (0)
  }
  gene <- read.table(normalfname, header=TRUE)
  colnames(gene) <- c("patient", "gene_n")
  return (gene)
}


readL1 <- function(L1filename, patientinfo)
{
  l1file = paste0(resultdir, "VSTcnts/", L1filename, ".txt")
  print(l1file)
  
  if (!file.exists(l1file)) {
    print ("TE file not found")
    return (0)
  }
  if (file.info(l1file)$size == 0) {
    print ("file size zero")
    return (0)
  }
  L1HS <- read.table(l1file, header=TRUE)
  colnames(L1HS)[2] <-"VSTcnts"
  L1HS <- merge(L1HS, patientinfo, by="patient")
  return (L1HS)
}

patientinfo <- read.table(file=paste0(datadir, "patient.info"), sep="\t", header = TRUE)

log2TPMsum <- read.table(file=paste0(generesultdir, "log2colsums.txt"), sep="\t", header = TRUE, row.names=1)
log2TPMsum <- cbind(substr(rownames(log2TPMsum), 1, 12), log2TPMsum)
colnames(log2TPMsum) <- c("patient", "log2TPMsum")

oldLINEVST <- read.table(file=paste0(generesultdir, "/WGCNA/consensus/eigengeneME1.txt"), sep="\t", header = FALSE)
colnames(oldLINEVST) <- c("patient", "oldLINEVST")

#instead of using a list of genes, we will use a list of L1 dups and one gene
#file.names <- dir(paste0(resultdir, "VSTcnts/"), pattern = paste0(L1family, ".txt"))
ChIPlist <- read.table(ChIPlistfile, header=FALSE, sep="\t", stringsAsFactor=FALSE)
#L1names <- substr(file.names, 1, nchar(file.names)-4)
L1names <- ChIPlist[,1]
ZFPgenes <- ChIPlist[,2]
L1filenames <- ChIPlist[,3]
#L1names = unlist(strsplit(L1names, ":"))[c(TRUE, FALSE, FALSE, FALSE)]
#stopifnot( L1names == unique(L1names))

missing_TEs = c()

cancertypes = list("BLCA", "BRCA", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), "LIHC", c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC")
cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))
dir.create(paste0(resultdir, "ZNF.ctrl.unbound"))
newdirs = paste(resultdir, "ZNF.ctrl.unbound", cancerdirs, sep="/")
sapply(newdirs, dir.create)


THCAradiation = read.table(file="data/THCA.radiation.txt", sep="\t", header=FALSE, row.names=1)
THCAradiation <- cbind.data.frame(substr(rownames(THCAradiation), 1, 12), THCAradiation)
colnames(THCAradiation) <- c("patient", "THCAradiation")


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

  filename = paste0(resultdir, "ZNF.ctrl.unbound/", cancerdir, "/", ChIPbasename, "_results_positive_lm.txt")
  if (file.exists(filename)) {
    print(paste0("skipping ", ChIPbasename))
    next
  }
  
  gene<- generatio(ZFPgenes[i])
  L1HS<-readL1(L1filenames[i], patientinfo)
  if (is.null(dim(gene))) {
    print("null dimensions for gene")
    next;
  }
  if (is.null(dim(L1HS))) {
    print("null dimensions for TE")
    missing_TEs = rbind(missing_TEs, c(L1names[i], ZFPgenes[i]))
    next;
  }

  data <- merge(L1HS, gene, by = "patient")
  data <- merge(data, log2TPMsum, by = "patient")
  data <- merge(data, oldLINEVST, by = "patient")
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
 
#
#
#  AICc = NULL
#  set.seed(0)
#  garbage <- rnorm(length(dataratio$VSTcnts))
#  if (cancertype == "THCA") {
#    if (length(summary(dataratio$nested.batch))>1) {
#      mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
#      mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
#      mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
#      mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
#      mod4 <- lm(VSTcnts ~ nested.batch, data=dataratio)
#      mod5 <- lm(VSTcnts ~ log2TPMsum + nested.batch, data=dataratio)
#      mod6 <- lm(VSTcnts ~ gene_n + nested.batch, data=dataratio)
#      mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch, data=dataratio)
#      mod8 <- lm(VSTcnts ~ THCAradiation, data=dataratio)
#      mod9 <- lm(VSTcnts ~ log2TPMsum + THCAradiation, data=dataratio)
#      mod10 <- lm(VSTcnts ~ gene_n + THCAradiation, data=dataratio)
#      mod11 <- lm(VSTcnts ~ gene_n + log2TPMsum + THCAradiation, data=dataratio)
#      mod12 <- lm(VSTcnts ~ nested.batch + THCAradiation, data=dataratio)
#      mod13 <- lm(VSTcnts ~ log2TPMsum + nested.batch + THCAradiation, data=dataratio)
#      mod14 <- lm(VSTcnts ~ gene_n + nested.batch + THCAradiation, data=dataratio)
#      mod15 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch + THCAradiation, data=dataratio)
#      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
#    } else {
#      mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
#      mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
#      mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
#      mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
#      mod4 <- lm(VSTcnts ~ THCAradiation, data=dataratio)
#      mod5 <- lm(VSTcnts ~ log2TPMsum + THCAradiation, data=dataratio)
#      mod6 <- lm(VSTcnts ~ gene_n + THCAradiation, data=dataratio)
#      mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + THCAradiation, data=dataratio)
#      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
#    }
#  } else {
#  if (length(summary(dataratio$nested.batch))>1) {
#    mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
#    mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
#    mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
#    mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
#    mod4 <- lm(VSTcnts ~ nested.batch, data=dataratio)
#    mod5 <- lm(VSTcnts ~ log2TPMsum + nested.batch, data=dataratio)
#    mod6 <- lm(VSTcnts ~ gene_n + nested.batch, data=dataratio)
#    mod7 <- lm(VSTcnts ~ gene_n + log2TPMsum + nested.batch, data=dataratio)
#    AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
#  } else {
#    mod0 <- lm(VSTcnts ~ garbage, data=dataratio)
#    mod1 <- lm(VSTcnts ~ log2TPMsum, data=dataratio)
#    mod2 <- lm(VSTcnts ~ gene_n, data=dataratio)
#    mod3 <- lm(VSTcnts ~ gene_n + log2TPMsum, data=dataratio)
#    AICc<-model.sel(mod0,mod1,mod2,mod3)
#  }
#  }
#

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


write.table(missing_TEs, file=paste0(resultdir, "ZNF.ctrl.unbound/", cancerdir, "/", ChIPbasename, "_missing_TEs.txt"), sep="\t", quote=FALSE)

testedidx = which(tested==1)
if (length(testedidx) > 0) {


generatio_qval = p.adjust(generatio_pval[testedidx])
filename = paste0(resultdir, "ZNF.ctrl.unbound/", cancerdir, "/", ChIPbasename,"_results_lm.txt")
results <- cbind.data.frame(L1names[testedidx], fixed_generatio[testedidx], bestmodelnames[testedidx], generatio_coef[testedidx], generatio_pval[testedidx], generatio_qval, partialeta2[testedidx])
colnames(results) <- c("L1names", "fixed_generatio", "bestmodelnames", "generatio_coef", "generatio_pval", "generatio_qval", "partialeta2")
results <- results[order(-results$partialeta2),]
write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)


filename = paste0(resultdir, "ZNF.ctrl.unbound/", cancerdir, "/", ChIPbasename,  "_results_positive_lm.txt")
positive<-as.data.frame(results[results$fixed_generatio==1,])
positive[,4] <- as.numeric(as.character(positive[,4]))
positive[,5] <- as.numeric(as.character(positive[,5]))
positive[,6] <- as.numeric(as.character(positive[,6]))
positive <- positive[order(-positive$partialeta2),]
write.table(positive, file=filename, sep="\t", row.names = FALSE, quote=FALSE)

}





coexpressed_plus = c()
coexpressed_minus = c()
fname = paste(resultdir, "ZNF.ctrl.unbound/", cancerdir, "/", ChIPbasename, "_results_positive_lm.txt", sep="")
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
  filename = paste0(resultdir, "ZNF.ctrl.unbound/", cancerdir, "/", ChIPbasename,  "_coexpressed_plus.txt")
  write.table(coexpressed_plus, file=filename, sep="\t", quote = FALSE, row.names = FALSE )
}
if (nrow(coexpressed_minus)>0) {
  filename = paste0(resultdir, "ZNF.ctrl.unbound/", cancerdir, "/", ChIPbasename,  "_coexpressed_minus.txt")
  write.table(coexpressed_minus, file=filename, sep="\t", quote = FALSE, row.names = FALSE )
}

