require(ggplot2)
library("scales")
#require(reshape2)




integer_breaks <- function(n = 3, ...) {
  breaker <- pretty_breaks(n, ...)
  function(x) {
    breaks <- breaker(x)
    newbreaks = breaks[breaks == floor(breaks)]
    if((length(newbreaks) > 3) || (breaks[length(breaks)]  != newbreaks[length(newbreaks)] )) {
      breaks = newbreaks
    }
    breaks
  }
}



generatio <- function(genename)
{
  normalfname <- paste("VSTcnts/", genename, ".txt", sep="")
  
  if (file.info(normalfname)$size == 0) {
    print ("file size zero")
    return (0)
  }
  gene <- read.table(normalfname, header=TRUE)
  colnames(gene) <- c("patient", "gene_n")
  return (gene)
}


gene2L1HS <- function(genename, cancertypes)
{
  gene<- generatio(genename)
  if (is.null(dim(gene))) next;
  
  data <- merge(L1HS, gene, by = "patient")
  
  idx <- lapply(data$broadertype, grepl, cancertypes)
  idx <- lapply(idx, any)
  idx <- unlist(idx)
  subdata <- data[idx,]
  
  return(subdata)
  
  
  
}


plotcorrelation <- function(resultdir, direction, genelist = NULL)
{

  coex <- read.table(paste0(resultdir, "coexpressed_", direction, ".txt"), header=TRUE, sep="\t")
  if (is.null(genelist)) {
    coex_dup <- read.table(paste0(resultdir, "coexpressed_", direction, ".duplicates.txt"))
    coex = cbind(coex_dup[,2], coex_dup)
    coex = coex[,-3]
    coex_g <- split(coex, coex[,1])
    
    eta2_g <- lapply(coex_g, "[",8)
    bar <- function(x) {return(sum(x[,1]))}
    med_eta2s <- unlist(lapply(eta2_g, bar))
    med_eta2s <- med_eta2s[order(-med_eta2s)]
    genelist = med_eta2s
  } 
  
  for (i in 1:12) {
    genename <- names(genelist[i])
    genesymbol = strsplit(genename, "[|]")[[1]][1]
    cancertypes=coex[coex[,2]==genename,1]
    if (length(cancertypes) == 0) next;
    subdata <- gene2L1HS(genename, cancertypes)
    
    
    npanels = length(levels(factor(subdata$broadertype)))
    
    #ggplot(subdata, aes(x=gene_n, y=VSTcnts, colour=broadertype) ) + geom_point() + scale_colour_manual(values = col_vector2) + theme_bw()
    ggplot(subdata) +
      geom_jitter(aes(x=gene_n, y=VSTcnts, colour=tissue)) + 
      geom_smooth(aes(x=gene_n, y=VSTcnts), method=lm, se=FALSE, color = "darkgrey") +
      facet_wrap(~broadertype, scales="free_x", nrow=1) +
      scale_x_continuous(breaks = integer_breaks()) +
      labs(x =  paste0("log2 ", genesymbol), y = "log2 L1HS 5'") +
      scale_colour_manual(values = col_vector[levels(factor(subdata$tissue))]) + theme_bw()
    
    ggsave(filename=paste0(resultdir, direction, ".", genesymbol, "2L1HS.pdf"), width = 2*npanels, height = 2, units = c("in"), dpi = 300)
    
    
    ggplot(subdata) +
      geom_jitter(aes(x=gene_n, y=oldLINEVST, colour=tissue)) + 
      geom_smooth(aes(x=gene_n, y=oldLINEVST), method=lm, se=FALSE, color = "darkgrey") +
      facet_wrap(~broadertype, scales="free_x", nrow=1) +
      scale_x_continuous(breaks = integer_breaks()) +
      labs(x =  paste0("log2 ", genesymbol), y = "log2 oldLINEs") +
      scale_colour_manual(values = col_vector[levels(factor(subdata$tissue))]) + theme_bw()
    
    ggsave(filename=paste0(resultdir, direction, ".", "controls", genesymbol, "2oldLINEs.pdf"), width = 2*npanels, height = 2, units = c("in"), dpi = 300)
    
    
  } 
}





rundir = "../result.rsem.TET/"
setwd(rundir)

L1HS <- read.table("L1HS5prime.VST.txt", sep="\t", header=TRUE, row.names=1)

oldLINEVST <- read.table(file="oldLINE.VST.txt", sep="\t", header = TRUE, row.names=1)
colnames(oldLINEVST) <- c("patient", "oldLINEVST")

broadertype = as.character(L1HS$tissue)
broadertype[broadertype=="COAD"] = "COADREAD"
broadertype[broadertype=="READ"] = "COADREAD"
broadertype[broadertype=="ESCA"] = "ESCASTAD"
broadertype[broadertype=="STAD"] = "ESCASTAD"
broadertype[broadertype=="KICH"] = "KICHKIRCKIRP"
broadertype[broadertype=="KIRC"] = "KICHKIRCKIRP"
broadertype[broadertype=="KIRP"] = "KICHKIRCKIRP"
broadertype[broadertype=="LUAD"] = "LUADLUSC"
broadertype[broadertype=="LUSC"] = "LUADLUSC"
L1HS <- cbind.data.frame(L1HS,broadertype)
L1HS <- merge(L1HS, oldLINEVST, by = "patient")





col_vector = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28",  "#800000", "#fffac8", "#aaffc3",  "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
col_vector2 = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#d2f53c", "#fabebe", "#e6beff", "#800000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")

names(col_vector) <- levels(factor(L1HS$tissue))
names(col_vector2) <- levels(factor(L1HS$broadertype))

resultdir = paste0(rundir, "gene2L1HS5prime/")


genelist <- read.table(paste0(resultdir, "REC.minus.new.txt"), header=TRUE)
genes <-  genelist[1:100,2]
names(genes) <-  genelist[1:100,1]
genes <- genes[!grepl(":", names(genes))]
plotcorrelation(resultdir, "minus", genes)
genelist <- read.table(paste0(resultdir, "REC.plus.new.txt"), header=TRUE)
genes <-  genelist[1:100,2]
names(genes) <-  genelist[1:100,1]
genes <- genes[!grepl(":", names(genes))]
plotcorrelation(resultdir, "plus", genes)

plotcorrelation(resultdir, "minus")
plotcorrelation(resultdir, "plus")
#resultdir = paste0(rundir, "results.VSTnormal/")
#plotcorrelation(resultdir, "minus")
#plotcorrelation(resultdir, "plus")

