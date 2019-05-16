
library("RDAVIDWebService")
library(org.Hs.eg.db)
library(annotate)
library("ReactomePA")


getClusterSymbols <- function (termCluster) {
  N= nrow(summary(termCluster))
  clustersymbols <- vector("list", N)
  for (i in 1:N ) {
    clusterTerms<-members(termCluster)[[i]]
    genes <- unique(unlist(strsplit(clusterTerms$Genes, ", ")))
    symbols <- getSYMBOL(genes, data='org.Hs.eg')
    clustersymbols[[i]] <- symbols
  }
  return (clustersymbols)
}


runDavid <- function(resultdir) {
  david<-DAVIDWebService(email="mira.han@unlv.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setAnnotationCategories(david, getDefaultCategoryNames(david))
  minuslist <- read.table(paste(resultdir,"duplicate.minus.rsquared",sep="/"), header=FALSE)
  ids = minuslist[,2]
  if(grepl('TET', resultdir)) {
    print (minuslist[,1])
    ids = select(org.Hs.eg.db, as.character(minuslist[,1]), "ENTREZID", "SYMBOL")[,2]
  }
    
  result<-addList(david, ids, idType="ENTREZ_GENE_ID", listName="minusduplicates", listType="Gene")
  termClusterMinus<-getClusterReport(david, type="Term")
  symbolfile = "termClusterSymbols.minus.tab"
  if (file.exists(symbolfile)) file.remove(symbolfile)
  if (length(members(termClusterMinus)) > 0) {
    symbolsMinus <- getClusterSymbols(termClusterMinus)
    symbolsMinus <- lapply(symbolsMinus, sort)
    lapply(symbolsMinus, write, paste(resultdir,symbolfile, sep="/"), append=TRUE, ncolumns=1000)
    getClusterReportFile(david, type="Term", fileName=paste(resultdir,"termClusterReport.minus.tab",sep="/"))
    getFunctionalAnnotationChartFile(david, fileName=paste(resultdir,"FAChart.minus.tab",sep="/"))
  }

  #davidGODag<-DAVIDGODag(members(termClusterMinus)[[1]],pvalueCutoff=0.1, "CC")
  #plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")

  pluslist <- read.table(paste(resultdir,"duplicate.plus.rsquared",sep="/"), header=FALSE)
  ids = pluslist[,2]
  if(grepl('TET', resultdir)) {
    print (pluslist[,1])
    ids = select(org.Hs.eg.db, as.character(pluslist[,1]), "ENTREZID", "SYMBOL")[,2]
  }
  result<-addList(david, ids, idType="ENTREZ_GENE_ID", listName="plusduplicates", listType="Gene")
  termClusterPlus<-getClusterReport(david, type="Term")
  symbolfile = "termClusterSymbols.plus.tab"
  if (file.exists(symbolfile)) file.remove(symbolfile)
  if (length(members(termClusterPlus)) > 0) { 
    symbolsPlus <- getClusterSymbols(termClusterPlus)
    symbolsPlus <- lapply(symbolsPlus, sort)
    lapply(symbolsPlus, write, paste(resultdir,symbolfile, sep="/"), append=TRUE, ncolumns=1000)
    getClusterReportFile(david, type="Term", fileName=paste(resultdir,"termClusterReport.plus.tab",sep="/"))
    getFunctionalAnnotationChartFile(david, fileName=paste(resultdir,"FAChart.plus.tab",sep="/"))
  }

}



setwd("../")
args = commandArgs(trailingOnly=TRUE)
resultdir = "result.rsem/gene2L1HS"
print (args)
if (length(args)>0) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
}
print(resultdir)
#runDavid(resultdir)


  pluslist <- read.table(paste(resultdir,"duplicate.plus.rsquared",sep="/"), header=FALSE)
  ids = pluslist[,2]
  x <- enrichPathway(gene=ids, pvalueCutoff=0.05, readable=T)
  if (! is.null(x) && nrow(x)) {
    write.table(as.data.frame(x), file = paste0(resultdir, paste("Reactome.plus.txt", sep=".") ), row.names = FALSE, quote = FALSE);
    dotplot(x, showCategory=15)
    ggsave(filename=paste0(resultdir, paste("Reactome.plus.pdf", sep=".")), width=16, height=(min(nrow(as.data.frame(x)),15)*0.32), units = "in", limitsize = FALSE )
    emapplot(x)
    ggsave(filename=paste0(resultdir, paste("Reactome.plus.em.pdf", sep=".")), width=8, height=12, units = "in")
  }
  minuslist <- read.table(paste(resultdir,"duplicate.minus.rsquared",sep="/"), header=FALSE)
  ids = minuslist[,2]
  x <- enrichPathway(gene=ids, pvalueCutoff=0.05, readable=T)
  if (! is.null(x) && nrow(x)) {
    write.table(as.data.frame(x), file = paste0(resultdir, paste("Reactome.minus.txt", sep=".") ), row.names = FALSE, quote = FALSE);
    dotplot(x, showCategory=15)
    ggsave(filename=paste0(resultdir, paste("Reactome.minus.pdf", sep=".")), width=16, height=(min(nrow(as.data.frame(x)),15)*0.32), units = "in", limitsize = FALSE )
    emapplot(x)
    ggsave(filename=paste0(resultdir, paste("Reactome.minus.em.pdf", sep=".")), width=8, height=12, units = "in")
  }

