#library("clusterProfiler")
library("pathview")
library("annotate")
library("org.Hs.eg.db")
#library(ReactomePA)



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
setwd(paste0("../", resultdir))

ranklist <- read.table("duplicate.rankedbyrsquared.rnk", header=FALSE)
rankminus <- ranklist[ranklist$V3 < 0,]
rankplus <- ranklist[ranklist$V3 > 0,]

#ranklist2 <- read.table("../results.oldLINE/duplicate.rankedbyrsquared.rnk", header=FALSE)
#rankminus2 <- ranklist2[ranklist2$V3 < 0,]
#rankplus2 <- ranklist2[ranklist2$V3 > 0,]
#
#x <- enrichPathway(gene=rankminus[,2], pvalueCutoff=0.04, readable=T)
#barplot(x, showCategory=15)
#dotplot(x, showCategory=15)
#x <- enrichPathway(gene=rankplus[,2], pvalueCutoff=0.04, readable=T)
#barplot(x, showCategory=15)
#dotplot(x, showCategory=15)
#
#rlist <- list(rankminus[,2], rankplus[,2], rankminus2[,2], rankplus2[,2])
#names(rlist) <- c("L1HS.negative", "L1HS.positive", "oldLINE.negative", "oldLINE.positive")
#
#res <- compareCluster(rlist, fun="enrichPathway")
#dotplot(res)

#en_david_KEGG <- enrichDAVID(gene = ranklist[,2], idType = "ENTREZ_GENE_ID", annotation="KEGG_PATHWAY", david.user = "mira.han@unlv.edu")
#x=compareCluster(rlist, fun="enrichDAVID", annotation="KEGG_PATHWAY", david.user="mira.han@unlv.edu")
#plot(x)

#en_david_GO <- enrichDAVID(gene = minuslist[,2], idType = "ENTREZ_GENE_ID", annotation=c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT"), david.user = "mira.han@unlv.edu")

#pdf("barplot.GO.minus.pdf")
#barplot(en_david_GO, drop=TRUE, showCategory=12)
#dev.off()

#pdf("dotplot.GO.minus.pdf")
#dotplot(en_david_GO)
#dev.off()

#pdf("enrichmap.GO.minus.pdf")
#enrichMap(en_david_GO)
#dev.off()


#pdf("GOterms.GO.minus.pdf")
#plotGOgraph(en_david_GO)
#dev.off()


genelist = ranklist[,2]
names(genelist) = ranklist[,1]
ids = select(org.Hs.eg.db, as.character(ranklist[,1]), "ENTREZID", "SYMBOL" )[,2]
names(genelist) = ids
hsa00190 <- pathview(gene.data  = genelist, pathway.id = "hsa00190", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa05012 <- pathview(gene.data  = genelist, pathway.id = "hsa05012", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa05016 <- pathview(gene.data  = genelist, pathway.id = "hsa05016", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa05010 <- pathview(gene.data  = genelist, pathway.id = "hsa05010", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa01100 <- pathview(gene.data  = genelist, pathway.id = "hsa01100", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa03010 <- pathview(gene.data  = genelist, pathway.id = "hsa03010", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa03040 <- pathview(gene.data  = genelist, pathway.id = "hsa03040", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04932 <- pathview(gene.data  = genelist, pathway.id = "hsa04932", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04210 <- pathview(gene.data  = genelist, pathway.id = "hsa04210", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04140 <- pathview(gene.data  = genelist, pathway.id = "hsa04140", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04141 <- pathview(gene.data  = genelist, pathway.id = "hsa04141", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04066 <- pathview(gene.data  = genelist, pathway.id = "hsa04066", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04150 <- pathview(gene.data  = genelist, pathway.id = "hsa04150", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04151 <- pathview(gene.data  = genelist, pathway.id = "hsa04151", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04064 <- pathview(gene.data  = genelist, pathway.id = "hsa04064", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04115 <- pathview(gene.data  = genelist, pathway.id = "hsa04115", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04144 <- pathview(gene.data  = genelist, pathway.id = "hsa04144", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
hsa04010 <- pathview(gene.data  = genelist, pathway.id = "hsa04010", species    = "hsa", low = list(gene = "red", cpd = "yellow"), high = list(gene = "green", cpd = "blue"))
