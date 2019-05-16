library("RColorBrewer")
library("DESeq2")
library("dendsort")
library("ComplexHeatmap")
library("circlize")
source("nmi.r")


setwd("../")
args = commandArgs(trailingOnly=TRUE)
resultdir = "result.rsem.TET.instance/"
print (args)
if (length(args)>0) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
  b_noovp = args[2]
}
print(resultdir)
setwd(resultdir)
print(b_noovp)

load(file = "ddsNew.RData")
load(file = "VSTcnts.DEGES.RData")

# filter ESCA STAD
cdata <- colData(ddsnew)
sampleidx = rep(TRUE, ncol(VSTcnts))
ESCASTADidx = (cdata$tissue == "ESCA") | (cdata$tissue == "STAD")
sampleidx = !ESCASTADidx
VSTcnts = VSTcnts[,sampleidx]
coldataNew = colData(ddsnew)
coldataNew = coldataNew[sampleidx,]

print(dim(VSTcnts))
TEs <- grepl(":", rownames(VSTcnts))
geneVSTcnts <- VSTcnts[!TEs,]
TEcnts <- VSTcnts[TEs,]
if (grepl("instance", resultdir) && (! is.na(b_noovp)) && (b_noovp == "noovp")) {
  filters = read.table("gene.noovp.100000.TEs.bed", sep="\t")
  filternames = as.character(filters[,4])
  print(dim(TEcnts))
  TEnames <- rownames(TEcnts)
  TElocus = matrix(unlist(strsplit(TEnames, ":")), ncol=5, byrow=TRUE)[,1]
  rownames(TEcnts) = TElocus
  names(TEnames) = TElocus
  filteredTEcnts <- TEcnts[filternames,]
  print(dim(filteredTEcnts))
  filteredTEnames <- TEnames[filternames]
  rownames(filteredTEcnts) = filteredTEnames
  #rm(VSTcnts)
  #VSTcnts = rbind(geneVSTcnts, filteredTEcnts)
} else {
  filteredTEcnts = TEcnts
}



lineVSTcnts <- filteredTEcnts[grepl("LINE", rownames(filteredTEcnts)), ]
print(dim(lineVSTcnts))
sineVSTcnts <- filteredTEcnts[grepl("SINE", rownames(filteredTEcnts)), ]
print(dim(sineVSTcnts))
dnaVSTcnts <- filteredTEcnts[grepl("DNA", rownames(filteredTEcnts)), ]
print(dim(dnaVSTcnts))
ltrVSTcnts <- filteredTEcnts[grepl("LTR", rownames(filteredTEcnts)), ]
print(dim(ltrVSTcnts))
hervVSTcnts <- filteredTEcnts[grepl("HERV", rownames(filteredTEcnts)), ]
print(dim(hervVSTcnts))
youngL1idx <- grepl("L1HS", rownames(filteredTEcnts))
youngL1idx <- youngL1idx | grepl("L1PA2_", rownames(filteredTEcnts))
youngL1idx <- youngL1idx | grepl("L1PA3_", rownames(filteredTEcnts))
youngL1VSTcnts <- filteredTEcnts[youngL1idx, ]
print(dim(youngL1VSTcnts))
print (rownames(youngL1VSTcnts))

# annotation
annotation_data <- as.data.frame(coldataNew$tissue)
rownames(annotation_data) <- colnames(geneVSTcnts)
colnames(annotation_data) <- "tissue"

broadertype = as.character(annotation_data$tissue)
broadertype[broadertype=="BLCA"] = "BLCAUCEC"
broadertype[broadertype=="UCEC"] = "BLCAUCEC"
broadertype[broadertype=="COAD"] = "COADREAD"
broadertype[broadertype=="READ"] = "COADREAD"
broadertype[broadertype=="ESCA"] = "ESCASTAD"
broadertype[broadertype=="STAD"] = "ESCASTAD"
broadertype[broadertype=="KICH"] = "KICHKIRCKIRP"
broadertype[broadertype=="KIRC"] = "KICHKIRCKIRP"
broadertype[broadertype=="KIRP"] = "KICHKIRCKIRP"
broadertype[broadertype=="LUAD"] = "LUADLUSC"
broadertype[broadertype=="LUSC"] = "LUADLUSC"
annotation_data <- cbind.data.frame(annotation_data,broadertype)

#callback = function(hc, ...){dendsort(hc)}


#colors
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28",  "#800000", "#fffac8", "#aaffc3",  "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
col_vector2 = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#d2f53c", "#fabebe", "#e6beff", "#800000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")

mycolors=list(tissue = col_vector[1:length(levels(annotation_data$tissue))], 
              broadertype = col_vector2[1:length(levels(annotation_data$broadertype))])
names(mycolors$tissue) <- levels(annotation_data$tissue) 
names(mycolors$broadertype) <- levels(annotation_data$broadertype) 


rnMI <- rep(0, 100)
for (i in 100) {
  randomtype1 <- sample(annotation_data$broadertype, length(annotation_data$broadertype), replace = FALSE)
  randomtype2 <- sample(annotation_data$broadertype, length(annotation_data$broadertype), replace = FALSE)
  nMI <- normalizedMI(randomtype1, randomtype2)
  rnMI[i] <- nMI 
}
print("nMI random")
print(mean(rnMI))



  #heatmap
  plotheatmap <- function(VSTcnts, filename, top) {

    varGenes <- rowVars(VSTcnts)
    topVarianceGenes <- head(order(varGenes, decreasing=T),top)
    mat <- VSTcnts[ topVarianceGenes, ]
    mat <- mat - rowMeans(mat)
    for (run in 1:10) {
      hc=hclust(d = dist(t(mat)), method = "ave")
      memb <- cutree(hc, k=length(levels(annotation_data$tissue)))

      cid <- as.data.frame(sort(memb))
      clusterid <- merge(annotation_data, cid, by="row.names")
      colnames(clusterid)[4] = "clusterid"
      cluster2gtype = rep(0, nrow(clusterid))
      clusterid$clusterid <- factor(clusterid$clusterid)

      for (i in 1:length(levels(clusterid$clusterid))) {
        gtype <- summary(clusterid[clusterid$clusterid==i,"broadertype"])
        typeidx = which(gtype==max(gtype))
        if (length(typeidx) > 1) {
          gtype <- summary(clusterid[clusterid$clusterid==i,"broadertype"])/summary(clusterid[,"broadertype"])
          typeidx = which(gtype==max(gtype))
        }
        cluster2gtype[clusterid$clusterid==i] = names(typeidx) # assign most frequent broadertype to the cluster
      }
      clusterid <- cbind.data.frame(clusterid, cluster2gtype)
      clusterid$cluster2gtype <- factor(clusterid$cluster2gtype, levels=levels(clusterid$broadertype))
      nMI <- normalizedMI(clusterid$broadertype, clusterid$cluster2gtype)
      print(paste0("normalized Mutual Information ", filename, ": "))
      print(nMI)

    }

    dend_col = dendsort(as.dendrogram(hc))
    dend_row = dendsort(as.dendrogram(hclust(dist(mat))))
    # select the 'contrast' you want
    pdf(file=filename, width=20, height=12)
    par(ps=3)


    ha = HeatmapAnnotation(annotation_data, col=mycolors, 
        #annotation_legend_param = list(type = list(grid_height=10)),
        )
    #pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
    ht = Heatmap(mat, cluster_rows = dend_row, row_dend_reorder = FALSE, 
              name = "ht", cluster_columns = dend_col, column_dend_reorder = FALSE, 
              top_annotation = ha,
              top_annotation_height = unit(8, "cm"),
              show_row_names = FALSE, show_column_names = FALSE,
              #col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
              col = colorRamp2(c(-5, 0, 5), c("blue", "#ffffbf", "red")),
              column_title = "")
    #upViewport(2)
    draw(ht, newpage = FALSE)
    dev.off()

  }


  plotheatmap(geneVSTcnts, paste0("geneheatmapTop150Var", b_noovp, "100k.pdf"), 150)
  plotheatmap(lineVSTcnts, paste0("LINEheatmapTop150Var", b_noovp, "100k.pdf"), 150)
  plotheatmap(sineVSTcnts, paste0("SINEheatmapTop150Var", b_noovp, "100k.pdf"), 150)
  plotheatmap(dnaVSTcnts, paste0("DNAheatmapTop150Var", b_noovp, "100k.pdf"), 150)
  plotheatmap(ltrVSTcnts, paste0("LTRheatmapTop150Var", b_noovp, "100k.pdf"), 150)
  plotheatmap(hervVSTcnts, paste0("HERVheatmapTop150Var", b_noovp, "100k.pdf"), 150)
  plotheatmap(youngL1VSTcnts, paste0("youngL1heatmapTop150Var", b_noovp, "100k.pdf"), 150)

