#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


library(WGCNA);
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "..";
setwd(workingDir); 
 
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#enableWGCNAThreads(nThreads = 2)
args = commandArgs(trailingOnly=TRUE)
datadir = "data/"
resultdir = "result.rsem.TET.instance/"
print (args)
if (length(args)>0) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
  VSTcntfile = args[2]
  setname = unlist(strsplit(VSTcntfile, "[.]"))[2]
  #setname = args[2]
  shortname = args[3]
  powerparam = as.numeric(args[4])
}
print(resultdir)
print(setname)
print(shortname)
print(powerparam)

dir.create(paste0(resultdir, "WGCNA")) 
dir.create(paste0(resultdir, "WGCNA/bytissue")) 
dir.create(paste0(resultdir, "WGCNA/bytissue/", shortname)) 
dir.create(paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots")) 



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

traitData = read.table(paste0(resultdir, "WGCNA/", "WGCNA.L1HS.clinical.txt"), header=TRUE, sep="\t", stringsAsFactors=TRUE);




  datExpr = NULL 
  datTraits = NULL
  inputfile = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("WGCNA.dataInput.",shortname,".RData"))
  if (file.exists(inputfile)) {
    load(inputfile)
  } else {

  #normalcnts <- read.table(paste0(resultdir, "VSTcnt.txt"), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
  normalcnts <- read.table(paste0(resultdir, VSTcntfile), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
  dim(normalcnts);
  names(normalcnts) = substr(names(normalcnts), 1, 12);

  samples = traitData[traitData$newtype==setname,]$patient
  datExpr0 = as.data.frame(t(normalcnts[,as.character(samples)]))
  names(datExpr0) = rownames(normalcnts);
  rownames(datExpr0) = samples;

  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK

  if (!gsg$allOK)
  {
  # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
       printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
       printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  write.table(datExpr0, file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("WGCNA.genes.",shortname,".txt")), quote=FALSE, row.names=TRUE, sep="\t")
#
#sampleTree = hclust(dist(datExpr0), method = "average");
## Plot the sample tree: Open a graphic output window of size 12 by 9 inches
## The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
##pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
#     cex.axis = 1.5, cex.main = 2)
#
## Plot a line to show the cut
#abline(h = 15, col = "red");
## Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
## clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
  datExpr = datExpr0
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)


  #traitData = read.table(paste0(resultdir, "WGCNA/", "WGCNA.L1HS.clinical.txt"), header=TRUE, sep="\t", stringsAsFactors=TRUE);
  #dim(traitData)
  #names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.

  Samples = rownames(datExpr);
  traitRows = match(Samples, traitData$patient);
  allTraits = traitData[traitRows,]
  rownames(allTraits) = traitData[traitRows, 1];
  datTraits = allTraits[,c("L1HS", "oldLINE")];

  write.table(allTraits, file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("WGCNA.trait.",shortname,".txt")), quote=FALSE, row.names=TRUE, sep="\t")
  collectGarbage();

# Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
  pdf(file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", "sample.clustering.",shortname,".pdf"), width = 12, height = 9);
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits), 
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
 
  save(datExpr, datTraits, allTraits, file = inputfile)

  }

  netfile = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("TCGA-networkConstruction-auto.",shortname,".RData"))
  if (file.exists(netfile)) {
    load(netfile)
  } else {
########################################################
## Choose a set of soft-thresholding powers
#  powers = c(c(1:10), seq(from = 12, to=20, by=2))
## Call the network topology analysis function
#  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
## Plot the results:
#  sizeGrWindow(9, 5)
#  par(mfrow = c(1,2));
#  cex1 = 0.9;
## Scale-free topology fit index as a function of the soft-thresholding power
#  pdf(file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", "scalefree.",shortname,".pdf"), width = 12, height = 9);
#  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#       main = paste("Scale independence"));
#  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#       labels=powers,cex=cex1,col="red");
## this line corresponds to using an R^2 cut-off of h
#  abline(h=0.90,col="red")
## Mean connectivity as a function of the soft-thresholding power
#  plot(sft$fitIndices[,1], sft$fitIndices[,5],
#       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#       main = paste("Mean connectivity"))
#  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#  dev.off()
#
##=====================================================================================

  TOMfile = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("TCGA.genes.",shortname,".TOM-block.1.RData"))
  if (file.exists(TOMfile)) {

      net = blockwiseModules(datExpr, power = powerparam, maxBlockSize = 30000,
                             networkType = "signed", TOMType = "signed", minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.30,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             loadTOM = TRUE,
                             saveTOMFileBase = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("TCGA.genes.",shortname,".TOM")), 
                             verbose = 3)
 
#      load(TOMfile)
#      dissTOM = 1-TOM
#      rm(TOM)
#      collectGarbage();
#
#      # Call the hierarchical clustering function
#      geneTree = hclust(as.dist(dissTOM), method = "average");
#      genetreefile = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("TCGA.genes.",shortname,".genetree.RData"))
#      save(geneTree, file=genetreefile)
#      # We like large modules, so we set the minimum module size relatively high:
#      minModuleSize = 30;
#      # Module identification using dynamic tree cut:
#      dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
#          deepSplit = 2, pamRespectsDendro = FALSE,
#          minClusterSize = minModuleSize);
#      table(dynamicMods)
#      # Convert numeric lables into colors
#      dynamicColors = labels2colors(dynamicMods)
#      table(dynamicColors)
#
#      # Calculate eigengenes
#      MEList = moduleEigengenes(datExpr, colors = dynamicColors)
#      MEs = MEList$eigengenes
#      # Calculate dissimilarity of module eigengenes
#      MEDiss = 1-cor(MEs);
#      # Cluster module eigengenes
#      METree = hclust(as.dist(MEDiss), method = "average");
#      # Plot the result
#      pdf(file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", "module.clusters.",shortname,".pdf"), width = 12, height = 9);
#      plot(METree, main = "Clustering of module eigengenes",
#            xlab = "", sub = "")
#      dev.off()
#
#      # Merge modules
#      MEDissThres = 0.25
#      # Plot the cut line into the dendrogram
#      abline(h=MEDissThres, col = "red")
#      # Call an automatic merging function
#      merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
#      # The merged module colors
#      mergedColors = merge$colors;
#      # Eigengenes of the new merged modules:
#      mergedMEs = merge$newMEs;
#      pdf(file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", "module.colors",shortname,".pdf"), width = 12, height = 9);
#      plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#            c("Dynamic Tree Cut", "Merged dynamic"),
#            dendroLabels = FALSE, hang = 0.03,
#            addGuide = TRUE, guideHang = 0.05)
#      dev.off()
#
#      # Rename to moduleColors
#      moduleColors = mergedColors
#      # Construct numerical labels corresponding to the colors
#      colorOrder = c("grey", standardColors(50));
#      moduleLabels = match(moduleColors, colorOrder)-1;
#      MEs = mergedMEs;
#
#      # Save module colors and labels for use in subsequent parts
#      save(MEs, moduleLabels, moduleColors, geneTree, file = netfile) 
#
    } else {

      net = blockwiseModules(datExpr, power = powerparam, maxBlockSize = 30000,
                             networkType = "signed", TOMType = "signed", minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.30,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE,
                             saveTOMFileBase = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("TCGA.genes.",shortname,".TOM")), 
                             verbose = 3)
    }
    #=====================================================================================
    # open a graphics window
      #sizeGrWindow(12, 9)
      pdf(file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", "module.colors",shortname,".pdf"), width = 24, height = 9);
    # Convert labels to colors for plotting
      mergedColors = labels2colors(net$colors)
    # Plot the dendrogram and the module colors underneath
      plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)
      dev.off()
      moduleLabels = net$colors
      moduleColors = labels2colors(net$colors)
      MEs = net$MEs;
      geneTree = net$dendrograms[[1]];

      collectGarbage();
      save(net, MEs, moduleLabels, moduleColors, geneTree, 
           file = netfile )

  }

# Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  rownames(MEs) <- rownames(datExpr)
# Recalculate MEs with color labels
#  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#  rownames(MEs0) <- rownames(datExpr)
#  MEs = orderMEs(MEs0)
#  stopifnot(rownames(MEs) == rownames(MEs0))
# save MEs as representative expression value
  write.table(MEs, paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("moduleEigengenes.",shortname,".txt")), quote=FALSE, row.names=TRUE, sep="\t")



source("src/modulecor.nocontrol.r")
library("RDAVIDWebService")
library("ReactomePA")
library(org.Hs.eg.db)
library(annotate)


  ginfofile = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("geneInfo.",shortname,".RData"))
  if (file.exists(ginfofile)) {
    load(ginfofile)
  } else {

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


#file = gzfile(description = "GeneAnnotation.csv.gz");
  annot = read.table(paste0(datadir, "TCGA.annot.uniq.txt"), header=TRUE, sep="\t");
# Match genes in the data set to the probe IDs in the annotation file 
  genenames = colnames(datExpr)
  geneidx = !grepl (":", genenames)
  annotnames = paste(annot$GENE, annot$ENTREZID, sep="|")
  gene2annot = match(genenames, annotnames)

# The following is the number or probes without annotation:
  sum(is.na(gene2annot[geneidx]))
# Should return 0.




# Create the starting data frame
  geneInfo = data.frame(Gene = genenames, GeneSymbol = annot$GENE[gene2annot],
               EntrezID = annot$ENTREZID[gene2annot],
                        moduleLabel = moduleLabels,
                        moduleColor = moduleColors
                        )
# Order modules by their significance for current_trait
#  modOrder = order(-abs(cor(MEs, current_trait, use = "p")));
# Add module membership information in the chosen order

  modNames = substring(names(MEs), 3)

  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");

  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo)
    geneInfo = data.frame(geneInfo, geneModuleMembership[, mod], 
                           MMPvalue[, mod]);
    names(geneInfo) = c(oldNames, paste("MM.", modNames[mod], sep=""),
                         paste("p.MM.", modNames[mod], sep=""))
  }
# Order the genes in the geneInfo variable first by module color
  geneOrder = order(geneInfo$moduleColor);
  geneInfo = geneInfo[geneOrder, ]


  save(geneInfo, file = ginfofile )

}

#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


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



  # Convert numerical lables to colors for labeling of modules in the plot
  MEnumeric = names(MEs)
  MEColors = labels2colors(as.numeric(substring(MEnumeric, 3)));
  MEColorNames = paste("ME", MEColors, sep="");
modulenumbers <- as.numeric(substring(MEnumeric, 3))
for( i in modulenumbers) {
  if (i == 0) {print ("skipping grey"); next;}
  reportFileName=paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste("termClusterReport",i,labels2colors(i),"txt",sep="."))
  if (file.exists(reportFileName)) {print (paste0("skipping ", reportFileName)); next}
  david<-DAVIDWebService(email="mira.han@unlv.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setTimeOut(david, 80000)
  setAnnotationCategories(david, getDefaultCategoryNames(david))
  #genelist <- geneInfo[geneInfo$moduleColor==color,"EntrezID"]
  genelist <- geneInfo[geneInfo$moduleLabel==i,"EntrezID"]
  result<-addList(david, genelist, idType="ENTREZ_GENE_ID", listName=paste("module", i, labels2colors(i)), listType="Gene")
  termCluster<-getClusterReport(david, type="Term")
  if (length(members(termCluster)) > 0) {
    symbols <- getClusterSymbols(termCluster)
    symbols <- lapply(symbols, sort)
    lapply(symbols, write, paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste("termClusterSymbols",i,labels2colors(i),"txt", sep=".")), append=TRUE, ncolumns=1000)
    Sys.sleep(1)
    getClusterReportFile(david, type="Term", fileName=reportFileName)
    Sys.sleep(1)
  }
}



for( i in modulenumbers) {
  if (i == 0) {print ("skipping grey"); next;}
  if (file.exists(paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", paste0("reactome.",shortname,".",i,".",labels2colors(i),".em.pdf")))) {print (paste0("skipping reactome ", labels2colors(i))); next}
  genelist <- geneInfo[geneInfo$moduleLabel==i,"EntrezID"]
  x <- enrichPathway(gene=genelist[!is.na(genelist)], pvalueCutoff=0.05, readable=T)
  Sys.sleep(1) 
  write.table(x, file = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste("reactome", i, labels2colors(i), "txt", sep=".")), quote=FALSE, row.names=TRUE, sep="\t")
  if (! is.null(x) && nrow(x)) {
    dotplot(x, showCategory=15)
    ggsave(filename=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", paste0("reactome.",shortname,".",i,".",labels2colors(i),".pdf")), width = 8, height = 12, units = "in");
    emapplot(x)
    ggsave(filename=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", paste0("reactome.",shortname,".",i,".",labels2colors(i),".em.pdf")), width = 8, height = 12, units = "in");
  }
}




  #moduleTraitCor = cor(MEs, datTraits, use = "p");
  #moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  module2L1HS = lm_module2L1HS(MEs, allTraits, paste0(resultdir, "WGCNA/bytissue/", shortname, "/"))
  module2oldLINE = lm_module2oldLINE(MEs, allTraits, paste0(resultdir, "WGCNA/bytissue/", shortname, "/"))
  moduleTraitCor = data.matrix(cbind(module2L1HS$moduleexp_coef, module2oldLINE$moduleexp_coef))
  rownames(moduleTraitCor) = colnames(MEs)
  colnames(moduleTraitCor) = colnames(datTraits)
  moduleTraitPvalue = data.matrix(cbind(module2L1HS$moduleexp_pval, module2oldLINE$moduleexp_pval))
  rownames(moduleTraitPvalue) = colnames(MEs)
  colnames(moduleTraitPvalue) = colnames(datTraits)


  #sizeGrWindow(10,6)

  pdf(file=paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", paste0("module.trait.cor.",shortname,".pdf")), width = 8, height = 12);
# Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                             signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
  zlimit = ceiling(max(abs(moduleTraitCor)))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = MEnumeric,
                 ySymbols = MEnumeric,
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-zlimit, zlimit),
                 main = paste("Module-trait relationships in", shortname))
  dev.off()



for (trait in 1:ncol(datTraits))
{
# Define variable current_trait containing the current_trait column of datTrait
  current_trait = as.data.frame(datTraits[,trait]);
  names(current_trait) = colnames(datTraits)[trait] 

  geneInfofile = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("geneInfo.", names(current_trait),".", shortname,".txt"))
  if (file.exists(geneInfofile)) {print (paste0("skipping ", names(current_trait))); next}

# names (colors) of the modules
  modNames = substring(names(MEs), 3)

  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");

  geneTraitSignificance = as.data.frame(cor(datExpr, current_trait, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

  names(geneTraitSignificance) = "GS";
  names(GSPvalue) = "p.GS"

  for (module in modNames) {
    column = match(module, modNames);
    moduleGenes = (moduleLabels==module);
    pdf(paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", "module.gene.sig.",shortname,".",names(current_trait),".",module,".pdf"))
    sizeGrWindow(7, 7);
    par(mfrow = c(1,1));
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", names(current_trait)),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()
  }


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


#names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


#names(datExpr)[moduleColors=="yellow"]




#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
  geneInfo0 = data.frame(Gene = genenames, GeneSymbol = annot$GENE[gene2annot],
               EntrezID = annot$ENTREZID[gene2annot],
                        moduleLabel = moduleLabels,
                        moduleColor = moduleColors,
                        geneTraitSignificance,
                        GSPvalue)
# Order modules by their significance for current_trait
  modOrder = order(-abs(cor(MEs, current_trait, use = "p")));
# Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS));
  geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


  write.table(geneInfo, file=geneInfofile, quote=FALSE, row.names=TRUE, sep="\t")
}




#=============================================
# plot heatmap
#==============================================

  TOMfile = paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste0("TCGA.genes.",shortname,".TOM-block.1.RData"))
  if (file.exists(TOMfile)) {
    load(TOMfile)
    tom = as.matrix(TOM)
    rm(TOM)
    collectGarbage()

    color1=moduleColors
    restGenes= (color1 != "grey")
    dissTom = 1 - tom[restGenes,restGenes]
    rm(tom)
    collectGarbage()

    # randomly plot 5000 genes
#    nSelect = 5000
#    # For reproducibility, we set the random seed
#    set.seed(10);
#    select = sample(restGenes, size = nSelect);
#    selectTom = dissTom[select, select];  
#    hier1=hclust(as.dist(selectTom), method="average" )

    hier1=hclust(as.dist(dissTom), method="average" )
    diag(dissTom) = NA; 
    #sizeGrWindow(7,7)
    pdf(paste0(resultdir, "WGCNA/bytissue/", shortname, "/Plots/", "TOMheatmap.pdf"))
    TOMplot(dissTom^4, hier1, as.character(color1[restGenes]),
         main = "TOM heatmap plot, module genes" )
    dev.off()


  }







