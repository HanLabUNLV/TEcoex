#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "..";
setwd(workingDir); 
# Load the WGCNA package
#library(DESeq2)
library(WGCNA);
library("ComplexHeatmap")
library("circlize")
library("RDAVIDWebService")
library("ReactomePA")
library(org.Hs.eg.db)
library(annotate)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads()
source("src/modulecor.nocontrol.r")
datadir = "data/"
resultdir = "result.rsem.TET/"
args = commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
}

dir.create(paste0(resultdir, "WGCNA")) 
dir.create(paste0(resultdir, "WGCNA/consensus/")) 
dir.create(paste0(resultdir, "WGCNA/consensus/Plots")) 


if (file.exists(paste0(resultdir, "WGCNA/consensus/Consensus-dataInput.RData"))) {
  load(paste0(resultdir, "WGCNA/consensus/Consensus-dataInput.RData")) 
  nSets = length(Traits);

} else {

#Read in the female liver data set
#normalcnts <- read.table(paste0(resultdir, "/normalizedcnts.gene.txt"), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
normalcnts <- read.table(paste0(resultdir, "VSTcnt.txt"), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
# Take a quick look at what is in the data set:
colnames(normalcnts) = substr(colnames(normalcnts), 1,12)
dim(normalcnts);
names(normalcnts);


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

traitData = read.table(paste0(resultdir, "L1HS5prime.VST.txt"), header=TRUE, sep="\t");
log2TPMsum <- read.table(file=paste0(resultdir, "log2colsums.txt"), sep="\t", header = TRUE, row.names=1)
log2TPMsum <- cbind(substr(rownames(log2TPMsum), 1, 12), log2TPMsum)
colnames(log2TPMsum) <- c("patient", "log2TPMsum")
traitData <- merge(traitData, log2TPMsum, by = "patient")

oldLINEVST <- read.table(file=paste0(resultdir, "oldLINE.VST.txt"), sep="\t", header = TRUE, row.names=1)
colnames(oldLINEVST) <- c("patient", "oldLINEVST")
traitData <- merge(traitData, oldLINEVST, by = "patient")

THCAradiation = read.table(file="data/THCA.radiation.txt", sep="\t", header=FALSE, row.names=1)
THCAradiation <- cbind.data.frame(substr(rownames(THCAradiation), 1, 12), THCAradiation)
colnames(THCAradiation) <- c("patient", "THCAradiation")
traitData$radiation = rep("NA", nrow(traitData))
rownames(traitData) = traitData$patient
commonpatients <- traitData[as.character(THCAradiation$patient),"patient"]
commonpatients <- commonpatients[!is.na(commonpatients)]
traitData[as.character(commonpatients),"radiation"] = THCAradiation[as.character(commonpatients),"THCAradiation"] 

traitData$newtype = traitData$tissue
traitData[traitData$tissue == "COAD", "newtype"] = "COADREAD"
traitData[traitData$tissue == "READ", "newtype"] = "COADREAD"
traitData[traitData$tissue == "ESCA", "newtype"] = "ESCASTAD"
traitData[traitData$tissue == "STAD", "newtype"] = "ESCASTAD"
traitData[traitData$tissue == "KICH", "newtype"] = "KICHKIRCKIRP"
traitData[traitData$tissue == "KIRC", "newtype"] = "KICHKIRCKIRP"
traitData[traitData$tissue == "KIRP", "newtype"] = "KICHKIRCKIRP"
traitData[traitData$tissue == "LUSC", "newtype"] = "LUADLUSC"
traitData[traitData$tissue == "LUAD", "newtype"] = "LUADLUSC"

colnames(traitData)[colnames(traitData) == "VSTcnts"] = "L1HS.5prime"
colnames(traitData)[colnames(traitData) == "VSTcnts.1"] = "L1HS"
colnames(traitData)[colnames(traitData) == "oldLINEVST"] = "oldLINE"
dim(traitData)
names(traitData)
print(summary(traitData))
write.table(traitData, paste0(resultdir, "WGCNA/", "WGCNA.L1HS.clinical.txt"), quote=FALSE, sep="\t", row.names=FALSE);


# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = levels(factor(traitData$newtype))
nSets = length(setLabels);
shortLabels = c("bladder", "breast", "colorectum", "esophagusstomach", "headneck", "kidney", "liver", "lung", "prostate", "thyroid", "endometrium")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

for (i in 1:nSets) {
samples = traitData[traitData$newtype==setLabels[i],]$patient
multiExpr[[i]] = list(data = as.data.frame(t(normalcnts[,samples])));
names(multiExpr[[i]]$data) = rownames(normalcnts);
rownames(multiExpr[[i]]$data) = samples;
}
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/SampleClustering.pdf"), width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();


# Choose the "base" cut height for the female data set
#baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
#cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# Re-plot the dendrograms including the cut lines
#pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
#par(mfrow=c(2,1))
#par(mar = c(0, 4, 2, 0))
#for (set in 1:nSets)
#{
#  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
#       xlab="", sub="", cex = 0.7);
#  abline(h=cutHeights[set], col = "red");
#}
#dev.off();

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================



#for (set in 1:nSets)
#{
#  # Find clusters cut by the line
#  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
#  # Keep the largest one (labeled by the number 1)
#  keep = (labels==1)
#  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
#}
#collectGarbage();
## Check the size of the leftover data
#exprSize = checkSets(multiExpr)
#exprSize


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


#traitData = read.table("WGCNA.L1HS.clinical.txt", header=TRUE, sep="\t");
dim(traitData)
names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, traitData$patient);
  allTraits = traitData[traitRows,]
  rownames(allTraits) = traitData[traitRows, 1];
  Traits[[set]] = list(data = allTraits);
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;



#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = paste0(resultdir, "WGCNA/consensus/Consensus-dataInput.RData"));

}



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

if (file.exists(paste0(resultdir, "WGCNA/consensus/Consensus-NetworkConstruction-auto.RData"))) {
  load(paste0(resultdir, "WGCNA/consensus/Consensus-NetworkConstruction-auto.RData"))
} else {


## Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
#sizeGrWindow(8, 6)
pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/scaleFreeAnalysis.pdf"), wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================



net = blockwiseConsensusModules(
        multiExpr, maxBlockSize = 25000,  power = 14, minModuleSize = 30, deepSplit = 2,
        pamRespectsDendro = FALSE, networkType = "signed", TOMType = "signed",
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = TRUE, 
        saveTOMFileBase = paste0(resultdir, "WGCNA/consensus/TCGA.genes.TOM"), 
        verbose = 5)



#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]; 

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


save(net, consMEs, moduleLabels, moduleColors, consTree, file = paste0(resultdir, "WGCNA/consensus/Consensus-NetworkConstruction-auto.RData"))


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


#sizeGrWindow(8,6);
pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/ConsensusDendrogram-auto.pdf"), wi = 16, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

dev.off()


}



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
moduleTraitEta2 = list();
# Calculate the correlations
nSets = length(Traits)
for (set in 1:nSets)
{
#  classic correlation
#  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data[,c("L1HS", "oldLINE")], use = "p");
#  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);

  traitmat <- Traits[[set]]$data
  traitmat$age = as.numeric(as.character(traitmat$age_at_diagnosis))/365
  Traits[[set]]$data = traitmat
# linear model with covariate
  rownames(consMEs[[set]]$data) = rownames(multiExpr[[set]]$data);
  module2L1HS = lm_module2trait(consMEs[[set]]$data, traitmat, paste0(resultdir, "WGCNA/consensus/"), "L1HS.5prime")
  module2oldLINE = lm_module2trait(consMEs[[set]]$data, traitmat,  paste0(resultdir, "WGCNA/consensus/"), "oldLINE")
  module2gender = lm_trait2module(consMEs[[set]]$data, traitmat,  paste0(resultdir, "WGCNA/consensus/"), "gender", TRUE)
  module2age = lm_trait2module(consMEs[[set]]$data, traitmat,  paste0(resultdir, "WGCNA/consensus/"), "age", FALSE)

  setModuleTraitCor = data.matrix(cbind(module2L1HS$moduleexp_coef, module2oldLINE$moduleexp_coef, module2gender$moduleexp_coef, module2age$moduleexp_coef))
  rownames(setModuleTraitCor) = colnames(consMEs[[set]]$data)
  colnames(setModuleTraitCor) = c("L1HS.5prime", "oldLINE", "gender", "age")
  
  setModuleTraitPvalue = data.matrix(cbind(module2L1HS$moduleexp_pval, module2oldLINE$moduleexp_pval, module2gender$moduleexp_pval, module2age$moduleexp_pval))
  rownames(setModuleTraitPvalue) = colnames(consMEs[[set]]$data)
  colnames(setModuleTraitPvalue) = c("L1HS.5prime", "oldLINE", "gender", "age")

  setModuleTraitEta2 = data.matrix(cbind(module2L1HS$partialeta2, module2oldLINE$partialeta2, module2gender$partialeta2, module2age$partialeta2))
  rownames(setModuleTraitEta2) = colnames(consMEs[[set]]$data)
  colnames(setModuleTraitEta2) = c("L1HS.5prime", "oldLINE", "gender", "age")

  moduleTraitCor[[set]] = setModuleTraitCor;
  moduleTraitPvalue[[set]] = setModuleTraitPvalue;
  moduleTraitEta2[[set]] = setModuleTraitEta2;
#
}


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# Convert numerical lables to colors for labeling of modules in the plot
MEnumeric = names(consMEs[[1]]$data) # modules are shared acorss sets
MEColors = labels2colors(as.numeric(substring(MEnumeric, 3)));
MEColorNames = paste("ME", MEColors, sep="");
modulenumbers <- as.numeric(substring(MEnumeric, 3))


for (set in 1:nSets) {

# Open a suitably sized window (the user should change the window size if necessary)
#sizeGrWindow(10,7)
pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/ModuleTraitRelationships-",shortLabels[set],".pdf"), wi = 10, he = 12);
# Plot the module-trait relationship table for set number 1
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                           signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
zlimit = ceiling(max(abs(moduleTraitEta2[[set]])))
colorlevel = sign(moduleTraitCor[[set]])*moduleTraitEta2[[set]]
#print(moduleTraitCor[[set]])
#print(sign(moduleTraitCor[[set]]))
#print(moduleTraitEta2[[set]])
#print(colorlevel)
labeledHeatmap(Matrix = colorlevel,
               xLabels = c("L1HS.5prime", "oldLINE", "gender", "age"),
               yLabels = MEnumeric,
               ySymbols = MEnumeric,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module--trait relationships in", shortLabels[set]))
dev.off();
}

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusEta2 = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# Find consensus negative correlations
negative = matrix(TRUE, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
for (set in 1:nSets) {
  negative = negative & ( moduleTraitCor[[set]] <= 0 );
}
args_cor <- vector("list", nSets)
args_pval <- vector("list", nSets)
args_eta2 <- vector("list", nSets)
for (set in 1:nSets) {
  args_cor[[set]] = moduleTraitCor[[set]][negative ]
  args_pval[[set]] = moduleTraitPvalue[[set]][negative ]
  args_eta2[[set]] = moduleTraitEta2[[set]][negative ]
}
#consensusCor[negative] = do.call(pmax, args_cor)
#consensusPvalue[negative] = do.call(pmax, args_pval)
args_cor_mat <- t(matrix(t(unlist(args_cor)), ncol=nSets))
args_cor_mat[args_cor_mat==0] = NA
#consensusCor[negative] = apply(args_cor_mat, 2, max, na.rm=TRUE)
consensusCor[negative] = apply(args_cor_mat, 2, median, na.rm=TRUE)

args_pval_mat <- t(matrix(t(unlist(args_pval)), ncol=nSets))
args_pval_mat[is.na(args_cor_mat)] = NA
#consensusPvalue[negative] = apply(args_pval_mat, 2, max, na.rm=TRUE)
consensusPvalue[negative] = apply(args_pval_mat, 2, median, na.rm=TRUE)

args_eta2_mat <- t(matrix(t(unlist(args_eta2)), ncol=nSets))
args_eta2_mat[is.na(args_cor_mat)] = NA
consensusEta2[negative] = apply(args_eta2_mat, 2, median, na.rm=TRUE)

# Find consensus positive correlations
positive = matrix(TRUE, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
for (set in 1:nSets) {
  positive = positive & ( moduleTraitCor[[set]] >= 0 );
}
args_cor <- vector("list", nSets)
args_pval <- vector("list", nSets)
args_eta2 <- vector("list", nSets)
for (set in 1:nSets) {
  args_cor[[set]] = moduleTraitCor[[set]][positive ]
  args_pval[[set]] = moduleTraitPvalue[[set]][positive ]
  args_eta2[[set]] = moduleTraitEta2[[set]][positive ]
}
#consensusCor[positive] = do.call(pmin, args_cor)
#consensusPvalue[positive] = do.call(pmax, args_pval)
args_cor_mat <- t(matrix(t(unlist(args_cor)), ncol=nSets))
args_cor_mat[args_cor_mat==0] = NA
#consensusCor[positive] = apply(args_cor_mat, 2, min, na.rm=TRUE)
consensusCor[positive] = apply(args_cor_mat, 2, median, na.rm=TRUE)

args_pval_mat <- t(matrix(t(unlist(args_pval)), ncol=nSets))
args_pval_mat[is.na(args_cor_mat)] = NA
#consensusPvalue[positive] = apply(args_pval_mat, 2, max, na.rm=TRUE)
consensusPvalue[positive] = apply(args_pval_mat, 2, median, na.rm=TRUE)

args_eta2_mat <- t(matrix(t(unlist(args_eta2)), ncol=nSets))
args_eta2_mat[is.na(args_cor_mat)] = NA
consensusEta2[positive] = apply(args_eta2_mat, 2, median, na.rm=TRUE)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


textMatrix =  paste(signif(consensusCor, 2), "\n(",
                           signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(consensusCor)
#sizeGrWindow(10,7)
pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/ModuleTraitRelationships-consensus.pdf"), wi = 10, he = 12);
par(mar = c(6, 8.8, 3, 2.2));
zlimit = ceiling(max(abs(consensusEta2), na.rm=TRUE))
colorlevel = sign(consensusCor)*consensusEta2
labeledHeatmap(Matrix = colorlevel,
               xLabels = c("L1HS.5prime", "oldLINE", "gender", "age"),
               yLabels = MEnumeric,
               ySymbols = MEnumeric,
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))
dev.off()

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


#file = gzfile(description = "GeneAnnotation.csv.gz");
annot = read.table(paste0(datadir, "TCGA.annot.uniq.txt"), header=TRUE, sep="\t");
# Match genes in the data set to the probe IDs in the annotation file 
genenames = names(multiExpr[[1]]$data)
geneidx = !grepl (":", genenames)
annotnames = paste(annot$GENE, annot$ENTREZID, sep="|")
gene2annot = match(genenames, annotnames)
# The following is the number or probes without annotation:
sum(is.na(gene2annot[geneidx]))
# Should return 0.

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
MET = consensusOrderMEs(consMEs.unord);
pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/EigengeneNetworks.pdf"), width= 12, height = 12);
par(cex = 0.9)

  for (set in 1:nSets) {
    greyLabel = 0
      labels = names(MET[[set]]$data)
      uselabels = labels[substring(labels, 3) != greyLabel]
      corME = cor(MET[[set]]$data[substring(labels, 
            3) != greyLabel, substring(labels, 3) != greyLabel], 
          use = "p")
      disME = as.dist(1 - corME)
      clust = fastcluster::hclust(disME, method = "average")
      main = setLabels[set]
      plotLabels = uselabels
      plot(clust, main = main, sub = "", xlab = "", labels = plotLabels, 
          ylab = "", ylim = c(0, 1))
  }
dev.off()

# plot heatmap
pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/EigengeneHeatmaps.pdf"), width= 12, height = 12);
for (set in 1:nSets) {
  corME = cor(MET[[set]]$data, use = "p")
  pME = corPvalueFisher(corME, nrow(MET[[set]]$data))
      labeledHeatmap(corME, names(MET[[set]]$data), 
      names(MET[[set]]$data), main = setLabels[[set]], invertColors = FALSE, 
      zlim = c(-1, 1), colorLabels = FALSE 
      )
}

dev.off();


# write moduleEigengenes
for (set in 1:nSets) {
  eigenGenes = MET[[set]]$data
  rownames(eigenGenes) = rownames(multiExpr[[set]]$data);
  write.table(eigenGenes, paste0(resultdir, "WGCNA/consensus/", paste0("moduleEigengenes.",setLabels[[set]],".txt")), quote=FALSE, row.names=TRUE, sep="\t")
}

GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data[,c("L1HS.5prime", "oldLINE")]);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


GS.metaZ = Reduce('+', lapply(GS, "[[", "Z"))/sqrt(nSets); #(GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
kME.metaZ = Reduce('+', lapply(kME, "[[", "Z"))/sqrt(nSets); #(kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


#GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
GSmat = do.call(rbind, lapply(GS, "[[", "cor")) 
GSmat = rbind(GSmat, do.call(rbind, lapply(GS, "[[", "p")))
GSmat = rbind(GSmat, GS.metaZ, GS.metaP); 
GSsmall = rbind(GS.metaZ, GS.metaP)
traitNames = c("L1HS.5prime", "oldLINE")
#nTraits = checkSets(Traits)$nGenes
nTraits = length(traitNames)
dim(GSmat) = c(nGenes, (2*nSets+2)*nTraits)
dim(GSsmall) = c(nGenes, 2*nTraits)
rownames(GSmat) = genenames;
rownames(GSsmall) = genenames;
colnames(GSmat) = spaste(
    c(paste0("GS.set", 1:11), paste0("p.GS.set", 1:11), "Z.GS.meta.", "p.GS.meta"), 
    rep(traitNames, rep((2*nSets+2), nTraits)))
colnames(GSsmall) = spaste(c("Z.GS.meta.", "p.GS.meta"),rep(traitNames, rep(2, nTraits)))



# Same code for kME:
#kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
kMEcor = do.call(rbind, lapply(kME, "[[", "cor")) 
kMEmat = rbind(kMEcor, do.call(rbind, lapply(kME, "[[", "p")))
kMEmat = rbind(kMEmat, kME.metaZ, kME.metaP); 
kMEsmall = rbind(kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEcor) = c(nGenes, nSets*nMEs)
dim(kMEmat) = c(nGenes, (2*nSets+2)*nMEs)
dim(kMEsmall) = c(nGenes, 2*nMEs)
rownames(kMEcor) = genenames;
rownames(kMEmat) = genenames;
rownames(kMEsmall) = genenames;
colnames(kMEcor) = outer(c(paste0("kME.set", 1:11)), MEnames, FUN = "paste0")[1:(nSets*nMEs)] 
colnames(kMEmat) = spaste(
    c(paste0("kME.set", 1:11), paste0("p.kME.set", 1:11), "Z.kME.meta.", "p.kME.meta"), 
    rep(MEnames, rep((2*nSets+2), nMEs)))
colnames(kMEsmall) = spaste(
    c("Z.kME.meta.", "p.kME.meta"), 
    rep(MEnames, rep(2, nMEs)))


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

info.GS = data.frame(Gene = genenames, GeneSymbol = annot$GENE[gene2annot],
             EntrezID = annot$ENTREZID[gene2annot],
             ModuleLabel = moduleLabels,
             ModuleColor = labels2colors(moduleLabels),
             GSsmall);
write.csv(info.GS, file = paste0(resultdir, "WGCNA/consensus/small.GS.csv"), row.names = FALSE, quote = FALSE);

info.kME = data.frame(Gene = genenames, GeneSymbol = annot$GENE[gene2annot],
             EntrezID = annot$ENTREZID[gene2annot],
             ModuleLabel = moduleLabels,
             ModuleColor = labels2colors(moduleLabels),
             kMEsmall);
info.kME.cor = data.frame(Gene = genenames, GeneSymbol = annot$GENE[gene2annot],
             EntrezID = annot$ENTREZID[gene2annot],
             ModuleLabel = moduleLabels,
             ModuleColor = labels2colors(moduleLabels),
             kMEcor);
write.csv(info.kME, file = paste0(resultdir, "WGCNA/consensus/small.kMEmat.csv"), row.names = FALSE, quote = FALSE);
write.csv(info.kME.cor, file = paste0(resultdir, "WGCNA/consensus/kME.membership.csv"), row.names = FALSE, quote = FALSE);


# plot genes correlated with TEs
plotGeneIdx = rep(FALSE, nGenes)
# find modules with more than 50% transposons 
# write a BED file of TSS's (+-1kb window) of genes in the TE modules
TSS = read.table(paste0(datadir, "TCGA.TSS.uniq.txt"), header=TRUE, sep="\t");
rownames(TSS) = paste(TSS[,1], TSS[,2], sep="|")

TEproportion = c()
TEmodules = c()
for( i in modulenumbers) {
  genesinmodule = (info.kME.cor$ModuleLabel == i)
  TEsinmodule = genesinmodule & grepl(":", info.kME.cor$Gene)
  TEproportion = c(TEproportion, sum(TEsinmodule)/sum(genesinmodule))
  if (sum(TEsinmodule) > sum(genesinmodule)*0.5) {
#  if (sum(TEsinmodule) > sum(genesinmodule)*0.05) {
    TEmodules = c(TEmodules, i)
  }
}
genesinTEmodules = info.kME.cor$ModuleLabel %in% TEmodules 
consistenttissues = 1:nSets
consistenttissues= consistenttissues[c(-1, -4, -5)] # let's remove bladder esophagus and head and neck
for( i in TEmodules) {
  importantgenes = rep(TRUE, nGenes)
  for (set in consistenttissues) {
    collabel = paste0("kME.set", set, "ME", i)
    FilterGenes= kMEcor[,collabel]>.8
    write.table(info.kME.cor[FilterGenes, 1:5], file = paste0(resultdir, "WGCNA/consensus/TEmodule.", i, ".", setLabels[set], ".highMMgenes.txt"), row.names = FALSE, quote = FALSE)

    FilterGenes= kMEcor[,collabel]>.6
    importantgenes = importantgenes & FilterGenes
  }
  write.table(info.kME.cor[importantgenes, 1:5], file = paste0(resultdir, "WGCNA/consensus/TEmodule.", i, ".consensus.highMMgenes.txt"), row.names = FALSE, quote = FALSE)
  plotGeneIdx = plotGeneIdx | importantgenes # importantgenes are cor > 0.6 across all tissues 


  importantTEs = importantgenes & grepl(":", info.kME.cor$Gene)
  importantgenes = importantgenes & !grepl(":", info.kME.cor$Gene)
  TSS.importantgenes = TSS[info.kME.cor[importantgenes,1],]
  TSS.importantgenes$start = TSS.importantgenes[,4]-1000
  TSS.importantgenes$end = TSS.importantgenes[,4]+1000
  TSS.importantgenes$score = rep(".", nrow(TSS.importantgenes))
  write.table(TSS.importantgenes[,c(3, 6, 7, 1, 8, 5)], file = paste0(resultdir, "WGCNA/consensus/TEmodule.", i, ".TSS.bed"), row.names = FALSE, quote = FALSE, sep="\t")

}
plotGeneIdx.pos = plotGeneIdx
write.table(info.kME.cor[plotGeneIdx.pos, 1:5], file = paste0(resultdir, "WGCNA/consensus/labeltree.pos.txt"), row.names = FALSE, quote = FALSE, sep="\t")

#for( i in c(TEmodules,36)) {  # now print bladder esophagus and head and neck including cluster 36 that is specific to these tissues
#  importantgenes = rep(TRUE, nGenes)
#  for (set in c(1,4,5)) {
#    collabel = paste0("kME.set", set, "ME", i)
#    FilterGenes= kMEcor[,collabel]>.8
#    write.table(info.kME.cor[FilterGenes, 1:5], file = paste0(resultdir, "WGCNA/consensus/TEmodule.", i, ".", setLabels[set], ".highMMgenes.txt"), row.names = FALSE, quote = FALSE)
#  }
#}


# find modules negatively correlated with the TEmodules
negTEmodules = c()
for (set in 1:nSets) {
  corME = cor(MET[[set]]$data, use = "p")
  for (i in TEmodules) {
    rowlabel = paste0("ME", i)
    mcor = corME[rowlabel,]
    negcoridx = mcor < -0.7
    negTEmodules = c(negTEmodules, modulenumbers[negcoridx])
  }
}
negrank = rank(table(negTEmodules))
negTEmodules = as.numeric(names(negrank[negrank == max(negrank)]))
for( i in negTEmodules) {
  importantgenes = rep(TRUE, nGenes)
  for (set in consistenttissues) {
    collabel = paste0("kME.set", set, "ME", i)
    FilterGenes= kMEcor[,collabel]>.8
    write.table(info.kME.cor[FilterGenes, 1:5], file = paste0(resultdir, "WGCNA/consensus/negTEmodule.", i, ".", setLabels[set], ".highMMgenes.txt"), row.names = FALSE, quote = FALSE)

    FilterGenes= kMEcor[,collabel]>.6
    importantgenes = importantgenes & FilterGenes
  }
  write.table(info.kME.cor[importantgenes, 1:5], file = paste0(resultdir, "WGCNA/consensus/negTEmodule.", i, ".consensus.highMMgenes.txt"), row.names = FALSE, quote = FALSE)
  plotGeneIdx = plotGeneIdx | importantgenes 

  importantTEs = importantgenes & grepl(":", info.kME.cor$Gene)
  importantgenes = importantgenes & !grepl(":", info.kME.cor$Gene)
  TSS.importantgenes = TSS[info.kME.cor[importantgenes,1],]
  TSS.importantgenes$start = TSS.importantgenes[,4]-1000
  TSS.importantgenes$end = TSS.importantgenes[,4]+1000
  TSS.importantgenes$score = rep(".", nrow(TSS.importantgenes))
  write.table(TSS.importantgenes[,c(3, 6, 7, 1, 8, 5)], file = paste0(resultdir, "WGCNA/consensus/negTEmodule.", i, ".TSS.bed"), row.names = FALSE, quote = FALSE, sep="\t")


}
plotGeneIdx.neg = plotGeneIdx
plotGeneIdx.neg[plotGeneIdx.pos] = FALSE
write.table(info.kME.cor[plotGeneIdx.neg, 1:5], file = paste0(resultdir, "WGCNA/consensus/labeltree.neg.txt"), row.names = FALSE, quote = FALSE, sep="\t")


pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/heatmap.pdf"), width = 12, height = 12);
power=6
color1=moduleColors
restGenes= (color1 != "grey")
restGenes = genesinTEmodules
for (set in 1:nSets)
{
  datExpr = multiExpr[[set]]$data
  # The following shows the correlations between the top genes
  plotGenes = names(datExpr)[plotGeneIdx]
  match1 = match(plotGenes, colnames(datExpr))
  match1 = match1[!is.na(match1)]
  nGenes = length(match1)

  datErest = datExpr[, match1]
  ADJ1 = adjacency(datErest, weights = NULL, power = 1, 
      type = "signed")
  corrmat = 2*ADJ1 - 1
  diss1 = 1 - ADJ1
  diag(diss1) = NA
  hier1 = hclust(as.dist(diss1), method = "ave")

  #mat = datErest
  #mat <- mat - rowMeans(mat)
  #diss1 = dist(t(mat))
  #hier1 = hclust(d = diss1, method = "ave")

  labeltree = names(data.frame(datErest)) # rownames and colnames
  labelrow = names(data.frame(datErest))
  labelrow[hier1$order[length(labeltree):1]] = labelrow[hier1$order]
  #write.table(labeltree, file = paste0(resultdir, "WGCNA/consensus/heatmap.labeltree.", setLabels[set], ".txt"))
  write.table(labelrow, file = paste0(resultdir, "WGCNA/consensus/heatmap.labelrow.", setLabels[set], ".txt"))

  TE = grepl(":", labeltree)
  TE[grepl(".DNA$", labeltree)] = 2
  TE[grepl(".HERV$", labeltree)] = 3
  TE[grepl(".LTR$", labeltree)] = 4
  TE[grepl(".SINE$", labeltree)] = 5
  TE[grepl(".LINE$", labeltree)] = 6
  TE = as.factor(TE)

  genes = rep(0, length(labeltree))
  labeltreegenenames = unlist(lapply(strsplit(labeltree, "[.]"), `[[`, 1))
#  ZFPs = read.table(paste0(datadir, "/ZNF/ZNFlist.txt"), header=FALSE)
#  ZFPs = ZFPs[[1]]
#  ZFPs = substr(ZFPs, 1, nchar(ZFPs)-4) 
#  ZFPs = sub("[|]", ".", ZFPs)
#  ZFPmatch = match(labeltree, ZFPs)
#  genes[!is.na(ZFPmatch)] = 1
  #ribosome = read.table(paste0(datadir, "/GO_ribosome.msigDB.txt"), header=TRUE, sep="\t")
  ribosome = read.table(paste0(datadir, "/GO_Struct_Ribosome.msigDB.txt"), header=TRUE, sep="\t")
  ribosome = ribosome[-1,]
  ribomatch = match(labeltreegenenames, ribosome)
  genes[!is.na(ribomatch)] = 3
  mitochondria = read.table(paste0(datadir, "/GO_mitochondrion.msigDB.txt"), header=TRUE, sep="\t")
  mitochondria = mitochondria[-1,]
  mitomatch = match(labeltreegenenames, mitochondria)
  genes[!is.na(mitomatch)] = 2
  immune = read.table(paste0(datadir, "/GO_IMMUNE_SYSTEM_PROCESS.msigDB.txt"), header=TRUE, sep="\t")
  immune = immune[-1,]
  immunematch = match(labeltreegenenames, immune)
  genes[!is.na(immunematch)] = 4
  genes = as.factor(genes)
  annotation_data <- cbind.data.frame(TE, genes)
  
  col_vector1 = c("lightgrey", "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#d2f53c", "#fabebe", "#e6    beff", "#800000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
  col_vector2 = c("lightgrey", "purple", "orange", "yellowgreen", "#e6beff", "#800000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
  mycolors=list(TE = col_vector1[1:length(levels(annotation_data$TE))],
                genes = col_vector2[1:length(levels(annotation_data$genes))])
  names(mycolors$TE) <- levels(annotation_data$TE) 
  names(mycolors$genes) <- levels(annotation_data$genes) 

  ha = HeatmapAnnotation(annotation_data, col=mycolors)

  dend = as.dendrogram(hier1)
  #ht = Heatmap(as.matrix(diss1), cluster_rows = dend, row_dend_reorder = FALSE,
  ht = Heatmap(as.matrix(corrmat), cluster_rows = dend, row_dend_reorder = FALSE,
               name = "ht", cluster_columns = dend, column_dend_reorder = FALSE,
               top_annotation = ha,
               top_annotation_height = unit(4, "cm"),
               show_row_names = FALSE, show_column_names = FALSE,
  #             col = colorRamp2(c(0, 0.5, 1), c("blue", "#ffffbf", "red")),
               col = colorRamp2(c(-1, 0, 1), c("blue", "#ffffbf", "red")),
               column_title = "")

  draw(ht)
}
dev.off()





interestingmodules = which((is.na(consensusCor[,1])&is.na(consensusCor[,2]))==FALSE)
#for( i in modulenumbers[interestingmodules]) {
for( i in modulenumbers) {
  if (i == 0) {print ("skipping grey"); next;}
  reportfileName=paste0(resultdir, "WGCNA/consensus/", paste("termClusterReport",i,labels2colors(i),"txt",sep="."))
  if (file.exists(reportfileName)) {print (paste0("skipping ", reportfileName)); next}

  david<-DAVIDWebService(email="mira.han@unlv.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setAnnotationCategories(david, getDefaultCategoryNames(david))
  genelist <- info.kME[info.kME$ModuleLabel==i,"EntrezID"]
  result<-addList(david, genelist, idType="ENTREZ_GENE_ID", listName=paste("module", i, labels2colors(i)), listType="Gene")
  termCluster<-getClusterReport(david, type="Term")
  if (length(members(termCluster)) > 0) {
    symbols <- getClusterSymbols(termCluster)
    symbols <- lapply(symbols, sort)
    lapply(symbols, write, paste0(resultdir, "WGCNA/consensus/", paste("termClusterSymbols",i,labels2colors(i),"txt", sep=".")), append=TRUE, ncolumns=1000)
    getClusterReportFile(david, type="Term", fileName=reportfileName)
  }
}


for( i in modulenumbers) {
  if (i == 0) {print ("skipping grey"); next;}
  reportFileName=paste0(resultdir, "WGCNA/consensus/", paste("Reactome",i,labels2colors(i),"txt",sep="."))
  if (file.exists(reportFileName)) {print (paste0("skipping ", reportFileName)); next}
  genelist <- info.kME[info.kME$ModuleLabel==i,"EntrezID"]
  genelist = genelist[!is.na(genelist)]
  x <- enrichPathway(gene=genelist, pvalueCutoff=0.05, readable=T)
  if (! is.null(x) && nrow(x)) {
    write.table(as.data.frame(x), file = paste0(resultdir, "WGCNA/consensus/", paste("Reactome",i,labels2colors(i),"txt", sep=".") ), row.names = FALSE, quote = FALSE);
    dotplot(x, showCategory=15)
    ggsave(filename=paste0(resultdir, "WGCNA/consensus/Plots/", paste("Reactome",i,labels2colors(i),"pdf", sep=".")), width=12, height=(min(nrow(as.data.frame(x)),15)*0.32), units = "in", limitsize = FALSE )
    emapplot(x)
    ggsave(filename=paste0(resultdir, "WGCNA/consensus/Plots/", paste("Reactome",i,labels2colors(i),"em.pdf", sep=".")), width=8, height=12, units = "in")
  }
}



############################################
# igraph
#graph<-wgcna2igraph(net = net.1, datExpr = datExpr.1,
#                    modules2plot = c("blue","green","turquoise","brown"),
#                    colors2plot = c("orange","darkred","cyan","cornflowerblue"),
#                    kME.threshold = 0.5, adjacency.threshold = 0.1,
#                    adj.power = pow, verbose = T,
#                    node.size = 0, frame.color = NA, node.color = NA,
#                    edge.alpha = .5, edge.width =1)
#plot(graph)


#=============================================
# plot heatmap
#==============================================

pdf(file = paste0(resultdir, "WGCNA/consensus/Plots/TOMplot.pdf"), width = 12, height = 12);
power=6
color1=moduleColors
  restGenes= (color1 != "grey")
  restGenes = genesinTEmodules
for (set in 1:nSets)
{
datExpr = multiExpr[[set]]$data
diss1=1-TOMsimilarityFromExpr( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
       main = setLabels[[set]] )
}
dev.off()


