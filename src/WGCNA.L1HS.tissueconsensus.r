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

colnames(traitData)[colnames(traitData) == "VSTcnts"] = "L1HS"
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


# Choose a set of soft-thresholding powers
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

# linear model with covariate
  rownames(consMEs[[set]]$data) = rownames(multiExpr[[set]]$data);
  module2L1HS = lm_module2L1HS(consMEs[[set]]$data, Traits[[set]]$data, paste0(resultdir, "WGCNA/consensus/"))
  module2oldLINE = lm_module2oldLINE(consMEs[[set]]$data, Traits[[set]]$data,  paste0(resultdir, "WGCNA/consensus/"))

  setModuleTraitCor = data.matrix(cbind(module2L1HS$moduleexp_coef, module2oldLINE$moduleexp_coef))
  rownames(setModuleTraitCor) = colnames(consMEs[[set]]$data)
  colnames(setModuleTraitCor) = c("L1HS", "oldLINE")
  
  setModuleTraitPvalue = data.matrix(cbind(module2L1HS$moduleexp_pval, module2oldLINE$moduleexp_pval))
  rownames(setModuleTraitPvalue) = colnames(consMEs[[set]]$data)
  colnames(setModuleTraitPvalue) = c("L1HS", "oldLINE")

  setModuleTraitEta2 = data.matrix(cbind(module2L1HS$partialeta2, module2oldLINE$partialeta2))
  rownames(setModuleTraitEta2) = colnames(consMEs[[set]]$data)
  colnames(setModuleTraitEta2) = c("L1HS", "oldLINE")

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
               xLabels = c("L1HS", "oldLINE"),
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
               xLabels = c("L1HS", "oldLINE"),
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
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data[,c("L1HS", "oldLINE")]);
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
traitNames = c("L1HS", "oldLINE")
#nTraits = checkSets(Traits)$nGenes
nTraits = length(traitNames)
dim(GSmat) = c(nGenes, (2*nSets+2)*nTraits)
rownames(GSmat) = genenames;
colnames(GSmat) = spaste(
    c(paste0("GS.set", 1:11), paste0("p.GS.set", 1:11), "Z.GS.meta.", "p.GS.meta"), 
    rep(traitNames, rep((2*nSets+2), nTraits)))
# Same code for kME:
#kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
kMEmat = do.call(rbind, lapply(kME, "[[", "cor")) 
kMEmat = rbind(kMEmat, do.call(rbind, lapply(kME, "[[", "p")))
kMEmat = rbind(kMEmat, kME.metaZ, kME.metaP); 
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, (2*nSets+2)*nMEs)
rownames(kMEmat) = genenames;
colnames(kMEmat) = spaste(
    c(paste0("kME.set", 1:11), paste0("p.kME.set", 1:11), "Z.kME.meta.", "p.kME.meta"), 
    rep(MEnames, rep((2*nSets+2), nMEs)))


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

info = data.frame(Gene = genenames, GeneSymbol = annot$GENE[gene2annot],
             EntrezID = annot$ENTREZID[gene2annot],
             ModuleLabel = moduleLabels,
             ModuleColor = labels2colors(moduleLabels),
             GSmat,
             kMEmat);
write.csv(info, file = paste0(resultdir, "WGCNA/consensus/consensusAnalysis-CombinedNetworkResults.csv"), row.names = FALSE, quote = FALSE);

interestingmodules = which((is.na(consensusCor[,1])&is.na(consensusCor[,2]))==FALSE)
#for( i in modulenumbers[interestingmodules]) {
for( i in modulenumbers) {
  if (i == 0) {print ("skipping grey"); next;}
  reportfileName=paste0(resultdir, "WGCNA/consensus/", paste("termClusterReport",i,labels2colors(i),"txt",sep="."))
  if (file.exists(reportfileName)) {print (paste0("skipping ", reportfileName)); next}

  david<-DAVIDWebService(email="mira.han@unlv.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setAnnotationCategories(david, getDefaultCategoryNames(david))
  genelist <- info[info$ModuleLabel==i,"EntrezID"]
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
  #reportFileName=paste0(resultdir, "WGCNA/bytissue/", shortname, "/", paste("termClusterReport",i,labels2colors(i),"txt",sep="."))
  #if (file.exists(reportFileName)) {print (paste0("skipping ", reportFileName)); next}
  genelist <- info[info$ModuleLabel==i,"EntrezID"]
  genelist = genelist[!is.na(genelist)]
  x <- enrichPathway(gene=genelist, pvalueCutoff=0.05, readable=T)
  if (! is.null(x) && nrow(x)) {
    write.table(as.data.frame(x), file = paste0(resultdir, "WGCNA/consensus/", paste("Reactome",i,labels2colors(i),"txt", sep=".") ), row.names = FALSE, quote = FALSE);
    dotplot(x, showCategory=15)
    ggsave(filename=paste0(resultdir, "WGCNA/consensus/Plots/", paste("Reactome",i,labels2colors(i),"pdf", sep=".")), width=16, height=(min(nrow(as.data.frame(x)),15)*0.32), units = "in", limitsize = FALSE )
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

power=6
color1=moduleColors
  restGenes= (color1 != "grey")
diss1=1-TOMsimilarityFromExpr( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
       main = "TOM heatmap plot, module genes" )


power=6
color1=moduleColors
restGenes= (color1 != "grey")
diss1=1-adjacency( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
       main = "Adjacency heatmap plot, module genes" )


sizeGrWindow(7,7)
topList=rank(NS1$p.Weighted,ties.method="first")<=30
gene.names= names(datExpr)[topList]
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                 networkType="signed", useTOM=FALSE,
                 power=1, main="signed correlations")

sizeGrWindow(7,7)
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                 networkType="unsigned", useTOM=FALSE,
                 power=1, main="signed correlations")

sizeGrWindow(7,7)
# The following shows the TOM heatmap in a signed network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                 networkType="signed", useTOM=TRUE,
                 power=12, main="C. TOM in a signed network")
# The following shows the TOM heatmap in a unsigned network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                 networkType="unsigned", useTOM=TRUE,
                 power=6, main="D. TOM in an unsigned network")


