library("LOLA")
library("GenomicRanges")

dbPath = "/data2/han_lab/database/LOLA/mobiledna/hg19/"
regionDB = loadRegionDB(dbPath)
names(regionDB)
regionSetBLCA = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/bladder/TEmodule.TSS.TEs.bed")
regionSetBRCA = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/breast/TEmodule.TSS.TEs.bed")
regionSetCOAD = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/colonrectum/TEmodule.TSS.TEs.bed")
regionSetUCEC = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/endometrium/TEmodule.TSS.TEs.bed")
regionSetHNSC = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/headneck/TEmodule.TSS.TEs.bed")
regionSetKICH = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/kidney/TEmodule.TSS.TEs.bed")
regionSetLIHC = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/liver/TEmodule.TSS.TEs.bed")
regionSetLUAD = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/lung/TEmodule.TSS.TEs.bed")
regionSetPRAD = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/prostate/TEmodule.TSS.TEs.bed")
regionSetTHCA = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/WGCNA/bytissue/thyroid/TEmodule.TSS.TEs.bed")
userSets = GRangesList(regionSetBLCA, regionSetBRCA, regionSetCOAD, regionSetUCEC, regionSetHNSC, regionSetKICH, regionSetLIHC, regionSetLUAD, regionSetPRAD, regionSetTHCA)
uniqTEs = readBed("/data2/han_lab/mobiledna/data/hg19_rmsk_uniq.bed")
expTEs = readBed("/data2/han_lab/mobiledna/data/TET/maxgt5.bed")
expuniqnoovpTEs = readBed("/data2/han_lab/mobiledna/result.rsem.TET.instance/gene.noovp.uniq.TEs.bed")
locResults = runLOLA(userSets, expuniqnoovpTEs, regionDB, cores=1)
locResults

#coreDB = regionDB = loadRegionDB("LOLACore/hg19")
coreDB = loadRegionDB("/data2/han_lab/database/LOLA/LOLACore/hg19")
testresults = runLOLA(userSets, expuniqnoovpTEs, coreDB, cores=1)
testresults
testresults[order(meanRnk, decreasing=FALSE),][1:20,]

