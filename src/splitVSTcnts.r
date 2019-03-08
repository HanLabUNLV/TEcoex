
resultdir = "../result.rsem.TET.instance/"
args = commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
}


significantloci = read.table(paste0(resultdir, "ZNF/significant_instances.noovp100000.txt"), header=FALSE, stringsAsFactors=FALSE)
normalcnts <- read.table(paste0(resultdir, "VSTcnt.txt"), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
gnames <- rownames(normalcnts)
TEidx <- grepl (":", gnames)
geneidx <-!TEidx
TEnames = gnames[TEidx]
TEnames = substr(TEnames,1,nchar(TEnames)-2)
rownames(normalcnts)[TEidx] <- TEnames
genecnts <- normalcnts[geneidx,]
sigTEs <- unique(significantloci[,1])
sigVSTcnts <- normalcnts[sigTEs,]
newcnts <- rbind.data.frame(genecnts, sigVSTcnts)
newcnts <- newcnts[!is.na(newcnts[,1]),]
cntnames = colnames(newcnts)


scnts <- newcnts[,(grepl("BLCA", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.BLCA.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("BRCA", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.BRCA.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("COAD", cntnames))|(grepl("READ", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.COADREAD.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("ESCA", cntnames))|(grepl("STAD", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.ESCASTAD.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("HNSC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.HNSC.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("KICH", cntnames))|(grepl("KIRC", cntnames))|(grepl("KIRP", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.KICHKIRCKIRP.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("LIHC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.LIHC.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("LUAD", cntnames))|(grepl("LUSC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.LUADLUSC.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("PRAD", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.PRAD.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("THCA", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.THCA.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("UCEC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.UCEC.sig.txt"), quote=FALSE, row.names=TRUE, sep="\t")
