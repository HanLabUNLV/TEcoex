
resultdir = "../result.rsem.TET.instance/"
b_noovp = ""
args = commandArgs(trailingOnly=TRUE)
if (length(args)>=1) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
  b_noovp = args[2]
}

print(resultdir)
print(b_noovp)

significantloci = read.table(paste0(resultdir, "ZNF/significant_instances.noovp100000.txt"), header=FALSE, stringsAsFactors=FALSE)
normalcnts <- read.table(paste0(resultdir, "VSTcnt.txt"), header=TRUE, row.names=1, sep="\t", check.names = FALSE)
gnames <- rownames(normalcnts)
TEidx <- grepl (":", gnames)
geneidx <-!TEidx
TEnames = gnames[TEidx]
TEnames = substr(TEnames,1,nchar(TEnames)-2)
rownames(normalcnts)[TEidx] <- TEnames
genecnts <- normalcnts[geneidx,]



if (grepl("instance", resultdir) && (! is.na(b_noovp)) && (b_noovp == "noovp")) {
  filters = read.table(paste0(resultdir, "gene.noovp.uniq.TEs.bed"), sep="\t")
  filternames = as.character(filters[,4])
  TEcnts <- normalcnts[TEidx,]
  TEnames <- rownames(TEcnts)
  TElocus = matrix(unlist(strsplit(TEnames, ":")), ncol=4, byrow=TRUE)[,1]
  rownames(TEcnts) = TElocus
  names(TEnames) = TElocus
  filteredTEcnts <- TEcnts[filternames,]
  filteredTEnames <- TEnames[filternames]
  rownames(filteredTEcnts) = filteredTEnames
  write.table(filteredTEcnts, paste0(resultdir, "VSTcnt.filtered.noovp.uniq.txt"), quote=FALSE, row.names=TRUE, sep="\t")
#  sigVSTcnts <- filteredTEcnts[sigTEs,]
  sigVSTcnts <- filteredTEcnts 
  sigVSTcnts <- sigVSTcnts[!is.na(sigVSTcnts[,1]),]
} else {
  sigTEs <- unique(significantloci[,1])
  sigVSTcnts <- normalcnts[sigTEs,]
  sigVSTcnts <- sigVSTcnts[!is.na(sigVSTcnts[,1]),]
}

newcnts <- rbind.data.frame(genecnts, sigVSTcnts)
newcnts <- newcnts[!is.na(newcnts[,1]),]
cntnames = colnames(newcnts)


scnts <- newcnts[,(grepl("BLCA", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.BLCA.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("BRCA", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.BRCA.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("COAD", cntnames))|(grepl("READ", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.COADREAD.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("ESCA", cntnames))|(grepl("STAD", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.ESCASTAD.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("HNSC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.HNSC.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("KICH", cntnames))|(grepl("KIRC", cntnames))|(grepl("KIRP", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.KICHKIRCKIRP.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("LIHC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.LIHC.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("LUAD", cntnames))|(grepl("LUSC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.LUADLUSC.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("PRAD", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.PRAD.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("THCA", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.THCA.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
scnts <- newcnts[,(grepl("UCEC", cntnames))]
write.table(scnts, paste0(resultdir, "VSTcnt.UCEC.sig", b_noovp, ".txt"), quote=FALSE, row.names=TRUE, sep="\t")
