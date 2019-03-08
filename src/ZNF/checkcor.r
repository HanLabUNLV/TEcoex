

args = commandArgs(trailingOnly = TRUE)
resultdir = args[1]
TEname = args[2]
ZNFname = args[3]
tissue = tail(args,n=-2)

if("COADREAD" %in% tissue) {
  tissue = c("COAD", "READ")
} else if ("ESCASTAD" %in% tissue) {
  tissue = c("ESCA", "STAD")
} else if ("KICHKIRCKIRP" %in% tissue) {
  tissue = c("KICH", "KIRC", "KIRP")
} else if ("LUADLUSC" %in% tissue) {
  tissue = c("LUAD", "LUSC")
}  

print(resultdir)
print(TEname)
print(ZNFname)
print(tissue)


setwd("../../")
ZNFfilename = paste0(resultdir, "/VSTcnts/", ZNFname, ".txt")
TEfilename = paste0(resultdir, ".instance/VSTcnts/", TEname, ".txt")
print(ZNFfilename)
print(TEfilename)

ZNF <- read.table(ZNFfilename, header=TRUE)
TE <- read.table(TEfilename, header = TRUE)
patient <- read.table("data/patient.info", sep="\t", header=TRUE)
dat <- cbind(patient$tissue, ZNF, TE)
b_tissue = patient$tissue %in% tissue
subdat <- dat[b_tissue,]

print(subdat)
pdf(paste0("result.rsem.TET.instance/ZNF/plots/", paste(TEname, ZNFname, tissue, "pdf",  sep=".")))
plot(subdat[,3], subdat[,5])
dev.off()
