L1HS5prime <- read.table("../data/TET/L1HS5prime.300.txt", row.names=1)
L1HS5prime.bowtie <- read.table("../data/TET.bowtie//L1HS5prime.300.txt", row.names=1)
L1HS5prime.comp <- merge(L1HS5prime, L1HS5prime.bowtie, by="row.names")
L1HS <- L1HS5prime.comp[grepl("TCGA", L1HS5prime.comp[,1]),]

patient.info <- read.table("../data/patient.info", header=TRUE, sep="\t")

rownames(L1HS) = L1HS[,1]
L1HS <- L1HS[patient.info[,1],]
L1HS <- cbind.data.frame(patient.info[,c(1:5, 9:10)], L1HS)
L1HS$age <- 2019 - as.numeric(as.character(L1HS$year_of_birth))

L1HS2 <- L1HS[L1HS$tissue != "BRCA",]
L1HS2 <- L1HS2[L1HS2$tissue != "PRAD",]



