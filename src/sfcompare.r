sf1old <- read.table("/data2/han_lab/mobiledna.oldresult/rsem/sizefactors.1.txt", sep="\t", row.names=1)
sf1new <- read.table("/data2/han_lab/mobiledna/result.rsem/sizefactors.1.txt", sep="\t")
sf1new <- cbind.data.frame(matrix(unlist(strsplit(matrix(unlist(strsplit(rownames(sf1new), "-")), ncol=3, byrow=TRUE)[,3],"[.]")), ncol=2, byrow=TRUE), sf1new)
colnames(sf1new) <- c("patient", "tissue", "sf")
sf1old <- cbind.data.frame(rownames(sf1old), sf1old)
colnames(sf1old) <- c("patient", "sf")
sf1old$patient = substr(sf1old$patient, 1, 4)
sf1comp <- merge(sf1old, sf1new, by="patient")
#plot(sf1comp$sf.y, sf1comp$sf.x)
cor.test(sf1comp$sf.x, sf1comp$sf.y)


sf2old <- read.table("/data2/han_lab/mobiledna.oldresult/rsem/sizefactors.2.txt", sep="\t", row.names=1)
sf2new <- read.table("/data2/han_lab/mobiledna/result.rsem/sizefactors.2.txt", sep="\t")
sf2new <- cbind.data.frame(matrix(unlist(strsplit(matrix(unlist(strsplit(rownames(sf2new), "-")), ncol=3, byrow=TRUE)[,3],"[.]")), ncol=2, byrow=TRUE), sf2new)
colnames(sf2new) <- c("patient", "tissue", "sf")
sf2old <- cbind.data.frame(rownames(sf2old), sf2old)
colnames(sf2old) <- c("patient", "sf")
sf2old$patient = substr(sf2old$patient, 1, 4)
sf2comp <- merge(sf2old, sf2new, by="patient")
#plot(sf2comp$sf.y, sf2comp$sf.x)
cor.test(sf2comp$sf.x, sf2comp$sf.y)
test <- merge(sf2old, sf2new, by="patient", all.x=TRUE, incomparables=NA)
missing <- test[is.na(test$tissue),]
patient2tissue <- read.table("Work/TE-pre/L1HSstats.DEGES/patientToCancerMappings.tsv", sep="\t")
patient2tissue <- cbind.data.frame(matrix(unlist(strsplit(as.character(patient2tissue$V1) , "-")), ncol=3, byrow=TRUE)[,3], patient2tissue)
colnames(patient2tissue)[1] <- "patient"
missingtissues <- merge(missing, patient2tissue, by="patient")

