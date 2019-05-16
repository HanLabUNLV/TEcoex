
chip <- read.table("../../data/chipTEenrichment.txt", sep="\t", header=TRUE, row.names=1)
coex <- read.table("../../result.rsem.TET/ZNF.noctrl/SigResults.txt", sep="\t", header=TRUE)

#TEsbound <- rownames(chip)
#TEsboundlist <- strsplit(TEsbound, "[.]")
#n.obs <- sapply(TEsboundlist, length)
#seq.max <- seq_len(max(n.obs))
#TEsboundmat <- t(sapply(TEsboundlist, "[", i = seq.max))
#TEsboundmat[is.na(TEsboundmat[,3]),3] = TEsboundmat[is.na(TEsboundmat[,3]),2]
#TEsboundmat <- TEsboundmat[,-1]
#TEsboundmat <- TEsboundmat[,-1]
#TEsboundmat[!is.na(TEsboundmat[,3]),2] = paste(TEsboundmat[!is.na(TEsboundmat[,3]),2], TEsboundmat[!is.na(TEsboundmat[,3]),3], sep="-")
#TEsboundmat <- TEsboundmat[,-3]
#TEsboundmat[TEsboundmat[,1]=="MaLR",1] = TEsboundmat[TEsboundmat[,1]=="MaLR",2]
#TEsboundmat[TEsboundmat[,1]=="Tigger",1] = TEsboundmat[TEsboundmat[,1]=="Tigger",2]
#TEsboundmat[TEsboundmat[,1]=="Charlie",1] = TEsboundmat[TEsboundmat[,1]=="Charlie",2]
#TEsboundmat[TEsboundmat[,1]=="Blackjack",1] = TEsboundmat[TEsboundmat[,1]=="Blackjack",2]
#TEsboundmat[TEsboundmat[,1]=="Tip100",1] = TEsboundmat[TEsboundmat[,1]=="Tip100",2]
#TEsboundmat[TEsboundmat[,1]=="Mariner",1] = TEsboundmat[TEsboundmat[,1]=="Mariner",2]
#TEsboundmat[TEsboundmat[,1]==TEsboundmat[,2],2] = NA
#TEsboundmat[!is.na(TEsboundmat[,2]),1] = paste(TEsboundmat[!is.na(TEsboundmat[,2]),1], TEsboundmat[!is.na(TEsboundmat[,2]),2], sep="-")
#rownames(chip) <- TEsboundmat[,1]
#
TEscoex <- as.character(coex[,1])
TEscoexlist <- strsplit(TEscoex, ":")
n.obs <- sapply(TEscoexlist, length)
seq.max <- seq_len(max(n.obs))
TEscoexmat <- t(sapply(TEscoexlist, "[", i = seq.max))
TEscoexmat <- TEscoexmat[,1]
coex = cbind.data.frame(TEscoexmat, coex)
coex[,1] = as.character(coex[,1])

TEscoexnames <- unique(TEscoexmat)
TEscoexnames <- TEscoexnames[order(TEscoexnames)]

joincoex = c()
joinchip = c()

ChIPZFPs = c()

for (i in 1:nrow(chip)) {
  TEname = rownames(chip)[i]
  coexpressed = coex[coex[,1]==TEname,]
  if (nrow(coexpressed) == 0) {
    next;
  }
  coexZFPs <-  matrix(unlist(strsplit(as.character(coexpressed[,3]), "[|]")), ncol=2, byrow=TRUE)
  coexZFPs <- coexZFPs[,1]
  print(coexZFPs)
  
  chipZFPs <- colnames(chip)[which(chip[TEname,]>20)]
  print(chipZFPs)
  for (j in 1:length(chipZFPs)) {
    ChIPZFPs = rbind(ChIPZFPs, c(TEname, chipZFPs[j]))
  }

  ZFPset <- intersect(coexZFPs, chipZFPs)
  if (length(ZFPset) > 0) { 
    joinchip <- rbind(joinchip, chip[TEname,])
    coexmatch <- apply(sapply(ZFPset, grepl, coexpressed[,3]), 1, any)
    joincoex <- rbind(joincoex, coexpressed[coexmatch,])
  } 

  
}

write.table(ChIPZFPs, "../../data/chipTEenrichment.pairs.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(joincoex, "../../result.rsem.TET/ZNF.noctrl/join.coex.txt", sep="\t", quote=FALSE)
write.table(joinchip, "../../result.rsem.TET/ZNF.noctrl/join.chip.txt", sep="\t", quote=FALSE)

