library(foreach)
library(doParallel)

#TCGA-ZU-A8S4-11A-60e3b4d8-081a-4247-8144-3cce35345dbf.ele.cntTable 

countTEs <- function(rawfile, discountfile, by, TEnamecols) {

  print (paste(rawfile, discountfile))
  rawcnts <- read.table(paste("../data/TET.bowtie/",as.character(rawfile), sep="/"), header=TRUE)
  discounts <- read.table(paste("../data/TET.bowtie/",as.character(discountfile), sep="/"), header=TRUE)
  #m <- matrix(unlist(strsplit(as.character(discounts[, 1]), ":")), ncol = TEnamecols, byrow = TRUE)
  if (nrow(discounts) != nrow(rawcnts)) {
    message("rownums do not match")
    rawcnts <- read.table(paste("../data/TET.bowtie/",as.character(rawfile), sep="/"), header=FALSE)
    #stop()
  }

  TEcnts <- cbind.data.frame(rawcnts[,2], discounts[, 2])
  TEcnts$diff <- TEcnts[,1]-TEcnts[,2]
  if (by == "ele") {
    #colnames(TEcnts) = c("TEname", "TEsubfamily", "TEfamily", "rawcnt", "discount", "diff")
  } else if (by == "instance") {
    #colnames(TEcnts) = c("TEintegrant", "TEname", "TEsubfamily", "TEfamily", "strand", "rawcnt", "discount", "diff")
  }
  colnames(TEcnts) = c("rawcnt", "discount", "diff")
  return (TEcnts)
} 

args = commandArgs(trailingOnly=TRUE)
rawfnames = args[1]
discountfnames = args[2]
print (paste(rawfnames, discountfnames))


  #cntby = "ele"
  cntby = "instance"
  TEnamecols = 0 
  if (cntby == "ele") {
    TEnamecols = 3
  }
  if (cntby == "instance") {
    TEnamecols = 5
  }

  #rawfile.names <- dir("../data/TET.bowtie/", pattern =paste0(".", cntby, ".cntTable2"), recursive=TRUE)
  rawfile.names <- read.table(rawfnames, stringsAsFactors = FALSE)
  #discountfile.names <- dir("../data/TET.bowtie/", pattern =paste0(".", cntby, ".cntTable"), recursive=TRUE)
  discountfile.names <- read.table(discountfnames, stringsAsFactors = FALSE)
  x <- grepl("discount", discountfile.names)
  discountidx <- x
  rawfiles <- rawfile.names[[1]]
  patientIDs <- substr(unlist(strsplit(rawfiles, "/"))[c(FALSE, TRUE)], 9, 12)
  rawfiles <- cbind.data.frame(patientIDs, rawfiles)
  colnames(rawfiles) <- c("patient", "filename")
  discountfiles <- discountfile.names[[1]][discountidx]
  patientIDs <- substr(unlist(strsplit(discountfiles, "/"))[c(FALSE, TRUE)], 9, 12)
  discountfiles <- cbind.data.frame(patientIDs, discountfiles)
  colnames(discountfiles) <- c("patient", "filename")
  stopifnot(as.character(rawfiles[,1]) == as.character(discountfiles[,1]))
  TEnames = read.table("../data/TEinstancenames.txt", header=TRUE, sep="\t")




  Maxdiff = as.data.frame(matrix(0, nrow=nrow(rawfiles), ncol=4*5))
  L1HSdiff = as.data.frame(matrix(0, nrow=nrow(rawfiles), ncol=4*5))

  for (i in 1:nrow(rawfiles)) {
  print (i)
    TEcnts <- cbind(TEnames, countTEs(rawfiles[i,2], discountfiles[i,2], cntby, TEnamecols))
    bylargediff = TEcnts[order(-TEcnts$diff),]
    L1HScnts <- TEcnts[TEcnts$TEname == "L1HS",]
    byL1HSdiff = L1HScnts[order(-L1HScnts$diff),]
    for (j in 0:4) {
      Maxdiff[i,1+(j*4)] = as.character(bylargediff[j+1,1])
      Maxdiff[i,2:4+(j*4)] = bylargediff[j+1,1:3+TEnamecols]
      L1HSdiff[i,1+(j*4)] = as.character(byL1HSdiff[j+1,1])
      L1HSdiff[i,2:4+(j*4)] = byL1HSdiff[j+1,1:3+TEnamecols]
    }
  }

  Maxdiff <- cbind.data.frame(rawfiles$patient,Maxdiff)
  L1HSdiff <- cbind.data.frame(rawfiles$patient,L1HSdiff)
  colnames(Maxdiff) = c("patient", rep(c("TEname", "rawcnt", "discount", "diff"), 5))
  colnames(L1HSdiff)= c("patient", rep(c("TEname", "rawcnt", "discount", "diff"), 5))


  #write.table(Maxdiff, file=paste("maxdiff",cntby,"txt", sep="." ), quote=FALSE, row.names = FALSE)
  #write.table(L1HSdiff, file=paste("L1HSdiff",cntby,"txt", sep="." ), quote=FALSE, row.names = FALSE)


  patient2cancer <- read.table("../data/patient.info", header=TRUE, sep="\t")
  rownames(patient2cancer) = substr(patient2cancer[,1], 9, 12)
  Maxdiff <- cbind.data.frame(Maxdiff, patient2cancer[Maxdiff$patient,"tissue"])
  L1HSdiff <- cbind.data.frame(L1HSdiff, patient2cancer[L1HSdiff$patient,"tissue"])

  write.table(Maxdiff, file=paste("maxdiff",cntby,rawfnames,"txt", sep="." ), quote=FALSE, row.names = FALSE)
  write.table(L1HSdiff, file=paste("L1HSdiff",cntby,rawfnames,"txt", sep="." ), quote=FALSE, row.names = FALSE)
