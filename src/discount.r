#TCGA-ZU-A8S4-11A-60e3b4d8-081a-4247-8144-3cce35345dbf.ele.cntTable 

countTEs <- function(rawfile, discountfile, by, TEnamecols) {

  rawcnts <- read.table(paste("../data/TET/",as.character(rawfile), sep="/"), header=FALSE)
  m <- matrix(unlist(strsplit(as.character(rawcnts[, 1]), ":")), ncol = TEnamecols, byrow = TRUE)
  TEcnts <- cbind.data.frame(m, rawcnts[, 2])
  

  discounts <- read.table(paste("../data/TET/",as.character(discountfile), sep="/"), header=TRUE)
  m <- matrix(unlist(strsplit(as.character(discounts[, 1]), ":")), ncol = TEnamecols, byrow = TRUE)
  if (any(m[,1]!=TEcnts[,1])) {
    stop()
  }
  TEcnts <- cbind.data.frame(TEcnts, discounts[, 2])
  TEcnts$diff <- TEcnts[,TEnamecols+1]-TEcnts[,TEnamecols+2]
  if (by == "ele") {
    colnames(TEcnts) = c("TEname", "TEsubfamily", "TEfamily", "rawcnt", "discount", "diff")
  } else if (by == "instance") {
    colnames(TEcnts) = c("TEintegrant", "TEname", "TEsubfamily", "TEfamily", "strand", "rawcnt", "discount", "diff")
  }

  return (TEcnts)
} 


  #cntby = "ele"
  cntby = "instance"
  TEnamecols = 0 
  if (cntby == "ele") {
    TEnamecols = 3
  }
  if (cntby == "instance") {
    TEnamecols = 5
  }

  rawfile.names <- dir("../data/TET/", pattern =paste0(".", cntby, ".cntTable2"), recursive=TRUE)
  discountfile.names <- dir("../data/TET/", pattern =paste0(".", cntby, ".cntTable"), recursive=TRUE)
  x <- grepl("discount", discountfile.names)
  discountidx <- x
  rawfiles <- rawfile.names
  patientIDs <- substr(unlist(strsplit(rawfiles, "/"))[c(FALSE, TRUE)], 9, 12)
  rawfiles <- cbind.data.frame(patientIDs, rawfiles)
  colnames(rawfiles) <- c("patient", "filename")
  discountfiles <- discountfile.names[discountidx]
  patientIDs <- substr(unlist(strsplit(discountfiles, "/"))[c(FALSE, TRUE)], 9, 12)
  discountfiles <- cbind.data.frame(patientIDs, discountfiles)
  colnames(discountfiles) <- c("patient", "filename")
  stopifnot(rawfiles[,1] == discountfiles[,1])


  Maxdiff = as.data.frame(matrix(0, nrow=nrow(rawfiles), ncol=4*5))
  L1HSdiff = as.data.frame(matrix(0, nrow=nrow(rawfiles), ncol=4*5))

  for (i in 1:nrow(rawfiles)) {
  print (i)
    TEcnts <- countTEs(rawfiles[i,2], discountfiles[i,2], cntby, TEnamecols)
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
  colnames(Maxdiff) = c("patient", "TEname", "rawcnt", "discount", "diff")
  colnames(L1HSdiff)= c("patient", "TEname", "rawcnt", "discount", "diff")

  patient2cancer <- read.table("../data/patient.info", header=TRUE)
  rownames(patient2cancer) = substr(patient2cancer[,1], 9, 12)
  Maxdiff <- cbind.data.frame(Maxdiff, patient2cancer[Maxdiff$patient,"tissue"])
  L1HSdiff <- cbind.data.frame(L1HSdiff, patient2cancer[L1HSdiff$patient,"tissue"])

  write.table(Maxdiff, file=paste("maxdiff",cntby,"txt", sep="." ), quote=FALSE)
  write.table(L1HSdiff, file=paste("L1HSdiff",cntby,"txt", sep="." ), quote=FALSE)
