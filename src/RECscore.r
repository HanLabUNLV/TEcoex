


######## from Sophia #########################

#rank function
rrfunction<- function(x,L)
{
  (x/L)-(1/(2*L))
}

#inverted rank
invrrfunction<- function(x,L)
{
  invr <- L-x+1
  if(any(invr<0, na.rm=TRUE)) {
    print(x)
    print(L)
    print(invr)
  }
  return(invr)
}


#HKNOTS
testpos<- function(r, L)
{
  rrw <- rrfunction(as.numeric(r),L)
  Hknotresult <- -2*sum(log(rrw), na.rm=TRUE)
  return(Hknotresult)
}

testneg<- function(r, L)
{
  invr <- invrrfunction(as.numeric(r),L)
  invrrw<- rrfunction(invr,L)
  Hknotinvresult<- -2*sum(log(invrrw), na.rm=TRUE)
  return(Hknotinvresult)
}

#################################


combine_ranks <- function(coefs_list) {
  ranklists = list()
  print (coefs_list)
  for (j in 1:length(coefs_list)) {
    coefs <- coefs_list[[j]]
    coefs <- coefs[order(-coefs[,"partialeta2"]),]
    ranks <- rank(coefs[,"partialeta2"], na.last=TRUE, ties.method=c("min"))
    ranks <- data.frame(cbind(rownames(coefs), ranks), row.names=rownames(coefs))
    ranks[,2] <- as.numeric(as.character(ranks[,2]))
    colnames(ranks) <- c("genenames", "rank")
    ranklists[[j]] <- ranks
  }
  
  allranks <- ranklists[[1]]
  for (i in 2:length(ranklists)) {
    allranks <- merge(allranks,ranklists[[i]], by="genenames", all=TRUE)
  }
  allranks
}



RECscore <- function(pminus, pplus) {
  REC=NULL
  if(pminus < pplus) {
    REC <- log(2*pminus, base=10)
  } else if (pminus > pplus) {
    REC <- log(2*pplus, base=10)
  } else {
    REC = 0
  }
  REC
}


calculate_REC <- function(genenames, generanks) {
  return_REC = data.frame(genenames, rep(NA,nrow(generanks)))
  L <- colSums(!is.na(generanks))
  print(L)
  for (i in 1:nrow(generanks)) {
    ranks <- generanks[i,]
    n = sum(!is.na(ranks))
    if (n < 5) { next; }
    print(i)
    Hknot=testpos(ranks, L)
    Hknotinv=testneg(ranks, L)
    tailedH<-rbind(Hknot, Hknotinv)
    df = n
    pminus<-1-pchisq(Hknot, df = 2*n)
    pplus<-1-pchisq(Hknotinv, df = 2*n)
    REC <- RECscore(pminus, pplus)
    return_REC[i,2] = REC
  }
  return_REC <- return_REC[order(return_REC[,2]),]
  return(return_REC)
}



save_RECscores <- function(result_dir) {

  cancertypes = list("BLCA", "BRCA", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), "LIHC",  c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC")
  cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))
  
  coefs_plus = list()
  coefs_minus = list()
  for (j in 1:length(cancertypes)) {
    type <- cancertypes[j]
    cancerdir <- cancerdirs[j]
    fname = paste(result_dir,cancerdir, "results_lm.txt", sep="/")
    
    if (file.exists(fname)) {
    print (fname)
      eval <- read.table(fname, sep="\t", header=TRUE, row.names = 1)
      eval <- eval[order(-eval$partialeta2),]
      minuscoef <- eval[eval$generatio_coef<0,]
      pluscoef <- eval[eval$generatio_coef>0,]
      coefs_plus[[j]]<- pluscoef[,c("generatio_coef",	"generatio_pval", "partialeta2")]
      coefs_minus[[j]]<- minuscoef[,c("generatio_coef",	"generatio_pval", "partialeta2")]
    }
  }
  
  
  minus_ranks <- combine_ranks(coefs_minus)
  minus_REC <- calculate_REC(minus_ranks[,1], minus_ranks[,-1])
  colnames(minus_REC) = c("gene", "RECscore")
  filename = paste(result_dir, "REC.minus.new.txt", sep="/")
  write.table(minus_REC, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  plus_ranks <- combine_ranks(coefs_plus)
  plus_REC <- calculate_REC(plus_ranks[,1], plus_ranks[,-1])
  colnames(plus_REC) = c("gene", "RECscore")
  filename = paste(result_dir, "REC.plus.new.txt", sep="/")
  write.table(plus_REC, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
}

#####################################################
args = commandArgs(trailingOnly=TRUE)
resultdir = "../result.rsem/gene2L1HS"
print (args)
if (length(args)>0) {
  resultdir = args[1]
  if (substr(resultdir, nchar(resultdir), nchar(resultdir)) != "/") {
    resultdir = paste0(resultdir, "/")
  }
}
print(resultdir)
save_RECscores(resultdir)
#resultdir <- "../result.rsem.TET.uniq/oldLINEcov"
#save_RECscores(resultdir)
