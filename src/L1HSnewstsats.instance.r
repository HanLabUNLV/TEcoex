
library(ggplot2)

rawcnts <- read.table("../data/TET/rawcnts.new.maxgt5.txt", header=TRUE, row.names=1)
L1HSpos <- read.table("../data/TET/maxgt5.bed")
L1HSidx <- grepl("L1HS", rownames(rawcnts))
L1HSposidx <- grepl("L1HS", L1HSpos[,4])
L1HScnts <- rawcnts[L1HSidx,]
dim(L1HScnts)
L1HSpos <- L1HSpos[L1HSposidx,]
L1HScntsnames <- unlist(lapply(strsplit(rownames(L1HScnts), ":"), "[[", 1)) 
rownames(L1HScnts) <- L1HScntsnames
L1HScntsreordered <- L1HScnts[as.character(L1HSpos[,4]),]
L1HSsums <- rowSums(L1HScntsreordered[,-1])
L1HSmax <- apply(L1HScntsreordered[,-1], 1, max)
L1HS <- cbind(L1HSpos, log2(L1HSsums), log2(L1HSmax))
L1HSmap <- read.table("../data/L1HS.rmsk.mappability.uniqpercent", header=FALSE)
rownames(L1HSmap) <- L1HSmap[,13]
L1HS <- cbind(L1HS, L1HSmap[as.character(L1HSpos[,4]),21])
L1HS <- cbind(L1HS, (L1HS[,3]-L1HS[,2])*L1HS[,9])
colnames(L1HS)[7:10] = c("L1HSsum", "L1HSmax", "uniqper", "uniqlen")

ggplot(data=L1HS, aes(x=uniqlen, y=L1HSsum)) + geom_point(alpha = 0.4) + theme_bw()
ggsave("mappability2cnt.png")
ggplot(data=L1HS, aes(x=uniqlen, y=L1HSmax)) + geom_point(alpha = 0.4) + theme_bw()
ggsave("mappability2max.png")


library (reshape)
bowtiecnts <- read.table("../data/TET.bowtie/L1HS.table")
bowtiecnts <- bowtiecnts[grepl("TCGA", bowtiecnts[,2]),]
patient <- gsub("-", ".", as.character(bowtiecnts[,2]))
bowtiecnts <- cbind.data.frame(bowtiecnts, patient)
md <- melt(bowtiecnts, id=(c("V3", "patient")), measure.vars="V8")
newcnt <- cast(md, V3 ~ patient)
rownames(newcnt) <- newcnt[,1]
newcnt2 <- newcnt[as.character(L1HS[,4]),]
bowtiesums <- rowSums(newcnt2[,-1])
bowtiemax <- apply(newcnt2[,-1], 1, max)
L1HS <- cbind.data.frame(L1HS, log2(bowtiesums), log2(bowtiemax))
colnames(L1HS)[11:12] = c("bowtiesum", "bowtiemax")

ggplot(data=L1HS, aes(x=uniqlen, y=bowtiesum)) + geom_point(alpha = 0.4) + theme_bw()
ggsave("mappability2cnt.bowtie.png")
ggplot(data=L1HS, aes(x=uniqlen, y=bowtiemax)) + geom_point(alpha = 0.4) + theme_bw()
ggsave("mappability2max.bowtie.png")

ggplot(data=L1HS, aes(x=bowtiesum, y=L1HSsum)) + geom_point(alpha = 0.4) + theme_bw() 
ggsave("L1HS.instance.sum.bowtie2star.png")
ggplot(data=L1HS, aes(x=bowtiemax, y=L1HSmax)) + geom_point(alpha = 0.4) + theme_bw() 
ggsave("L1HS.instance.max.bowtie2star.png")

L1HS$bowtiemax[L1HS$bowtiemax==-Inf] = -1
L1HS$bowtiesum[L1HS$bowtiesum==-Inf] = -1
test <- lm(L1HSmax ~ uniqlen, L1HS)
print(anova(test))
#Response: L1HSmax
#           Df Sum Sq Mean Sq F value  Pr(>F)  
#uniqlen     1   2.76 2.76431  2.9692 0.0.437 .
#Residuals 608 566.05 0.93101              
#
test <- lm(L1HSsum ~ uniqlen, L1HS)
print(anova(test))
test2 <- lm(bowtiemax ~ uniqlen, L1HS)
print(anova(test2))
#Response: bowtiemax
#           Df Sum Sq Mean Sq F value    Pr(>F)    
#uniqlen     1  75.55  75.549  57.867 1.072e-13 ***
#Residuals 608 793.78   1.306  
#
test2 <- lm(bowtiesum ~ uniqlen, L1HS)
print(anova(test2))





