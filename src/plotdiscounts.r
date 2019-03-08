library("ggplot2")
library("fields")


args = commandArgs(trailingOnly=TRUE)
ele_name="L1HS"
print (args)
if (length(args)>0) {
  ele_name = args[1]
}
print(ele_name)



raw_ele <- read.table(paste0("../result.rsem.TET/", ele_name, ".rawcnts.txt"), sep="\t")
patient <- substr(unlist(strsplit(as.character(raw_ele[,1]), "/"))[c(FALSE, FALSE, FALSE, FALSE, TRUE)], 9, 12)
rownames(raw_ele) = patient
discount_ele <- read.table(paste0("../result.rsem.TET/", ele_name, ".discount.txt"), sep="\t")
patient <- substr(unlist(strsplit(as.character(discount_ele[,1]), "/"))[c(FALSE, FALSE, FALSE, FALSE, TRUE)], 9, 12)
rownames(discount_ele) = patient
p2c <- read.table("../data/patient.info", header=TRUE, sep="\t")
rownames(p2c) <- substr(p2c$patient, 9, 12)

r2d <- cbind.data.frame(rownames(raw_ele), p2c[rownames(raw_ele),2], log2(raw_ele[rownames(raw_ele),2]), log2(discount_ele[rownames(raw_ele),2]))
colnames(r2d) = c("patient", "tissue", "raw", "discount")
rownames(r2d) = r2d$patient



uniq_raw_ele <- read.table(paste0("../result.rsem.TET.uniq/", ele_name, ".rawcnts.txt"), sep="\t")
patient <- substr(unlist(strsplit(as.character(uniq_raw_ele[,1]), "/"))[c(FALSE, FALSE, FALSE, FALSE, TRUE)], 9, 12)
rownames(uniq_raw_ele) = patient
uniq_discount_ele <- read.table(paste0("../result.rsem.TET.uniq/", ele_name, ".discount.txt"), sep="\t")
patient <- substr(unlist(strsplit(as.character(uniq_discount_ele[,1]), "/"))[c(FALSE, FALSE, FALSE, FALSE, TRUE)], 9, 12)
rownames(uniq_discount_ele) = patient

uniqpatient <- as.character(patient)
r2d <- cbind.data.frame(r2d[uniqpatient,], log2(uniq_raw_ele[,2]), log2(uniq_discount_ele[,2]))
colnames(r2d) = c("patient", "tissue", "raw", "discount", "uniq.raw", "uniq.discount")



bowtie_raw_ele <- read.table(paste0("../result.rsem.TET.bowtie/", ele_name, ".rawcnts.txt"), sep="\t")
patient <- substr(unlist(strsplit(as.character(bowtie_raw_ele[,1]), "/"))[c(FALSE, FALSE, FALSE, FALSE, TRUE)], 9, 12)
rownames(bowtie_raw_ele) = patient
bowtie_discount_ele <- read.table(paste0("../result.rsem.TET.bowtie/", ele_name, ".discount.txt"), sep="\t")
patient <- substr(unlist(strsplit(as.character(bowtie_discount_ele[,1]), "/"))[c(FALSE, FALSE, FALSE, FALSE, TRUE)], 9, 12)
rownames(bowtie_discount_ele) = patient

bowtiepatient <- as.character(patient)
r2d <- cbind.data.frame(r2d[bowtiepatient,], log2(bowtie_raw_ele[,2]), log2(bowtie_discount_ele[,2]))
colnames(r2d) = c("patient", "tissue", "raw", "discount", "uniq.raw", "uniq.discount", "bowtie.raw", "bowtie.discount")


cancertypes = c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

r2ddata <- r2d[r2d$tissue %in% cancertypes,]
r2ddata$tissue <- droplevels(r2ddata$tissue)

myColors = c("#e6194b", "#3cb44b", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#ffe119", "#000000")
names(myColors) <- levels(r2ddata$tissue)
colScale <- scale_colour_manual(name = "tissue",values = myColors)

# Add density for each point
p = ggplot(r2ddata,  aes(raw, discount, color=tissue
  #, alpha = 1/density
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8))  +
  theme_minimal()
ggsave(paste0(ele_name, ".raw2discount.pdf"))


p = ggplot(r2ddata,  aes(uniq.raw, uniq.discount, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name,".uniqraw2uniqdiscount.pdf"))



p = ggplot(r2ddata,  aes(raw, uniq.raw, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name, ".raw2uniqraw.pdf"))



p = ggplot(r2ddata,  aes(discount, uniq.discount, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name, ".discount2uniqdiscount.pdf"))



p = ggplot(r2ddata,  aes(bowtie.raw, bowtie.discount, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name,".bowtieraw2bowtiediscount.pdf"))


p = ggplot(r2ddata,  aes(raw, bowtie.raw, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name, ".raw2bowtieraw.pdf"))



p = ggplot(r2ddata,  aes(discount, bowtie.discount, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name, ".discount2bowtiediscount.pdf"))



p = ggplot(r2ddata,  aes(uniq.raw, bowtie.raw, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name, ".uniqraw2bowtieraw.pdf"))



p = ggplot(r2ddata,  aes(uniq.discount, bowtie.discount, color=tissue
  #, alpha=1/density4
  ))
p + geom_point(shape=16, size=2) + colScale + 
#   scale_alpha(range = c(.25, .8)) + 
  theme_minimal()
ggsave(paste0(ele_name, ".uniqdiscount2bowtiediscount.pdf"))



