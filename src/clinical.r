library(survival)
library(ggplot2)


L1HS <- read.table("../result.rsem.TET/L1HS5prime.VST.txt", sep="\t", header=TRUE)
clinical <- L1HS
#clinical <- L1HS[,c(2,3,4, 19, 20, 22, 24, 31, 1)]
ME1 <- read.table("../result.rsem.TET/WGCNA/consensus/eigengeneME1.txt", header=FALSE, sep="\t")
rownames(ME1) = ME1[,1]
ME1 <- ME1[as.character(clinical$patient),]
all_clin  <- cbind.data.frame(clinical, ME1[,2])
#colnames(all_clin)[c(7,6,8, 9, 10)] <- c('new_tumor_days', 'death_days', 'followUp_days', 'L1HS5prime', 'ME1')
radiation <- read.table("../data/THCA.radiation.txt")
all_clin$radiation = as.character(all_clin$treatment_or_therapy)
rownames(radiation) = radiation[,1]
all_clin[all_clin$tissue=="THCA","radiation"]  = as.character(radiation[rownames(all_clin[all_clin$tissue=="THCA",]),2])
colnames(all_clin)[38] = "ME1"
ggplot(all_clin[all_clin$tissue == "THCA",], aes(x=radiation, y=VSTcnts)) + geom_boxplot()
ggsave("box.pdf")
m = lm(VSTcnts ~ radiation, data=all_clin[all_clin$tissue == "THCA",])
summary(m)
m = lm(ME1 ~ radiation, data=all_clin[all_clin$tissue == "THCA",])
summary(m)
ggplot(all_clin[all_clin$tissue == "THCA",], aes(x=radiation, y=ME1)) + geom_boxplot()
ggsave("boxME1.pdf")
normcnt <- read.table("../result.rsem.TET/normalizedcnt.txt", header=TRUE, sep="\t")
TEnames <- read.table("../result.rsem.TET/WGCNA/consensus/turquoise.1.TE.names")
turquoisecnts <- normcnt[as.character(TEnames[,1]),]
turquoiseSums <- colSums(turquoisecnts)
all_clin_THCA = cbind.data.frame(all_clin[all_clin$tissue == "THCA",],log2(turquoiseSums[all_clin$tissue == "THCA"]))
colnames(all_clin_THCA)[ncol(all_clin_THCA)] = "turquoiseSums"
ggplot(all_clin_THCA, aes(x=radiation, y=turquoiseSums)) + geom_boxplot()
ggsave("boxturquoise.pdf")

# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                 as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

# create vector for death censoring
table(all_clin$vital_status)
#   -- alive  dead 
#   11   433   253 

all_clin$death_event <- ifelse(all_clin$vital_status == 'alive', 0,1)

#finally add row.names to clinical
rownames(all_clin) <- all_clin$patient

# run survival

s <- survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~all_clin$L1HS5prime)
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~all_clin$L1HS5prime), error = function(e) return(NA))


all_clin$age_at_diagnosis = as.numeric(as.character(ll_clin$age_at_diagnosis))/365
all_clin$age_range = all_clin$age_at_diagnosis
all_clin$age_range[!is.na(all_clin$age_at_diagnosis) & all_clin$age_at_diagnosis < 40] = 2
all_clin$age_range[!is.na(all_clin$age_at_diagnosis) & all_clin$age_at_diagnosis >= 40 & all_clin$age_at_diagnosis < 60] = 3
all_clin$age_range[!is.na(all_clin$age_at_diagnosis) & all_clin$age_at_diagnosis >= 60] = 4
#all_clin$age_range[!is.na(all_clin$age_at_diagnosis) & all_clin$age_at_diagnosis >= 60 & all_clin$age_at_diagnosis < 80] = 4
#all_clin$age_range[!is.na(all_clin$age_at_diagnosis) & all_clin$age_at_diagnosis >= 80 ] =  5
#all_clin$age_range = factor(all_clin$age_range, levels=2:5, ordered=TRUE)

all_clin$L1HSactive = all_clin$VSTcnts
all_clin$L1HSactive[all_clin$VSTcnts >= 5] = "high"
all_clin$L1HSactive[all_clin$VSTcnts >= 3 & all_clin$VSTcnts < 5] = "iNA"
all_clin$L1HSactive[all_clin$VSTcnts < 3] = "low"



tissuetypes = levels(all_clin$tissue)
for (t in tissuetypes) {
# Build plot
ggplot(all_clin[all_clin$tissue==t,], aes(age_range, fill = L1HSactive)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent)
ggsave(paste0("tissue", t, ".pdf"))
}

for (t in tissuetypes) {
  ggplot(all_clin[all_clin$tissue==t,], aes(x=age_at_diagnosis, y=ME1)) + geom_point()
  ggsave(paste0("ME12agetissue", t, ".pdf"))
  m <- lm (ME1 ~ age_at_diagnosis, data= all_clin[all_clin$tissue==t,])
  #m <- lm (ME1 ~ age_range, data= all_clin[all_clin$tissue==t,])
  print(summary(m))
}

library(MASS)
m <- polr(factor(L1HSactive) ~ age_range , data = all_clin, Hess=TRUE)
summary(m)
