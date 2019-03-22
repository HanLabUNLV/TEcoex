#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)



lm_module2L1HS <- function(eigenGenes, traits, resultdir) {
  #eigenGenes  dim = (samples, modules)
  modules = colnames(eigenGenes)  
  cancertype = levels(factor(traits$newtype))
  cancertypename = traits$newtype[1]
  bestmodelnames <- rep("mod", length(modules))
  fixed_moduleexp <- rep(0, length(modules))
  moduleexp_coef <- rep(0,length(modules))
  moduleexp_pval <- rep(1,length(modules))
  partialeta2 <- rep(0, length(modules))
  tested <- rep(0, length(modules))
  
  for(m in 1:length(modules)){
    print (m)
    print (modules[m])
    gene<- eigenGenes[,m]
    gene = cbind.data.frame(rownames(eigenGenes), gene)
    colnames(gene) = c("patient", "modulePC1")

    module2trait <- merge(traits, gene, by = "patient")
    
    if (nrow(module2trait)<8) next;
    #if (all(colQuantiles(as.matrix(module2trait$modulePC1))["25%"]<log2(2))) next;
    module2trait$newtype <- factor(module2trait$newtype)
    module2trait$newtype <- droplevels(module2trait$newtype)
    module2trait$BatchId <- factor(module2trait$BatchId)
    if (length(cancertype)>1) {
      module2trait$BatchId <- droplevels(module2trait$BatchId)
      module2trait$nested.batch=as.numeric(module2trait$BatchId)
    }
    module2trait$nested.batch = factor(module2trait$nested.batch)
    
    
    AICc = NULL    
    set.seed(0)
    garbage <- rnorm(length(module2trait$L1HS))  
    if (cancertype == "THCA") {
      if (length(summary(module2trait$nested.batch))>1) {
        mod0 <- lm(L1HS ~ garbage, data=module2trait)
        mod1 <- lm(L1HS ~ log2TPMsum, data=module2trait)
        mod2 <- lm(L1HS ~ modulePC1, data=module2trait)
        mod3 <- lm(L1HS ~ modulePC1 + log2TPMsum, data=module2trait)
        mod4 <- lm(L1HS ~ nested.batch, data=module2trait)
        mod5 <- lm(L1HS ~ log2TPMsum + nested.batch, data=module2trait)
        mod6 <- lm(L1HS ~ modulePC1 + nested.batch, data=module2trait)
        mod7 <- lm(L1HS ~ modulePC1 + log2TPMsum + nested.batch, data=module2trait)
        mod8 <- lm(L1HS ~ radiation, data=module2trait)
        mod9 <- lm(L1HS ~ log2TPMsum + radiation, data=module2trait)
        mod10 <- lm(L1HS ~ modulePC1 + radiation, data=module2trait)
        mod11 <- lm(L1HS ~ modulePC1 + log2TPMsum + radiation, data=module2trait)
        mod12 <- lm(L1HS ~ nested.batch + radiation, data=module2trait)
        mod13 <- lm(L1HS ~ log2TPMsum + nested.batch + radiation, data=module2trait)
        mod14 <- lm(L1HS ~ modulePC1 + nested.batch + radiation, data=module2trait)
        mod15 <- lm(L1HS ~ modulePC1 + log2TPMsum + nested.batch + radiation, data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
      } else {
        mod0 <- lm(L1HS ~ garbage, data=module2trait)
        mod1 <- lm(L1HS ~ log2TPMsum, data=module2trait)
        mod2 <- lm(L1HS ~ modulePC1, data=module2trait)
        mod3 <- lm(L1HS ~ modulePC1 + log2TPMsum, data=module2trait)
        mod4 <- lm(L1HS ~ radiation, data=module2trait)
        mod5 <- lm(L1HS ~ log2TPMsum + radiation, data=module2trait)
        mod6 <- lm(L1HS ~ modulePC1 + radiation, data=module2trait)
        mod7 <- lm(L1HS ~ modulePC1 + log2TPMsum + radiation, data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      }      
    } else {
    if (length(summary(module2trait$nested.batch))>1) {
      mod0 <- lm(L1HS ~ garbage, data=module2trait)
      mod1 <- lm(L1HS ~ log2TPMsum, data=module2trait)
      mod2 <- lm(L1HS ~ modulePC1, data=module2trait)
      mod3 <- lm(L1HS ~ modulePC1 + log2TPMsum, data=module2trait)
      mod4 <- lm(L1HS ~ nested.batch, data=module2trait)
      mod5 <- lm(L1HS ~ log2TPMsum + nested.batch, data=module2trait)
      mod6 <- lm(L1HS ~ modulePC1 + nested.batch, data=module2trait)
      mod7 <- lm(L1HS ~ modulePC1 + log2TPMsum + nested.batch, data=module2trait)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
    } else {
      mod0 <- lm(L1HS ~ garbage, data=module2trait)
      mod1 <- lm(L1HS ~ log2TPMsum, data=module2trait)
      mod2 <- lm(L1HS ~ modulePC1, data=module2trait)
      mod3 <- lm(L1HS ~ modulePC1 + log2TPMsum, data=module2trait)
      AICc<-model.sel(mod0,mod1,mod2,mod3)
    }
    }
    
    bestmodel <- eval(getCall(AICc, 1))
    bestmodelnames[m] <- rownames(AICc)[1]
    sm <- summary(bestmodel)
    fixed_moduleexp[m] <- any(rownames(sm$coefficients)=="modulePC1")
    if (fixed_moduleexp[m]) {
      #plot(bestmodel$model[,"modulePC1"], resid(bestmodel), ylab="residuals", xlab="modulePC1")
      moduleexp_coef[m] <- sm$coefficients["modulePC1","Estimate"]
      moduleexp_pval[m] <- sm$coefficients["modulePC1","Pr(>|t|)"]
      if (modules[m] != "L1HS:L1:LINE" && !any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2[m] <- effectsize$Effects["modulePC1","pEta-sqr"]
      } else {
        partialeta2[m] <- NA
      }
    }
    tested[m] <- 1
    
  }

  testedidx = which(tested==1)
  if (length(testedidx) > 0) {
  
  moduleexp_qval = qvalue(moduleexp_pval[testedidx],lambda=0)$qvalues 
  filename = paste(resultdir,  paste0("L1HS_lm_", cancertypename, ".txt"), sep="/")
  results <- cbind.data.frame(modules[testedidx], fixed_moduleexp[testedidx], bestmodelnames[testedidx], moduleexp_coef[testedidx], moduleexp_pval[testedidx], moduleexp_qval, partialeta2[testedidx])
  colnames(results) <- c("modules", "fixed_moduleexp", "bestmodelnames", "moduleexp_coef", "moduleexp_pval", "moduleexp_qval", "partialeta2")
  results <- results[order(-results$partialeta2),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
  filename = paste(resultdir,  paste0("L1HS_positive_lm_", cancertypename, ".txt"), sep="/")
  positive<-as.data.frame(results[results$fixed_moduleexp==1,])
  positive[,4] <- as.numeric(as.character(positive[,4]))
  positive[,5] <- as.numeric(as.character(positive[,5]))
  positive[,6] <- as.numeric(as.character(positive[,6]))
  positive <- positive[order(-positive$partialeta2),]
  write.table(positive, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  }
  results <- cbind.data.frame(moduleexp_coef, moduleexp_pval, partialeta2)
  rownames(results) = modules
  return (results)
}



lm_module2oldLINE <- function(eigenGenes, traits, resultdir) {
  #eigenGenes  dim = (samples, modules)
  modules = colnames(eigenGenes)  
  cancertype = levels(factor(traits$newtype)) 
  cancertypename = traits$newtype[1]
  bestmodelnames <- rep("mod", length(modules))
  fixed_moduleexp <- rep(0, length(modules))
  moduleexp_coef <- rep(0,length(modules))
  moduleexp_pval <- rep(1,length(modules))
  partialeta2 <- rep(0, length(modules))
  tested <- rep(0, length(modules))
  
  for(m in 1:length(modules)){
    print (m)
    print (modules[m])
    gene<- eigenGenes[,m]
    gene = cbind.data.frame(rownames(eigenGenes), gene)
    colnames(gene) = c("patient", "modulePC1")

    module2trait <- merge(traits, gene, by = "patient")
    
    if (nrow(module2trait)<8) next;
    #if (all(colQuantiles(as.matrix(module2trait$modulePC1))["25%"]<log2(2))) next;
    module2trait$newtype <- factor(module2trait$newtype)
    module2trait$newtype <- droplevels(module2trait$newtype)
    module2trait$BatchId <- factor(module2trait$BatchId)
    if (length(cancertype)>1) {
      module2trait$BatchId <- droplevels(module2trait$BatchId)
      module2trait$nested.batch=as.numeric(module2trait$BatchId)
    }
    module2trait$nested.batch = factor(module2trait$nested.batch)
    
    
    AICc = NULL    
    set.seed(0)
    garbage <- rnorm(length(module2trait$oldLINE))  
    if (cancertype == "THCA") {
      if (length(summary(module2trait$nested.batch))>1) {
        mod0 <- lm(oldLINE ~ garbage, data=module2trait)
        mod1 <- lm(oldLINE ~ log2TPMsum, data=module2trait)
        mod2 <- lm(oldLINE ~ modulePC1, data=module2trait)
        mod3 <- lm(oldLINE ~ modulePC1 + log2TPMsum, data=module2trait)
        mod4 <- lm(oldLINE ~ nested.batch, data=module2trait)
        mod5 <- lm(oldLINE ~ log2TPMsum + nested.batch, data=module2trait)
        mod6 <- lm(oldLINE ~ modulePC1 + nested.batch, data=module2trait)
        mod7 <- lm(oldLINE ~ modulePC1 + log2TPMsum + nested.batch, data=module2trait)
        mod8 <- lm(oldLINE ~ radiation, data=module2trait)
        mod9 <- lm(oldLINE ~ log2TPMsum + radiation, data=module2trait)
        mod10 <- lm(oldLINE ~ modulePC1 + radiation, data=module2trait)
        mod11 <- lm(oldLINE ~ modulePC1 + log2TPMsum + radiation, data=module2trait)
        mod12 <- lm(oldLINE ~ nested.batch + radiation, data=module2trait)
        mod13 <- lm(oldLINE ~ log2TPMsum + nested.batch + radiation, data=module2trait)
        mod14 <- lm(oldLINE ~ modulePC1 + nested.batch + radiation, data=module2trait)
        mod15 <- lm(oldLINE ~ modulePC1 + log2TPMsum + nested.batch + radiation, data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
      } else {
        mod0 <- lm(oldLINE ~ garbage, data=module2trait)
        mod1 <- lm(oldLINE ~ log2TPMsum, data=module2trait)
        mod2 <- lm(oldLINE ~ modulePC1, data=module2trait)
        mod3 <- lm(oldLINE ~ modulePC1 + log2TPMsum, data=module2trait)
        mod4 <- lm(oldLINE ~ radiation, data=module2trait)
        mod5 <- lm(oldLINE ~ log2TPMsum + radiation, data=module2trait)
        mod6 <- lm(oldLINE ~ modulePC1 + radiation, data=module2trait)
        mod7 <- lm(oldLINE ~ modulePC1 + log2TPMsum + radiation, data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      }      
    } else {
    if (length(summary(module2trait$nested.batch))>1) {
      mod0 <- lm(oldLINE ~ garbage, data=module2trait)
      mod1 <- lm(oldLINE ~ log2TPMsum, data=module2trait)
      mod2 <- lm(oldLINE ~ modulePC1, data=module2trait)
      mod3 <- lm(oldLINE ~ modulePC1 + log2TPMsum, data=module2trait)
      mod4 <- lm(oldLINE ~ nested.batch, data=module2trait)
      mod5 <- lm(oldLINE ~ log2TPMsum + nested.batch, data=module2trait)
      mod6 <- lm(oldLINE ~ modulePC1 + nested.batch, data=module2trait)
      mod7 <- lm(oldLINE ~ modulePC1 + log2TPMsum + nested.batch, data=module2trait)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
    } else {
      mod0 <- lm(oldLINE ~ garbage, data=module2trait)
      mod1 <- lm(oldLINE ~ log2TPMsum, data=module2trait)
      mod2 <- lm(oldLINE ~ modulePC1, data=module2trait)
      mod3 <- lm(oldLINE ~ modulePC1 + log2TPMsum, data=module2trait)
      AICc<-model.sel(mod0,mod1,mod2,mod3)
    }
    }
    
    bestmodel <- eval(getCall(AICc, 1))
    bestmodelnames[m] <- rownames(AICc)[1]
    sm <- summary(bestmodel)
    fixed_moduleexp[m] <- any(rownames(sm$coefficients)=="modulePC1")
    if (fixed_moduleexp[m]) {
      #plot(bestmodel$model[,"modulePC1"], resid(bestmodel), ylab="residuals", xlab="modulePC1")
      moduleexp_coef[m] <- sm$coefficients["modulePC1","Estimate"]
      moduleexp_pval[m] <- sm$coefficients["modulePC1","Pr(>|t|)"]
      if (modules[m] != "L1HS:L1:LINE" && !any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2[m] <- effectsize$Effects["modulePC1","pEta-sqr"]
      } else {
        partialeta2[m] <- NA
      }
    }
    tested[m] <- 1
    
  }

  testedidx = which(tested==1)
  if (length(testedidx) > 0) {
  
  moduleexp_qval = qvalue(moduleexp_pval[testedidx],lambda=0)$qvalues 
  filename = paste(resultdir,  paste0("oldLINE_lm_", cancertypename, ".txt"), sep="/")
  results <- cbind.data.frame(modules[testedidx], fixed_moduleexp[testedidx], bestmodelnames[testedidx], moduleexp_coef[testedidx], moduleexp_pval[testedidx], moduleexp_qval, partialeta2[testedidx])
  colnames(results) <- c("modules", "fixed_moduleexp", "bestmodelnames", "moduleexp_coef", "moduleexp_pval", "moduleexp_qval", "partialeta2")
  results <- results[order(-results$partialeta2),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
  filename = paste(resultdir,  paste0("oldLINE_positive_lm_", cancertypename, ".txt"), sep="/")
  positive<-as.data.frame(results[results$fixed_moduleexp==1,])
  positive[,4] <- as.numeric(as.character(positive[,4]))
  positive[,5] <- as.numeric(as.character(positive[,5]))
  positive[,6] <- as.numeric(as.character(positive[,6]))
  positive <- positive[order(-positive$partialeta2),]
  write.table(positive, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  }
  results <- cbind.data.frame(moduleexp_coef, moduleexp_pval, partialeta2)
  rownames(results) = modules
  return (results)
}

