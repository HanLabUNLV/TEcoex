#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)



lm_module2trait <- function(eigenGenes, traits, resultdir, responsevar) {
  print(responsevar)
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
    module2trait <- module2trait[! is.na(module2trait[[responsevar]]),]    
    if (nrow(module2trait)<8) next;
    
    
    AICc = NULL    
    set.seed(0)
    garbage <- rnorm(nrow(module2trait))  
    if (cancertype == "THCA") {
      if (length(summary(module2trait$nested.batch))>1) {
        mod0 <- lm(as.formula(paste(responsevar, "~ garbage")), data=module2trait)
        mod1 <- lm(as.formula(paste(responsevar, "~ log2TPMsum")), data=module2trait)
        mod2 <- lm(as.formula(paste(responsevar, "~ modulePC1")), data=module2trait)
        mod3 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum")), data=module2trait)
        mod4 <- lm(as.formula(paste(responsevar, "~ nested.batch")), data=module2trait)
        mod5 <- lm(as.formula(paste(responsevar, "~ log2TPMsum + nested.batch")), data=module2trait)
        mod6 <- lm(as.formula(paste(responsevar, "~ modulePC1 + nested.batch")), data=module2trait)
        mod7 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum + nested.batch")), data=module2trait)
        mod8 <- lm(as.formula(paste(responsevar, "~ radiation")), data=module2trait)
        mod9 <- lm(as.formula(paste(responsevar, "~ log2TPMsum + radiation")), data=module2trait)
        mod10 <- lm(as.formula(paste(responsevar, "~ modulePC1 + radiation")), data=module2trait)
        mod11 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum + radiation")), data=module2trait)
        mod12 <- lm(as.formula(paste(responsevar, "~ nested.batch + radiation")), data=module2trait)
        mod13 <- lm(as.formula(paste(responsevar, "~ log2TPMsum + nested.batch + radiation")), data=module2trait)
        mod14 <- lm(as.formula(paste(responsevar, "~ modulePC1 + nested.batch + radiation")), data=module2trait)
        mod15 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum + nested.batch + radiation")), data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
      } else {
        mod0 <- lm(as.formula(paste(responsevar, "~ garbage")), data=module2trait)
        mod1 <- lm(as.formula(paste(responsevar, "~ log2TPMsum")), data=module2trait)
        mod2 <- lm(as.formula(paste(responsevar, "~ modulePC1")), data=module2trait)
        mod3 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum")), data=module2trait)
        mod4 <- lm(as.formula(paste(responsevar, "~ radiation")), data=module2trait)
        mod5 <- lm(as.formula(paste(responsevar, "~ log2TPMsum + radiation")), data=module2trait)
        mod6 <- lm(as.formula(paste(responsevar, "~ modulePC1 + radiation")), data=module2trait)
        mod7 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum + radiation")), data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      }      
    } else {
    if (length(summary(module2trait$nested.batch))>1) {
      mod0 <- lm(as.formula(paste(responsevar, "~ garbage")), data=module2trait)
      mod1 <- lm(as.formula(paste(responsevar, "~ log2TPMsum")), data=module2trait)
      mod2 <- lm(as.formula(paste(responsevar, "~ modulePC1")), data=module2trait)
      mod3 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum")), data=module2trait)
      mod4 <- lm(as.formula(paste(responsevar, "~ nested.batch")), data=module2trait)
      mod5 <- lm(as.formula(paste(responsevar, "~ log2TPMsum + nested.batch")), data=module2trait)
      mod6 <- lm(as.formula(paste(responsevar, "~ modulePC1 + nested.batch")), data=module2trait)
      mod7 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum + nested.batch")), data=module2trait)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
    } else {
      mod0 <- lm(as.formula(paste(responsevar, "~ garbage")), data=module2trait)
      mod1 <- lm(as.formula(paste(responsevar, "~ log2TPMsum")), data=module2trait)
      mod2 <- lm(as.formula(paste(responsevar, "~ modulePC1")), data=module2trait)
      mod3 <- lm(as.formula(paste(responsevar, "~ modulePC1 + log2TPMsum")), data=module2trait)
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
      if ( !any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
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
  filename = paste(resultdir,  paste0(responsevar, "_lm_", cancertypename, ".txt"), sep="/")
  results <- cbind.data.frame(modules[testedidx], fixed_moduleexp[testedidx], bestmodelnames[testedidx], moduleexp_coef[testedidx], moduleexp_pval[testedidx], moduleexp_qval, partialeta2[testedidx])
  colnames(results) <- c("modules", "fixed_moduleexp", "bestmodelnames", "moduleexp_coef", "moduleexp_pval", "moduleexp_qval", "partialeta2")
  results <- results[order(-results$partialeta2),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
  filename = paste(resultdir,  paste0(responsevar, "_positive_lm_", cancertypename, ".txt"), sep="/")
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




lm_trait2module <- function(eigenGenes, traits, resultdir, responsevar, as_factor) {
  print(responsevar)
  #eigenGenes  dim = (samples, modules)
  modules = colnames(eigenGenes)  
  cancertype = levels(factor(traits$newtype))
  cancertypename = traits$newtype[1]
  bestmodelnames <- rep("mod", length(modules))
  fixed_var <- rep(0, length(modules))
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
    module2trait <- module2trait[! is.na(module2trait[[responsevar]]),]    
    if (nrow(module2trait)<8) next;
    if (as_factor) {
      module2trait[responsevar] = as.factor(module2trait[[responsevar]])
      if (length(levels(module2trait[[responsevar]]))<2) next;
    } 
    
    
    AICc = NULL    
    set.seed(0)
    garbage <- rnorm(nrow(module2trait))  
    if (cancertype == "THCA") {
      if (length(summary(module2trait$nested.batch))>1) {
        mod0 <- lm(as.formula(paste("modulePC1 ~  garbage")), data=module2trait)
        mod1 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum")), data=module2trait)
        mod2 <- lm(as.formula(paste("modulePC1 ~ ", responsevar)), data=module2trait)
        mod3 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum")), data=module2trait)
        mod4 <- lm(as.formula(paste("modulePC1 ~  nested.batch")), data=module2trait)
        mod5 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum + nested.batch")), data=module2trait)
        mod6 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + nested.batch")), data=module2trait)
        mod7 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum + nested.batch")), data=module2trait)
        mod8 <- lm(as.formula(paste("modulePC1 ~  radiation")), data=module2trait)
        mod9 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum + radiation")), data=module2trait)
        mod10 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + radiation")), data=module2trait)
        mod11 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum + radiation")), data=module2trait)
        mod12 <- lm(as.formula(paste("modulePC1 ~  nested.batch + radiation")), data=module2trait)
        mod13 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum + nested.batch + radiation")), data=module2trait)
        mod14 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + nested.batch + radiation")), data=module2trait)
        mod15 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum + nested.batch + radiation")), data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
      } else {
        mod0 <- lm(as.formula(paste("modulePC1 ~  garbage")), data=module2trait)
        mod1 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum")), data=module2trait)
        mod2 <- lm(as.formula(paste("modulePC1 ~ ", responsevar)), data=module2trait)
        mod3 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum")), data=module2trait)
        mod4 <- lm(as.formula(paste("modulePC1 ~  radiation")), data=module2trait)
        mod5 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum + radiation")), data=module2trait)
        mod6 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + radiation")), data=module2trait)
        mod7 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum + radiation")), data=module2trait)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      }      
    } else {
    if (length(summary(module2trait$nested.batch))>1) {
      mod0 <- lm(as.formula(paste("modulePC1 ~  garbage")), data=module2trait)
      mod1 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum")), data=module2trait)
      mod2 <- lm(as.formula(paste("modulePC1 ~ ", responsevar)), data=module2trait)
      mod3 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum")), data=module2trait)
      mod4 <- lm(as.formula(paste("modulePC1 ~  nested.batch")), data=module2trait)
      mod5 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum + nested.batch")), data=module2trait)
      mod6 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + nested.batch")), data=module2trait)
      mod7 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum + nested.batch")), data=module2trait)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
    } else {
      mod0 <- lm(as.formula(paste("modulePC1 ~  garbage")), data=module2trait)
      mod1 <- lm(as.formula(paste("modulePC1 ~  log2TPMsum")), data=module2trait)
      mod2 <- lm(as.formula(paste("modulePC1 ~ ", responsevar)), data=module2trait)
      mod3 <- lm(as.formula(paste("modulePC1 ~ ", responsevar, " + log2TPMsum")), data=module2trait)
      AICc<-model.sel(mod0,mod1,mod2,mod3)
    }
    }
    
    bestmodel <- eval(getCall(AICc, 1))
    bestmodelnames[m] <- rownames(AICc)[1]
    sm <- summary(bestmodel)
    fixed_var[m] <- any(rownames(sm$coefficients)==responsevar)
    if (fixed_var[m]) {
      moduleexp_coef[m] <- sm$coefficients[responsevar,"Estimate"]
      moduleexp_pval[m] <- sm$coefficients[responsevar,"Pr(>|t|)"]
      if ( !any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2[m] <- effectsize$Effects[responsevar,"pEta-sqr"]
      } else {
        partialeta2[m] <- NA
      }
    }
    tested[m] <- 1
    
  }

  testedidx = which(tested==1)
  if (length(testedidx) > 0) {
  
  moduleexp_qval = qvalue(moduleexp_pval[testedidx],lambda=0)$qvalues 
  filename = paste(resultdir,  paste0(responsevar, "_lm_", cancertypename, ".txt"), sep="/")
  results <- cbind.data.frame(modules[testedidx], fixed_var[testedidx], bestmodelnames[testedidx], moduleexp_coef[testedidx], moduleexp_pval[testedidx], moduleexp_qval, partialeta2[testedidx])
  colnames(results) <- c("modules", "fixed_var", "bestmodelnames", "moduleexp_coef", "moduleexp_pval", "moduleexp_qval", "partialeta2")
  results <- results[order(-results$partialeta2),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
  filename = paste(resultdir,  paste0(responsevar, "_positive_lm_", cancertypename, ".txt"), sep="/")
  positive<-as.data.frame(results[results$fixed_var==1,])
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


