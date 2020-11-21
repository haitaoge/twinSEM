library(OpenMx)
library(dplyr)

source("miFunctions.R")
mxOption(NULL, "Default optimizer", "SLSQP")

# in this script: two rater rater bias models
# rater code: "mf" is mom, dad | "kf" is kid, dad | "mk" is mom, kid

twoACErbPath <- function(data, vars, raters) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  
  # variables for analysis
  nv        <- 2                         # number of variables
  ntv       <- nv*2                      # number of total variables
  selVars   <- paste(vars, c(rep("1",nv), rep("2",nv)), sep="")
  
  # data for analysis
  dat.mz <- select(mzdata, all_of(selVars)) 
  dat.dz <- select(dzdata, all_of(selVars))
  
  # statistics
  meanlab <- colMeans(dat.dz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  svE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # create path coef matrices
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = labVars("mean", vars), name = "meanP")
  
  coefA <- mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = TRUE, 
                    values = svA, labels = "a11", name = "a")
  coefC <- mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = TRUE, 
                    values = svC, labels = "c11", name = "c")
  coefE <- mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = TRUE, 
                    values = svE, labels = "e11", name = "e")
  # create cov matrices
  covA <- mxAlgebra(expression = a %*% t(a), name = "A")
  covC <- mxAlgebra(expression = c %*% t(c), name = "C")
  covE <- mxAlgebra(expression = e %*% t(e), name = "E")
  
  # rater specific matrices
  momB <- mxMatrix(type = "Full", nrow = 2, ncol = 2, 
                   free = c(T,F,T,F), values = .2, 
                   labels = c('momB',NA,'momB',NA), name = "mB")
  dadB <- mxMatrix(type = "Full", nrow = 2, ncol = 2, 
                   free = c(F,T,F,T), values = .2, 
                   labels = c(NA,'dadB',NA,'dadB'), name = "fB")
  kid1B <- mxMatrix(type = "Full", nrow = 2, ncol = 2,
                    free = c(T,F,T,F), values = .2,
                    labels = c('kid1B',NA,'kid1B',NA), name = "k1B")
  kid2B <- mxMatrix(type = "Full", nrow = 2, ncol = 2,
                    free = c(F,T,F,T), values = .2,
                    labels = c(NA,'kid2B',NA,'kid2B'), name = "k2B")
  
  if(raters == c("mf")) {allB <- mxAlgebra(expression = rbind(mB, fB), name = "aB")}
  if(raters == c("mk")) {allB <- mxAlgebra(expression = rbind(mB, k2B), name = "aB")}
  if(raters == c("kf")) {allB <- mxAlgebra(expression = rbind(k1B, fB), name = "aB")}
  
  
  momF <- mxMatrix(type = "Diag", nrow = 2, ncol = 2, 
                   free = TRUE, values = .2, 
                   labels = "momf", name = "mF")
  dadF <- mxMatrix(type = "Diag", nrow = 2, ncol = 2, 
                   free = TRUE, values = .2, 
                   labels = "dadf", name = "fF")
  kidF <- mxMatrix(type = "Diag", nrow = 2, ncol = 2,
                   free = TRUE, values = .2,
                   labels = "kidf", name = "kF")
  space <- mxMatrix(type = "Zero", nrow = 2, ncol = 2, name = "z")
  # residuals
  if(raters == c("mf")) {allF <- mxAlgebra(expression = rbind(cbind(mF,z), cbind(z,fF)), name = "aF")}
  if(raters == c("mk")) {allF <- mxAlgebra(expression = rbind(cbind(mF,z), cbind(z,kF)), name = "aF")}
  if(raters == c("kf")) {allF <- mxAlgebra(expression = rbind(cbind(kF,z), cbind(z,fF)), name = "aF")}
  
  
  matI <- mxMatrix(type = "Iden", nrow = 2, ncol = 2, name = "I")
  matS <- mxMatrix(type = "Diag", nrow = 2, ncol = 2, free = F,
                   values = 1, labels = "s", name = "S")
  matIS <- mxAlgebra(expression = rbind(I, S), name = "IS")
  
  # covariances
  covR <- mxAlgebra(expression = aB %*% t(aB), name = "R")
  covJ <- mxAlgebra(expression = aF %*% t(aF), name = "J")
  
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")
  
  expCovMZ <- mxAlgebra(expression = R + J + IS%&%rbind(cbind(V, cMZ), 
                                                        cbind(cMZ, V)), 
                        name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = R + J + IS%&%rbind(cbind(V, cDZ), 
                                                        cbind(cDZ, V)), 
                        name = "expCovDZ")
  
  # create data objects for multiple groups
  datMZ <- mxData(observed = dat.mz, type = "raw")
  datDZ <- mxData(observed = dat.dz, type = "raw")
  funML <- mxFitFunctionML()
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = selVars)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = selVars)
  
  # create model objects for multiple groups
  if(raters == c("mf")) {parsR <- c(momB, dadB, momF, dadF)}
  if(raters == c("mk")) {parsR <- c(momB, kid2B, momF, kidF)}
  if(raters == c("kf")) {parsR <- c(kid1B, dadB, kidF, dadF)}
  
  pars <- list(meanP, coefA, coefC, coefE, covA, covC, covE, covP,
               parsR, allB, space, allF, matI, matS, matIS, covR, covJ)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, datMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, datDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))
  
  # create algebra for variance components
  rowUV <- rep('UV', 1)
  colUV <- rep(c('VA','VC','VE'),each = 1)
  estUV <- mxAlgebra(expression = cbind(A,C,E), 
                         name = "UV", dimnames = list(rowUV,colUV))
  rowSV <- rep('SV', 1)
  colSV <- rep(c('SA','SC','SE'),each = 1)
  estSV <- mxAlgebra(expression = cbind(A/V,C/V,E/V), 
                         name = "SV", dimnames = list(rowSV,colSV))
  # CIs
  ciUV <- mxCI(c("UV"))
  ciSV <- mxCI(c("SV"))
  
  #build model with confidence intervals
  modelACErb <- mxModel("raterbiasACEc", pars, modelMZ, modelDZ, multi, estUV, estSV, ciUV, ciSV)
  fitACErb <- mxTryHard(modelACE, intervals=TRUE) #intervals = estimate CIs
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor)
  summary <- summary(fitACE)
  returns <- list(model = modelACErb, fit = fitACErb, s = summary, cors = correlations)
  return(returns)
}

twoACErbVC <- function(data, vars, raters) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  
  # variables for analysis
  nv        <- 2                         # number of variables
  ntv       <- nv*2                      # number of total variables
  selVars   <- paste(vars, c(rep("1",nv), rep("2",nv)), sep="")
  
  # data for analysis
  dat.mz <- select(mzdata, all_of(selVars)) 
  dat.dz <- select(dzdata, all_of(selVars))
  
  # statistics
  meanlab <- colMeans(dat.dz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  svE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # create path coef matrices
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = labVars("mean", vars), name = "meanP")
  
  covA <- mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = TRUE, 
                   values = svA, labels = "VA", name = "A") 
  covC <- mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = TRUE, 
                   values = svC, labels = "VC", name = "C")
  covE <- mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = TRUE, 
                   values = svE, labels = "VE", name = "E")
  
  # rater specific matrices
  momB <- mxMatrix(type = "Full", nrow = 2, ncol = 2, 
                   free = c(T,T,F,F), values = 0, 
                   labels = c('momB','momB',NA,NA), name = "mB")
  dadB <- mxMatrix(type = "Full", nrow = 2, ncol = 2, 
                   free = c(F,F,T,T), values = 0, 
                   labels = c(NA,NA,'dadB','dadB'), name = "fB")
  kid1B <- mxMatrix(type = "Full", nrow = 2, ncol = 2,
                    free = c(T,T,F,F), values = 0,
                    labels = c('kid1B','kid1B',NA,NA), name = "k1B")
  kid2B <- mxMatrix(type = "Full", nrow = 2, ncol = 2,
                    free = c(F,F,T,T), values = 0,
                    labels = c(NA,NA,'kid2B','kid2B'), name = "k2B")
  
  if(raters == c("mf")) {allB <- mxAlgebra(expression = rbind(mB, fB), name = "aB")}
  if(raters == c("mk")) {allB <- mxAlgebra(expression = rbind(mB, k2B), name = "aB")}
  if(raters == c("kf")) {allB <- mxAlgebra(expression = rbind(k1B, fB), name = "aB")}
  
  
  momF <- mxMatrix(type = "Diag", nrow = 2, ncol = 2, 
                   free = TRUE, values = .2, 
                   labels = "momf", name = "mF")
  dadF <- mxMatrix(type = "Diag", nrow = 2, ncol = 2, 
                   free = TRUE, values = .2, 
                   labels = "dadf", name = "fF")
  kidF <- mxMatrix(type = "Diag", nrow = 2, ncol = 2,
                   free = TRUE, values = .2,
                   labels = "kidf", name = "kF")
  space <- mxMatrix(type = "Zero", nrow = 2, ncol = 2, name = "z")
  allF <- mxAlgebra(expression = rbind(cbind(mF, z),
                                       cbind(z, fF)), name = "aF") # residuals
  # residuals
  if(raters == c("mf")) {allF <- mxAlgebra(expression = rbind(cbind(mF,z), cbind(z,fF)), name = "aF")}
  if(raters == c("mk")) {allF <- mxAlgebra(expression = rbind(cbind(mF,z), cbind(z,kF)), name = "aF")}
  if(raters == c("kf")) {allF <- mxAlgebra(expression = rbind(cbind(kF,z), cbind(z,fF)), name = "aF")}
  
  matI <- mxMatrix(type = "Iden", nrow = 2, ncol = 2, name = "I")
  matS <- mxMatrix(type = "Diag", nrow = 2, ncol = 2, free = F,
                   values = 1, labels = "s", name = "S")
  matIS <- mxAlgebra(expression = rbind(I, S), name = "IS")
  
  # covariances
  covR <- mxAlgebra(expression = aB %*% t(aB), name = "R")
  covJ <- mxAlgebra(expression = aF %*% t(aF), name = "J")
  
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")
  
  expCovMZ <- mxAlgebra(expression = R + J + IS%&%rbind(cbind(V, cMZ), 
                                                        cbind(cMZ, V)), 
                        name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = R + J + IS%&%rbind(cbind(V, cDZ), 
                                                        cbind(cDZ, V)), 
                        name = "expCovDZ")
  
  # create data objects for multiple groups
  datMZ <- mxData(observed = dat.mz, type = "raw")
  datDZ <- mxData(observed = dat.dz, type = "raw")
  funML <- mxFitFunctionML()
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = selVars)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = selVars)
  
  # create model objects for multiple groups
  if(raters == c("mf")) {parsR <- c(momB, dadB, momF, dadF)}
  if(raters == c("mk")) {parsR <- c(momB, kid2B, momF, kidF)}
  if(raters == c("kf")) {parsR <- c(kid1B, dadB, kidF, dadF)}
  
  pars <- list(meanP, covA, covC, covE, covP,
               parsR, allB, space, allF, matI, matS, matIS, covR, covJ)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, datMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, datDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))
  
  # create algebra for variance components
  rowUV <- rep('UV', 1)
  colUV <- rep(c('VA','VC','VE'),each = 1)
  estUV <- mxAlgebra(expression = cbind(A,C,E), 
                     name = "UV", dimnames = list(rowUV,colUV))
  rowSV <- rep('SV', 1)
  colSV <- rep(c('SA','SC','SE'),each = 1)
  estSV <- mxAlgebra(expression = cbind(A/V,C/V,E/V), 
                     name = "SV", dimnames = list(rowSV,colSV))
  # CIs
  ciUV <- mxCI(c("UV"))
  ciSV <- mxCI(c("SV"))
  
  #build model with confidence intervals
  modelACErb <- mxModel("raterbiasACEc", pars, modelMZ, modelDZ, multi, estUV, estSV, ciUV, ciSV)
  fitACErb <- mxTryHard(modelACE, intervals = TRUE) #intervals = estimate CIs
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor)
  summary <- summary(fitACE)
  returns <- list(model = modelACErb, fit = fitACErb, s = summary, cors = correlations)
  return(returns)
}