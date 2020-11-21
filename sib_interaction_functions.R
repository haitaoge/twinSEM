library(OpenMx)
library(dplyr)

source("miFunctions.R")
mxOption(NULL, "Default optimizer", "SLSQP")

#-------MZ/DZ only----------

twoACEsib <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # select data for analysis
  dat.mz <- select(mzdata, all_of(varnames))
  dat.dz <- select(dzdata, all_of(varnames))
  
  # statistics for starting values
  meanlab <- colMeans(dat.mz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svpA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svpC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  ifelse(is.nan(svpC) == TRUE, svpC <- 0, svpC <- svpC)
  svpE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # create mean matrices
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "mean", name = "meanP")
  
  # create path coef matrices
  covA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpA, labels = "VA", name = "A")
  covC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpC, labels = "VC", name = "C")
  covE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpE, labels = "VE", name = "E")
  
  # create sib interaction matrices
  matI <- mxMatrix(type = "Iden", nrow = 2, ncol = 2, name = "I")
  matB <- mxMatrix(type = "Symm", free = c(F,T,T,F), labels = c(NA,'b','b',NA),
                   nrow = 2, ncol = 2, name = "B")
  
  # create cov matrices
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")
  
  expCovMZ <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(A+C+E, A+C), 
                                                          cbind(A+C, A+C+E)),
                        name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(A+C+E, 0.5%x%A+C), 
                                                          cbind(0.5%x%A+C, A+C+E)),
                        name = "expCovDZ")
  
  # create data objects for multiple groups
  dataMZ <- mxData(observed = dat.mz, type = "raw")
  dataDZ <- mxData(observed = dat.dz, type = "raw")
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # model objects for multiple groups
  pars <- list(meanP, matI, matB, covA, covC, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))
  
  # algebra for variance components
  rowVC	<- rep('VC', nv)
  colVC <- rep(c('A','C','E','SA','SC','SE'), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, C, E, A/V, C/V, E/V), 
                     name = "VC", dimnames = list(rowVC, colVC))
  
  # ci's
  ciACE <- mxCI("VC[1,1:3]")
  
  # build model 
  modelACE <- mxModel("sibACE", pars, modelMZ, modelDZ, multi, estVC, ciACE)
  
  # run
  fitACE <- mxTryHard(modelACE, intervals = TRUE)
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor)
  returns <- list(model = modelACE, fit = fitACE, cors = correlations)
  return(returns)
}

#------three groups-------------

threeACEsib <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  adoptdata <- subset(data, zyg1 == 0)
  
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # select data for analysis
  dat.mz <- select(mzdata, all_of(varnames))
  dat.dz <- select(dzdata, all_of(varnames))
  dat.adopt <- select(adoptdata, all_of(varnames))
  
  # statistics for starting values
  meanlab <- colMeans(dat.mz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  adoptcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adoptcor <- cor(dat.adopt, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svpA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svpC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  ifelse(is.nan(svpC) == TRUE, svpC <- 0, svpC <- svpC)
  svpE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # create mean matrices
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "mean", name = "meanP")
  
  # create path coef matrices
  covA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpA, labels = "VA", name = "A")
  covC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpC, labels = "VC", name = "C", ubound = 1.195)
  covE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpE, labels = "VE", name = "E")
  
  # create sib interaction matrices
  matI <- mxMatrix(type = "Iden", nrow = 2, ncol = 2, name = "I")
  matB <- mxMatrix(type = "Symm", free = c(F,T,T,F), labels = c(NA,'b','b',NA),
                   nrow = 2, ncol = 2, name = "B")
  
  # create cov matrices
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")
  covAdopt <- mxAlgebra(expression = C, name = "cAdopt")
  
  expCovMZ <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(A+C+E, A+C), 
                                                          cbind(A+C, A+C+E)),
                        name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(A+C+E, 0.5%x%A+C), 
                                                          cbind(0.5%x%A+C, A+C+E)),
                        name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(A+C+E, C), 
                                                             cbind(C, A+C+E)),
                           name = "expCovAdopt")
  
  # create data objects for multiple groups
  dataMZ <- mxData(observed = dat.mz, type = "raw")
  dataDZ <- mxData(observed = dat.dz, type = "raw")
  dataAdopt <- mxData(observed = dat.adopt, type = "raw")
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt", means = "meanP", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # model objects for multiple groups
  pars <- list(meanP, matI, matB, covA, covC, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, dataAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ", "Adopt"))
  
  # algebra for variance components
  rowVC	<- rep('VC', nv)
  colVC <- rep(c('A','C','E','SA','SC','SE'), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, C, E, A/V, C/V, E/V), 
                     name = "VC", dimnames = list(rowVC, colVC))
  
  # ci's
  ciACE <- mxCI("VC[1,1:3]")
  
  # build model 
  modelACE <- mxModel("oneACEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciACE)
  
  # run
  fitACE <- mxTryHard(modelACE, intervals = TRUE)
  
  # returns
  covs <- list(mz = mzcov, dz = dzcov, adopt = adoptcov)
  correlations <- list(mz = mzcor, dz = dzcor, adopt = adoptcor)
  ACEvc <- list(modelACE, fitACE, correlations, covs)
  names(ACEvc) <- c("modelACE", "fitACE", "cors","covs")
  return(ACEvc)
}

threeACDEsib <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  adoptdata <- subset(data, zyg1 == 0)
  
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # select data for analysis
  dat.mz <- select(mzdata, all_of(varnames))
  dat.dz <- select(dzdata, all_of(varnames))
  dat.adopt <- select(adoptdata, all_of(varnames))
  
  # statistics for starting values
  meanlab <- colMeans(dat.mz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  adoptcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adoptcor <- cor(dat.adopt, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svpA <- .4
  svpC <- .2
  ifelse(is.nan(svpC) == TRUE, svpC <- 0, svpC <- svpC)
  svpD <- .2
  svpE <- .4
  
  # create mean matrices
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "mean", name = "meanP")
  
  # create path coef matrices
  covA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpA, labels = "VA", name = "A")
  covC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpC, labels = "VC", name = "C")
  covD <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE,
                   values = svpD, labels = "VD", name = "D")
  covE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                   values = svpE, lbound = 0.0001, labels = "VE", name = "E")
  
  # create sib interaction matrices
  matI <- mxMatrix(type = "Iden", nrow = 2, ncol = 2, name = "I")
  matB <- mxMatrix(type = "Symm", free = c(F,T,T,F), labels = c(NA,'b','b',NA),
                   nrow = 2, ncol = 2, name = "B")
  
  # create cov matrices
  covP <- mxAlgebra(expression = A + C + D + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C + D, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C + 0.25%x%D, name = "cDZ")
  covAdopt <- mxAlgebra(expression = C, name = "cAdopt")
  
  expCovMZ <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(V, cMZ), 
                                                          cbind(cMZ, V)),
                        name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(V, cDZ), 
                                                          cbind(cDZ, V)),
                        name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = solve(I-B) %&% rbind(cbind(V, cAdopt), 
                                                             cbind(cAdopt, V)),
                           name = "expCovAdopt")
  
  # create data objects for multiple groups
  dataMZ <- mxData(observed = dat.mz, type = "raw")
  dataDZ <- mxData(observed = dat.dz, type = "raw")
  dataAdopt <- mxData(observed = dat.adopt, type = "raw")
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt", means = "meanP", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # model objects for multiple groups
  pars <- list(meanP, matI, matB, covA, covC, covD, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, dataAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ", "Adopt"))
  
  # algebra for variance components
  rowVC	<- rep('cVC', nv)
  colVC <- rep(c('A','C','D', 'E','SA','SC','SD', 'SE'), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, C, D, E, A/V, C/V, D/V, E/V), 
                     name = "cVC", dimnames = list(rowVC, colVC))
  
  # ci's
  ciACE <- mxCI("cVC[1,1:4]")
  
  # build model 
  modelACDE <- mxModel("oneACDEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciACE)
  
  # run
  fitACDE <- mxTryHard(modelACE, intervals = TRUE)
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor, adopt = adoptcor)
  cov <- list(mz = mzcov, dz = dzcov, adopt = adoptcov)
  returns <- list(model = modelACDE, fit = fitACDE, cov = cov, cors = correlations)
  return(returns)
}

summit <- function(models, names, n) {
  newl <- c()
  for(i in models) {
    sum <- summary(i)
    cVC <- i$cVC
    vcov <- cov2cor(vcov(i))
    new <- list(sum = sum, cVC = cVC, vcov = vcov)
    newl <- append(newl, list(new))
  }
  names(newl) <- names
  if(n == 4) {AIC <- AICtab(models[[1]], models[[2]], models[[3]], models[[4]], 
                            weights = T, base = T, logLik = T)}
  if(n == 2) {AIC <- AICtab(models[[1]], models[[2]], weights = T, base = T, logLik = T)}
  compare <- mxCompare(models)
  
  returns <- list(s = newl, AIC = AIC, comp = compare)
  return(returns)
}