library(OpenMx)
library(dplyr)

# models in this function script
# saturated models two/three groups
# two group ACE/ADE VC and cholesky
# three group ACE/ADE/ACDE VC and cholesky 
#   (why is there ACE/ADE when you could just run submodels? good question)
# submodels for ACDE by cholesky and VC and umx just for good measure

#---------saturated model-------------
twoSAT <- function(vars, data) {
  nv <- 1
  ntv <- 2
  varnames <- paste(vars, sep = "",
                   c(rep(1, nv), rep(2, nv)))
  
  # data for analysis
  mzdata <- subset(data, zyg1 == 1, varnames)
  dzdata <- subset(data, zyg1 == 2, varnames)
  
  svMe <- .1                        # start value for means
  svVa <- .8                        # start value for variance
  lbVa <- .0001                     # lower bound for variance
  
  # Create Algebra for expected Mean Matrices
  meanMZ <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, 
                     labels = c("mMZ1","mMZ2"), name = "meanMZ")
  meanDZ <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, 
                     labels = c("mDZ1","mDZ2"), name = "meanDZ")
  
  # Create Algebra for expected Variance/Covariance Matrices
  covMZ <- mxMatrix(type = "Symm", nrow = ntv, ncol = ntv, free = TRUE, values = valDiag(svVa,ntv), 
                    lbound = valDiag(lbVa,ntv), labels = c("vMZ1","cMZ21","vMZ2"), name = "covMZ")
  covDZ <- mxMatrix(type = "Symm", nrow = ntv, ncol = ntv, free = TRUE, values = valDiag(svVa,ntv), 
                    lbound = valDiag(lbVa,ntv), labels = c("vDZ1","cDZ21","vDZ2"), name = "covDZ")
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData(observed = mzdata, type = "raw")
  dataDZ <- mxData(observed = dzdata, type = "raw")
  
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal(covariance = "covMZ", means = "meanMZ", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "covDZ", means = "meanDZ", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  modelMZ <- mxModel(meanMZ, covMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(meanDZ, covDZ, dataDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ","DZ"))
  
  # Create Confidence Interval Objects
  ciCov <- mxCI(c('MZ.covMZ','DZ.covDZ'))
  ciMean <- mxCI(c('MZ.meanMZ','DZ.meanDZ'))
  
  # Build Saturated Model with Confidence Intervals
  modelSAT <- mxModel("oneSATc", modelMZ, modelDZ, multi, ciCov, ciMean)
  fitSAT <- mxTryHard(modelSAT, intervals = F)
  expected <- mxGetExpected(fitSAT, c("means","covariance"))
  
  returns <- list(model = modelSAT, fit = fitSAT, exp = expected)
  return(returns)
}

submodelsSATtwo <- function(modelfit) {
  svMe <- .1                        # start value for means
  svVa <- .8                        # start value for variance
  
  # Constrain expected Means to be equal across Twin Order
  modelEMO <- mxModel(modelfit, name = "oneEMOc")
  modelEMO <- omxSetParameters(modelEMO, label = c("mMZ1","mMZ2"), free = TRUE, 
                               values = svMe, newlabels = 'mMZ')
  modelEMO <- omxSetParameters(modelEMO, label = c("mDZ1","mDZ2"), free = TRUE, 
                               values = svMe, newlabels = 'mDZ')
  fitEMO <- mxTryHard(modelEMO, intervals = F)
  
  # Constrain expected Means and Variances to be equal across Twin Order
  modelEMVO <- mxModel(fitEMO, name = "oneEMVOc")
  modelEMVO <- omxSetParameters(modelEMVO, label = c("vMZ1","vMZ2"), free = TRUE, 
                                values = svVa, newlabels = 'vMZ')
  modelEMVO <- omxSetParameters(modelEMVO, label = c("vDZ1","vDZ2"), free = TRUE, 
                                values = svVa, newlabels = 'vDZ')
  fitEMVO <- mxTryHard(modelEMVO, intervals = F )
  
  # Constrain expected Means and Variances to be equal across Twin Order and Zygosity
  modelEMVZ <- mxModel(fitEMVO, name = "oneEMVZc")
  modelEMVZ <- omxSetParameters(modelEMVZ, label = c("mMZ","mDZ"), free = TRUE, 
                                values = svMe, newlabels = 'mZ' )
  modelEMVZ <- omxSetParameters(modelEMVZ, label = c("vMZ","vDZ"), free = TRUE, 
                                values = svVa, newlabels = 'vZ' )
  fitEMVZ <- mxTryHard(modelEMVZ, intervals = F )
  
  # Print Comparative Fit Statistics
  compare <- mxCompare(modelfit, c(fitEMO, fitEMVO, fitEMVZ))
  
  AIC <- bbmle::AICtab(modelfit, fitEMVO, fitEMVZ, fitEMO, weights = T, base = T, logLik = T)
  
  returns <- list(modelEMO = modelEMO, fitEMO = fitEMO, modelEMVO = modelEMVO, 
                  fitEMVO = fitEMVO, modelEMVZ = modelEMVZ, fitEMVZ = fitEMVZ, 
                  compare = compare, AIC = AIC)
  return(returns)
}

threeSAT <- function(vars, data) {
  nv <- 1
  ntv <- 2
  varnames <- paste(vars, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # data for analysis
  mzdat <- subset(data, zyg1 == 1, varnames)
  dzdat <- subset(data, zyg1 == 2, varnames)
  adopdat <- subset(data, zyg1 == 0, varnames)
  
  svMe      <- .1                        # start value for means
  svVa      <- .8                        # start value for variance
  lbVa      <- .0001                     # lower bound for variance
  
  # Create Algebra for expected Mean Matrices
  meanMZ    <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                        values = svMe, labels = c("mMZ1","mMZ2"), name = "meanMZ" )
  meanDZ    <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                        values = svMe, labels = c("mDZ1","mDZ2"), name = "meanDZ" )
  meanAdopt <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                        values = svMe, labels = c("mAD1", "mAD2"), name = "meanAdopt")
  
  # Create Algebra for expected Variance/Covariance Matrices
  covMZ <- mxMatrix(type = "Symm", nrow = ntv, ncol = ntv, free = TRUE, values = valDiag(svVa,ntv), 
                    lbound = valDiag(lbVa,ntv), labels = c("vMZ1","cMZ21","vMZ2"), name = "covMZ")
  covDZ <- mxMatrix(type = "Symm", nrow = ntv, ncol = ntv, free = TRUE, values = valDiag(svVa,ntv), 
                    lbound = valDiag(lbVa,ntv), labels = c("vDZ1","cDZ21","vDZ2"), name = "covDZ")
  covAdopt <- mxMatrix(type = "Symm", nrow = ntv, ncol = ntv, free = TRUE, values = valDiag(svVa,ntv), 
                       lbound = valDiag(lbVa,ntv), labels = c("vAD1","cAD21","vAD2"), name = "covAdopt")
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData(observed = mzdat, type = "raw")
  dataDZ <- mxData(observed = dzdat, type = "raw")
  dataAdopt <- mxData(observed = adopdat, type = "raw")
  
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal(covariance = "covMZ", means = "meanMZ", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "covDZ", means = "meanDZ", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "covAdopt", means = "meanAdopt", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  modelMZ <- mxModel(meanMZ, covMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(meanDZ, covDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(meanAdopt, covAdopt, dataAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ","DZ", "Adopt"))
  
  # Create Confidence Interval Objects
  ciCov <- mxCI(c('MZ.covMZ','DZ.covDZ', 'Adopt.covAdopt'))
  ciMean <- mxCI(c('MZ.meanMZ','DZ.meanDZ', 'Adopt.meanAdopt'))
  
  # Build Saturated Model with Confidence Intervals
  modelSAT <- mxModel("oneSATc", modelMZ, modelDZ, modelAdopt, multi, ciCov, ciMean)
  fitSAT <- mxTryHard(modelSAT, intervals = F)
  expected <- mxGetExpected(fitSAT, c("means","covariance"))
  
  returns <- list(model = modelSAT, fit = fitSAT, exp = expected)
  return(returns)
}

submodelsSAT <- function(modelfit) {
  svMe <- .1                        # start value for means
  svVa <- .8                        # start value for variance
  
  # Constrain expected Means to be equal across Twin Order
  modelEMO <- mxModel(modelfit, name = "oneEMOc" )
  modelEMO <- omxSetParameters(modelEMO, label = c("mMZ1","mMZ2"), free = TRUE, 
                               values = svMe, newlabels = 'mMZ')
  modelEMO <- omxSetParameters(modelEMO, label = c("mDZ1","mDZ2"), free = TRUE, 
                               values = svMe, newlabels = 'mDZ')
  modelEMO <- omxSetParameters(modelEMO, label = c("mAD1","mAD2"), free = TRUE, 
                               values = svMe, newlabels = 'mAD')
  fitEMO <- mxTryHard(modelEMO, intervals = F)
  
  # Constrain expected Means and Variances to be equal across Twin Order
  modelEMVO <- mxModel(fitEMO, name = "oneEMVOc" )
  modelEMVO <- omxSetParameters(modelEMVO, label = c("vMZ1","vMZ2"), free = TRUE, 
                                values = svVa, newlabels = 'vMZ')
  modelEMVO <- omxSetParameters(modelEMVO, label = c("vDZ1","vDZ2"), free = TRUE, 
                                values = svVa, newlabels = 'vDZ')
  modelEMVO <- omxSetParameters(modelEMVO, label = c("vAD1","vAD2"), free = TRUE, 
                                values = svVa, newlabels = 'vAD')
  fitEMVO <- mxTryHard(modelEMVO, intervals = F )
  
  # Constrain expected Means and Variances to be equal across Twin Order and Zygosity
  modelEMVZ <- mxModel(fitEMVO, name = "oneEMVZc")
  modelEMVZ <- omxSetParameters(modelEMVZ, label = c("mMZ","mDZ",'mAD'), free = TRUE, 
                                values = svMe, newlabels = 'mZ')
  modelEMVZ <- omxSetParameters(modelEMVZ, label = c("vMZ","vDZ",'vAD'), free = TRUE, 
                                values = svVa, newlabels = 'vZ')
  fitEMVZ <- mxTryHard(modelEMVZ, intervals = F)
  
  # Print Comparative Fit Statistics
  compare <- mxCompare(modelfit, c(fitEMO, fitEMVO, fitEMVZ))
  AIC <- bbmle::AICtab(modelfit, fitEMVO, fitEMVZ, fitEMO, weights = T, base = T, logLik = T)
  
  returns <- list(modelEMO = modelEMO, fitEMO = fitEMO, modelEMVO = modelEMVO, 
                  fitEMVO = fitEMVO, modelEMVZ = modelEMVZ, fitEMVZ = fitEMVZ, 
                  compare = compare, AIC = AIC)
  return(returns)
}

#------two groups--------

twoACEpath <- function(data, var) {
  mzdat <- subset(data, zyg1 == 1)
  dzdat <- subset(data, zyg1 == 2)
  
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
  coefA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpA, labels = "a11", name = "a")
  coefC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpC, labels = "c11", name = "c")
  coefE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpE, labels = "e11", name = "e")
  
  # create cov matrices
  covA <- mxAlgebra(expression = a %*% t(a), name = "A")
  covC <- mxAlgebra(expression = c %*% t(c), name = "C")
  covE <- mxAlgebra(expression = e %*% t(e), name = "E")
  
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")
  
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(cMZ, V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(cDZ, V)), name = "expCovDZ")
  
  # create data objects for multiple groups
  dataMZ <- mxData(observed = dat.mz, type = "raw")
  dataDZ <- mxData(observed = dat.dz, type = "raw")
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # model objects for multiple groups
  pars <- list(meanP, coefA, coefC, coefE, covA, covC, covE, covP)
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
  modelACE <- mxModel("oneACEc", pars, modelMZ, modelDZ, multi, estVC, ciACE)
  
  # run
  fitACE <- mxRun(modelACE, intervals = TRUE)
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor)
  returns <- list(model = modelACE, fit = fitACE, cors = correlations)
  return(returns)
}

twoACEVC <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # data for analysis
  dat.mz <- select(mzdata, all_of(varnames)) 
  dat.dz <- select(dzdata, all_of(varnames))
  
  # statistics
  meanlab <- colMeans(dat.dz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  
  # set starting values
  # note: if you get NA errors it sets starting values to 0
  svmu <- meanlab[1]
  svA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  ifelse(is.nan(svC) == TRUE, svC <- 0, svC <- svC)
  svE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # mean and variance components
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "x1", name = "meanP")
  
  covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svA, labels = "VA", name = "A")
  covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svC, labels = "VC", name = "C")
  covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svE, labels = "VE", name = "E")
  
  # covariances
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")

  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name = "expCovDZ")
  
  # create data objects for multiple groups
  datMZ <- mxData(observed = dat.mz, type = "raw")
  datDZ <- mxData(observed = dat.dz, type = "raw")
  funML <- mxFitFunctionML()
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  
  # create model objects for multiple groups
  pars <- list(meanP, covA, covC, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, datMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, datDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))
  
  # create algebra for variance components
  rowVC <- rep("cVC", nv)
  colVC <- rep(c("A", "C", "E", "SA", "SC", "SE"), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, C, E, A/V, C/V, E/V),
                     name = "cVC", dimnames = list(rowVC, colVC))
  # CIs
  ciACE <- mxCI("cVC[1,1:3]")
  
  #build model with confidence intervals
  modelACE <- mxModel("oneACEc", pars, modelMZ, modelDZ, multi, estVC, ciACE)
  fitACE <- mxRun(modelACE, intervals=TRUE) #intervals = estimate CIs
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor)
  returns <- list(model = modelACE, fit = fitACE, cors = correlations)
  return(ACEvar)
}

#------path coef functions--------

threeACEpath <- function(data, var) {
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
  adopcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adopcor <- cor(dat.adopt, use = "complete")
  
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
  coefA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpA, labels = "a11", name = "a")
  coefC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpC, labels = "c11", name = "c")
  coefE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpE, labels = "e11", name = "e")
  
  # create cov matrices
  covA <- mxAlgebra(expression = a %*% t(a), name = "A")
  covC <- mxAlgebra(expression = c %*% t(c), name = "C")
  covE <- mxAlgebra(expression = e %*% t(e), name = "E")
  
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")
  covAdopt <- mxAlgebra(expression = C, name = "cAdopt")
  
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(cMZ, V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(cDZ, V)), name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(cAdopt, V)), 
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
  pars <- list(meanP, coefA, coefC, coefE, covA, covC, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, dataAdopt, expAdopt,
                        funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ","Adopt"))
  
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
  fitACE <- mxRun(modelACE, intervals = TRUE)
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor, adopt = adopcor)
  returns <- list(model = modelACE, fit = fitACE, cors = correlations)
  return(ACEpath)
}

threeADEpath <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  adoptdata <- subset(data, zyg1 == 0)
  
  # variables for analysis
  var <- variable
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # select data for analysis - assuming zyg is coded like this
  # change if zyg is coded differently
  dat.mz <- select(mzdata, all_of(varnames))
  dat.dz <- select(dzdata, all_of(varnames))
  dat.adopt <- select(adoptdata, all_of(varnames))
  
  # statistics for starting values
  meanlab <- colMeans(dat.mz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  adopcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adopcor <- cor(dat.adopt, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svpA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svpD <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  ifelse(is.nan(svpD) == TRUE, svpD <- 0, svpD <- svpD)
  svpE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # create mean matrices
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "mean", name = "meanP")
  
  # create path coef matrices
  coefA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpA, labels = "a11", name = "a")
  coefD <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpD, labels = "d11", name = "d")
  coefE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svpE, labels = "e11", name = "e")
  
  # create cov matrices
  covA <- mxAlgebra(expression = a %*% t(a), name = "A")
  covD <- mxAlgebra(expression = d %*% t(d), name = "D")
  covE <- mxAlgebra(expression = e %*% t(e), name = "E")
  
  covP <- mxAlgebra(expression = A + D + E, name = "V")
  covMZ <- mxAlgebra(expression = A + D, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + 0.25%x%D, name = "cDZ")
  covAdopt <- mxAlgebra(expression = 0, name = "cAdopt")
  
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(cMZ, V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(cDZ, V)), name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(cAdopt, V)),
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
  pars <- list(meanP, coefA, coefD, coefE, covA, covD, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, dataAdopt, expAdopt,
                        funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ","Adopt"))
  
  # algebra for variance components
  rowVC	<- rep('VC', nv)
  colVC <- rep(c('A','D','E','SA','SD','SE'), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, D, E, A/V, D/V, E/V), 
                     name = "VC", dimnames = list(rowVC, colVC))
  
  # ci's
  ciADE <- mxCI("VC[1,1:3]")
  
  # build model 
  modelADE <- mxModel("oneADEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciADE)
  
  # run
  fitADE <- mxRun(modelADE, intervals = TRUE)
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor, adopt = adopcor)
  returns <- list(model = modelADE, fit = fitADE, cors = correlations)
  return(returns)
}

threeACDEpath <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  adoptdata <- subset(data, zyg1 == 0)
  
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # select data for analysis - assuming zyg is coded like this
  # change if zyg is coded differently
  dat.mz <- select(mzdata, all_of(varnames))
  dat.dz <- select(dzdata, all_of(varnames))
  dat.adopt <- select(adoptdata, all_of(varnames))
  
  # statistics for starting values
  meanlab <- colMeans(dat.mz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  adopcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adopcor <- cor(dat.adopt, use = "complete")
  
  # set starting values
  # note: if you get NA errors it sets the starting values to 0
  # please do not be alarmed by the warnings
  svmu <- meanlab[1]
  svpA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svpC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  ifelse(is.nan(svpC) == TRUE, svpC <- 0, svpC <- svpC)
  svpD <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  ifelse(is.nan(svpD) == TRUE, svpD <- 0, svpD <- svpD)
  svpE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # create mean matrices
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv,
                    free = TRUE, values = svmu,
                    labels = "mean", name = "meanP")
  
  # create path coef matrices
  coefA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv,
                    free = TRUE, values = svpA,
                    labels = "a11", name = "a")
  coefC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv,
                    free = TRUE, values = svpC,
                    labels = "c11", name = "c")
  coefD <- mxMatrix(type = "Lower", nrow = nv, ncol = nv,
                    free = TRUE, values = svpD,
                    labels = "d11", name = "d")
  coefE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv,
                    free = TRUE, values = svpE,
                    labels = "e11", name = "e")
  
  # create cov matrices
  covA <- mxAlgebra(expression = a %*% t(a), name = "A")
  covC <- mxAlgebra(expression = c %*% t(c), name = "C")
  covD <- mxAlgebra(expression = d %*% t(d), name = "D")
  covE <- mxAlgebra(expression = e %*% t(e), name = "E")
  
  covP <- mxAlgebra(expression = A + C + D + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C + D, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C + 0.25%x%D, name = "cDZ")
  covAdopt <- mxAlgebra(expression = C, name = "cAdopt")
  
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(cMZ, V)),
                        name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(cDZ, V)),
                        name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(cAdopt, V)),
                           name = "expCovAdopt")
  
  # create data objects for multiple groups
  dataMZ <- mxData(observed = dat.mz, type = "raw")
  dataDZ <- mxData(observed = dat.dz, type = "raw")
  dataAdopt <- mxData(observed = dat.adopt, type = "raw")
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ",
                               means = "meanP",
                               dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ",
                               means = "meanP",
                               dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt",
                                  means = "meanP",
                                  dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # model objects for multiple groups
  pars <- list(meanP, coefA, coefC, coefD, coefE, covA, covC, covD, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, dataAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ","Adopt"))
  
  # algebra for variance components
  rowVC	<- rep('VC', nv)
  colVC <- rep(c('A','C','D','E','SA','SC','SD','SE'), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, C, D, E, A/V, C/V, D/V, E/V), 
                     name = "VC", dimnames = list(rowVC, colVC))
  
  # ci's
  ciACDE <- mxCI("VC[1,1:4]")
  
  # build model 
  modelACDE <- mxModel("oneACDEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciACDE)
  
  # run
  fitACDE <- mxRun(modelACDE, intervals = TRUE)
  
  # returns
  correlations <- rbind(mz = mzcor, dz = dzcor, adopt = adopcor)
  ACDEpath <- list(modelACDE, fitACDE, correlations)
  names(ACDEpath) <- c("model", "fit", "cors")
  return(ACDEpath)
}

#-----var component functions--------

threeACEVC <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  adoptdata <- subset(data, zyg1 == 0)
  
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # data for analysis
  dat.mz <- select(mzdata, all_of(varnames)) 
  dat.dz <- select(dzdata, all_of(varnames))
  dat.adopt <- select(adoptdata, all_of(varnames))
  
  # statistics
  meanlab <- colMeans(dat.dz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  adopcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adopcor <- cor(dat.adopt, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  svE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # mean and variance components
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "mean", name = "meanP")
  
  covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svA, labels = "VA", name = "A")
  covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svC, labels = "VC", name = "C")
  covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svE, labels = "VE", name = "E")
  
  # covariances
  covP <- mxAlgebra(expression = A + C + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C, name = "cDZ")
  covAdopt <- mxAlgebra(expression = C, name = "cAdopt")

  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(t(cAdopt), V)), name = "expCovAdopt")
  
  # create data objects for multiple groups
  datMZ <- mxData(observed = dat.mz, type = "raw")
  datDZ <- mxData(observed = dat.dz, type = "raw")
  datAdopt <- mxData(observed = dat.adopt, type = "raw")
  
  funML <- mxFitFunctionML()
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt", means = "meanP", dimnames = varnames)
  
  # create model objects for multiple groups
  pars <- list(meanP, covA, covC, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, datMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, datDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, datAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ", "Adopt"))
  
  # create algebra for variance components
  rowVC <- rep("cVC", nv)
  colVC <- rep(c("A", "C", "E", "SA", "SC", "SE"), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, C, E, A/V, C/V, E/V),
                      name = "cVC", dimnames = list(rowVC, colVC))
  # CIs
  ciACE <- mxCI("cVC[1,1:3]")
  
  #build model with confidence intervals
  modelACE <- mxModel("oneACEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciACE)
  fitACE <- mxRun(modelACE, intervals=TRUE) #intervals = estimate CIs
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor, adopt = adopcor)
  returns <- list(model = modelACE, fit = fitACE, cors = correlations)
  return(returns)
}

threeADEVC <- function(data, var) {
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  adoptdata <- subset(data, zyg1 == 0)
  
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  # data for analysis
  dat.mz <- select(mzdata, all_of(varnames)) 
  dat.dz <- select(dzdata, all_of(varnames))
  dat.adopt <- select(adoptdata, all_of(varnames))
  
  # statistics
  meanlab <- colMeans(dat.dz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  adopcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adopcor <- cor(dat.adopt, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svD <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  svE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # mean and variance components
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "x1", name = "meanP")
  
  covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svA, labels = "VA", name = "A")
  covD <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svD, labels = "VD", name = "D")
  covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svE, labels = "VE", name = "E")
  
  # covariances
  covP <- mxAlgebra(expression = A + D + E, name = "V")
  covMZ <- mxAlgebra(expression = A + D, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + 0.25%x%D, name = "cDZ")
  covAdopt <- mxAlgebra(expression = 0, name = "cAdopt")
  
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(t(cAdopt), V)), name = "expCovAdopt")
  
  # create data objects for multiple groups
  datMZ <- mxData(observed = dat.mz, type = "raw")
  datDZ <- mxData(observed = dat.dz, type = "raw")
  datAdopt <- mxData(observed = dat.adopt, type = "raw")
  
  funML <- mxFitFunctionML()
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt", means = "meanP", dimnames = varnames)
  
  # create model objects for multiple groups
  pars <- list(meanP, covA, covD, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, datMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, datDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, datAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ", "Adopt"))
  
  # create algebra for variance components
  rowVC <- rep("cVC", nv)
  colVC <- rep(c("A", "D", "E", "SA", "SD", "SE"), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, D, E, A/V, D/V, E/V),
                     name = "cVC", dimnames = list(rowVC, colVC))
  # CIs
  ciADE <- mxCI("cVC[1,1:3]")
  
  #build model with confidence intervals
  modelADE <- mxModel("oneADEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciADE)
  fitADE <- mxRun(modelADE, intervals=TRUE) #intervals = estimate CIs
  
  # returns
  correlations <- rbind(mz = mzcor, dz = dzcor, adopt = adopcor)
  returns <- list(model = modelADE, fit = fitADE, cors = correlations)
  return(returns)
}

threeACDEVC <- function(data, var) {
  # variables for analysis
  nv <- 1
  ntv <- 2
  varnames <- paste(var, sep = "",
                    c(rep(1, nv), rep(2, nv)))
  
  mzdata <- subset(data, zyg1 == 1)
  dzdata <- subset(data, zyg1 == 2)
  adoptdata <- subset(data, zyg1 == 0)
  
  # data for analysis
  dat.mz <- select(mzdata, all_of(varnames)) 
  dat.dz <- select(dzdata, all_of(varnames))
  dat.adopt <- select(adoptdata, all_of(varnames))
  
  # statistics
  meanlab <- colMeans(dat.dz, na.rm = TRUE)
  mzcov <- cov(dat.mz, use = "complete")
  dzcov <- cov(dat.dz, use = "complete")
  adopcov <- cov(dat.adopt, use = "complete")
  mzcor <- cor(dat.mz, use = "complete")
  dzcor <- cor(dat.dz, use = "complete")
  adopcor <- cor(dat.adopt, use = "complete")
  
  # set starting values
  # note: if you get NA errors just set the starting values to 0
  svmu <- meanlab[1]
  svA <- sqrt(2*(mzcov[1,2] - dzcov[1,2]))
  svC <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  svD <- sqrt(2*dzcov[1,2]- mzcov[1,2])
  svE <- sqrt(mzcov[2,2] - mzcov[1,2])
  
  # mean and variance components
  meanP <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svmu, labels = "x1", name = "meanP")
  
  covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svA, labels = "VA", name = "A")
  covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svC, labels = "VC", name = "C")
  covD <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svD, labels = "VD", name = "D")
  covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = svE, labels = "VE", name = "E")
  
  # covariances
  covP <- mxAlgebra(expression = A + C + D + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C + D, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C + 0.25%x%D, name = "cDZ")
  covAdopt <- mxAlgebra(expression = C, name = "cAdopt")
  
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(t(cAdopt), V)), name = "expCovAdopt")
  
  # create data objects for multiple groups
  datMZ <- mxData(observed = dat.mz, type = "raw")
  datDZ <- mxData(observed = dat.dz, type = "raw")
  datAdopt <- mxData(observed = dat.adopt, type = "raw")
  
  funML <- mxFitFunctionML()
  
  # create expectation objects for multiple groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanP", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanP", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt", means = "meanP", dimnames = varnames)
  
  # create model objects for multiple groups
  pars <- list(meanP, covA, covC, covD, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, datMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, datDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, datAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ", "Adopt"))
  
  # create algebra for variance components
  rowVC <- rep("cVC", nv)
  colVC <- rep(c("A", "C", "D", "E", "SA", "SC", "SD", "SE"), each = nv)
  estVC <- mxAlgebra(expression = cbind(A, C, D, E, A/V, C/V, D/V, E/V),
                     name = "cVC", dimnames = list(rowVC, colVC))
  # CIs
  ciACDE <- mxCI("cVC[1,1:4]")
  
  #build model with confidence intervals
  modelACDE <- mxModel("oneACDEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciACDE)
  fitACDE <- mxTryHard(modelACDE, intervals=TRUE) #intervals = estimate CIs
  
  # returns
  correlations <- list(mz = mzcor, dz = dzcor, adopt = adopcor)
  covs <- list(mz = mzcov, dz = dzcov, adopt = adopcov)
  returns <- list(model = modelACDE, fit = fitACDE, covs = covs, cors = correlations)
  return(returns)
}
  
#------submodels--------

submodelsACDEpath <- function(modelfit) {
  modelACE <- mxModel(modelfit, name="ACEd")
  modelACE <- omxSetParameters(modelACE, labels = "d11", free = FALSE, values = 0)
  fitACE <- mxRun(modelACE, intervals=TRUE)
  
  modelADE <- mxModel(modelfit, name = "ADEd")
  modelADE <- omxSetParameters(modelADE, labels = "c11", free = FALSE, values = 0)
  fitADE <- mxRun(modelADE, intervals = TRUE)
  
  modelCDE <- mxModel(modelfit, name = "CDEd")
  modelCDE <- omxSetParameters(modelCDE, labels = "a11", free = FALSE, values = 0)
  fitCDE <- mxRun(modelCDE, intervals = TRUE)
  
  modelAE <- mxModel(modelfit, name = "AEd")
  modelAE <- omxSetParameters(modelAE, labels = c("c11", "d11"), free = FALSE, values = 0)
  fitAE <- mxRun(modelAE, intervals = TRUE)
  
  modelDE <- mxModel(modelfit, name = "DEd")
  modelDE <- omxSetParameters(modelDE, labels = c("a11", "c11"), free = FALSE, values = 0)
  fitDE <- mxRun(modelDE, intervals = TRUE)
  
  modelCE <- mxModel(modelfit, name = "CEd")
  modelCE <- omxSetParameters(modelCE, labels = c("a11", "d11"), free = FALSE, values = 0)
  fitCE <- mxRun(modelCE, intervals = TRUE)
  
  modelE <- mxModel(modelfit, name = "Ed")
  modelE <- omxSetParameters(modelE, labels = c("a11", "c11", "d11"), free = FALSE, values = 0)
  fitE <- mxRun(modelE, intervals = TRUE)
  
  compare <- mxCompare(modelfit, 
                       c(fitACE, fitADE, fitCDE, fitAE, 
                         fitDE, fitCE, fitE))
  returns <- list(modelACE = modelACE, fitACE = fitACE, 
                  modelADE = modelADE, fitADE = fitADE, 
                  modelCDE = modelCDE, fitCDE = fitCDE, 
                  modelAE = modelAE, fitAE = fitAE, 
                  modelDE = modelDE, fitDE = fitDE, 
                  modelCE = modelCE, fitCE = fitCE, 
                  modelE = modelE, fitE = fitE, 
                  compare = compare)

  return(returns)
}

submodelsACDEVC <- function(modelfit) {
  modelACE <- mxModel(modelfit, name="ACEd")
  modelACE <- omxSetParameters(modelACE, labels = "VD", free = FALSE, values = 0)
  fitACE <- mxRun(modelACE, intervals=TRUE)
  
  modelADE <- mxModel(modelfit, name = "ADEd")
  modelADE <- omxSetParameters(modelADE, labels = "VC", free = FALSE, values = 0)
  fitADE <- mxRun(modelADE, intervals = TRUE)
  
  modelCDE <- mxModel(modelfit, name = "CDEd")
  modelCDE <- omxSetParameters(modelCDE, labels = "VA", free = FALSE, values = 0)
  fitCDE <- mxRun(modelCDE, intervals = TRUE)
  
  modelAE <- mxModel(modelfit, name = "AEd")
  modelAE <- omxSetParameters(modelAE, labels = c("VC", "VD"), free = FALSE, values = 0)
  fitAE <- mxRun(modelAE, intervals = TRUE)
  
  modelCE <- mxModel(modelfit, name = "CEd")
  modelCE <- omxSetParameters(modelCE, labels = c("VA", "VD"), free = FALSE, values = 0)
  fitCE <- mxRun(modelCE, intervals = TRUE)
  
  modelE <- mxModel(modelfit, name = "Ed")
  modelE <- omxSetParameters(modelE, labels = c("VA", "VC", "VD"), free = FALSE, values = 0)
  fitE <- mxRun(modelE, intervals = TRUE)
  
  compare <- mxCompare(modelfit, 
                       c(fitACE, fitADE, 
                         #fitCDE, 
                         fitAE, 
                         fitCE, fitE))
  returns <- list(modelACE = modelACE, fitACE = fitACE, 
                  modelADE = modelADE, fitADE = fitADE, 
                  modelCDE = modelCDE, fitCDE = fitCDE, 
                  modelAE = modelAE, fitAE = fitAE, 
                  modelCE = modelCE, fitCE = fitCE, 
                  modelE = modelE, fitE = fitE, 
                  compare = compare)
  
  return(returns)
}

umxsubmodelsACE <- function(model) {
  AE <- umxModify(model, update = "C_r1c1", comparison = T, name = "AE")
  CE <- umxModify(model, update = "A_r1c1", comparison = T, name = "CE")
  E <- umxModify(model, update = c("A_r1c1", "C_r1c1"), comparison = T, name = "E")
  compare <- umxCompare(model, c(AE, CE, E))
  return(list(AE = AE, CE = CE, E = E, compare = compare))
}
