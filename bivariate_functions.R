library(OpenMx)
library(dplyr)
source("miFunctions.R")

# in this function script: bivariate ACDE via VC/cholesky with submodels

bivariateACDEVC <- function(data, vars) {
  # Select Variables for Analysis
  nv <- 2                         # number of variables
  ntv <- nv*2                      # number of total variables
  varnames <- paste(vars, c(rep("1",nv), rep("2",nv)), sep="")
  
  # Select Data for Analysis
  mzdat <- subset(data, zyg1 == 1, varnames)
  dzdat <- subset(data, zyg1 == 2, varnames)
  adopdat <- subset(data, zyg1 == 0, varnames)
  
  # Generate Descriptive Statistics
  colMeans(mzdat, na.rm = TRUE)
  colMeans(dzdat, na.rm = TRUE)
  cov(mzdat, use = "complete")
  cov(dzdat, use = "complete")
  mzcor <- cor(mzdat, use = "complete")
  dzcor <- cor(dzdat, use = "complete")
  adoptcor <- cor(adopdat, use = "complete")
  
  # Set Starting Values 
  svMe <- c(.01, .1)                # start value for means
  svPa <- .4                        # start value for path coefficient
  svPc <- .2
  svPe <- .4                        # start value for path coefficient for e
  
  #-----PREPARE MODEL-------
  
  # Create Algebra for expected Mean Matrices
  meanG <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svMe, labels = labVars("mean", vars), name = "meanG")
  
  # Create Matrices for Variance Components
  covA <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = valDiag(svPa,nv), labels = labLower("VA",nv), name = "VA") 
  covC <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = valDiag(svPc,nv), labels = labLower("VC",nv), name = "VC")
  covD <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = valDiag(svPc,nv), labels = labLower("VD",nv), name = "VD")
  covE <- mxMatrix(type = "Symm", nrow = nv, ncol = nv, free = TRUE, 
                   values = valDiag(svPe,nv), labels = labLower("VE",nv), name = "VE")
  
  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covP <- mxAlgebra(expression = VA + VC + VD + VE, name = "V")
  covMZ <- mxAlgebra(expression = VA + VC + VD, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%VA + 0.25%x%VD + VC, name = "cDZ")
  covAdopt <- mxAlgebra(expression = VC, name = "cAdopt")
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(t(cAdopt), V)), name = "expCovAdopt")
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData(observed = mzdat, type = "raw")
  dataDZ <- mxData(observed = dzdat, type = "raw")
  dataAdopt <- mxData(observed = adopdat, type = "raw")
  
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanG", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanG", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt", mean = "meanG", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  pars <- list(meanG, covA, covC, covD, covE, covP)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, dataAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ","DZ","Adopt"))
  
  # Create Algebra for Standardization
  matI <- mxMatrix(type="Iden", nrow = nv, ncol = nv, name = "I")
  invSD <- mxAlgebra(expression = solve(sqrt(I*V)), name = "iSD")
  
  # Calculate genetic and environmental correlations
  corA <- mxAlgebra(expression = solve(sqrt(I*VA))%&%VA, name ="rA") #cov2cor()
  corC <- mxAlgebra(expression = solve(sqrt(I*VC))%&%VC, name ="rC")
  corD <- mxAlgebra(expression = solve(sqrt(I*VD))%&%VD, name ="rD")
  corE <- mxAlgebra(expression = solve(sqrt(I*VE))%&%VE, name ="rE")
  
  # Calculate Phenotypic Correlation
  corP <- mxAlgebra(expression = solve(sqrt(I*V))%*%V%*%solve(sqrt(I*V)), name = "rP")
  
  # Create Algebra for Variance Components
  rowUV <- rep('UV', nv)
  colUV <- rep(c('VA','VC','VD','VE'),each = nv)
  estUV <- mxAlgebra(expression = cbind(VA,VC,VD,VE), 
                     name = "UV", dimnames = list(rowUV,colUV))
  rowSV <- rep('SV', nv)
  colSV <- rep(c('SA','SC','SD','SE'),each = nv)
  estSV <- mxAlgebra(expression = cbind(VA/V,VC/V,VD/V,VE/V), 
                     name = "SV", dimnames = list(rowSV,colSV))
  
  # Create Confidence Interval Objects
  ciUV <- mxCI(c("UV"))
  ciSV <- mxCI(c("SV"))
  
  # Build Model with Confidence Intervals
  calc <- list(matI, invSD, corA, corC, corD, corE, corP, estUV, estSV, ciUV, ciSV)
  modelACDE  <- mxModel("twoACDEvc", pars, modelMZ, modelDZ, modelAdopt, multi, calc)
  
  # run mdatr
  fitACDE <- mxRun(modelACDE, intervals=T)
  cors <- list(mz = mzcor, dz = dzcor, adopt = adoptcor)
  returns <- list(model = modelACDE, fit = fitACDE, cors = cors)
  return(returns)
}

biACDEVCsubs <- function(fit, omit) {
  nv <- 2
  # Run ADE omitting bpd and item c | ACE v ADE
  if('ADEf' %in% omit) {}
  else {
    modelADEf <- mxModel(fit, name = "twoADEvcf")
    modelADEf <- omxSetParameters(modelADEf, labels = labLower("VC",nv), free = FALSE, values = 0)
    fitADEf <- mxRun(modelADEf, intervals = T)
    compnames <- c(fitADEf)
    list <- c(modelADEf, fitADEf)
    listnames <- c("modelADEf", "fitADEf")
  }

  # Run ADE omitting bpd c, keep item c | ADE v ACDE
  if('ADEp' %in% omit) {}
  else {
    modelADEp <- mxModel(fit, name = "twoADEvcp")
    modelADEp <- omxSetParameters(modelADEp, labels = c("VC21","VC22"), free = FALSE, values = 0)
    fitADEp <- mxRun(modelADEp, intervals = T )
    compnames <- c(compnames, fitADEp)
    list <- c(list, modelADEp, fitADEp)
    listnames <- c(listnames, "modelADEp", "fitADEp")
  }
  
  # Run ACE omitting bpd d and item d | ACE v ACE
  if('ACEf' %in% omit) {}
  else {
    modelACEf <- mxModel(fit, name = "twoACEvcf")
    modelACEf <- omxSetParameters(modelACEf, labels = labLower("VD",nv), free = FALSE, values = 0)
    fitACEf <- mxRun(modelACEf, intervals = T)
    compnames <- c(compnames, fitACEf)
    list <- c(list, modelACEf, fitACEf)
    listnames <- c(listnames, "modelACEf", "fitACEf")
  }
  
  # Run ACE omitting bpd d, keep item d | ACE v ACDE
  if('ACEp' %in% omit) {}
  else {
    modelACEp <- mxModel(fit, name = "twoACEvcp")
    modelACEp <- omxSetParameters(modelACEp, labels = c("VD21","VD22"), free = FALSE, values = 0)
    fitACEp <- mxRun(modelACEp, intervals = T)
    compnames <- c(compnames, fitACEp)
    list <- c(list, modelACEp, fitACEp)
    listnames <- c(listnames, "modelACEp", "fitACEp")
  }
  
  # Run AE omitting bpd c/d, keeping item c/d | AE v ACDE
  if('AEpp' %in% omit) {}
  else {
    modelAEpp <- mxModel(fit, name = "twoAEvcpp")
    modelAEpp <- omxSetParameters(modelAEpp, labels = c("VC21","VC22"), free = FALSE, values = 0)
    modelAEpp <- omxSetParameters(modelAEpp, labels = c("VD21","VD22"), free = FALSE, values = 0)
    fitAEpp <- mxRun(modelAEpp, intervals = T)
    compnames <- c(compnames, fitAEpp)
    list <- c(list, modelAEpp, fitAEpp)
    listnames <- c(listnames, "modelAEpp", "fitAEpp")
  }
  
  # Run AE omitting bpd c/d and item d, keeping item c | AE vs ACE
  if('AEpf' %in% omit) {}
  else {
    modelAEpf <- mxModel(fit, name = "twoAEvcpf")
    modelAEpf <- omxSetParameters(modelAEpf, labels = c("VC21", "VC22"), free = FALSE, values = 0)
    modelAEpf <- omxSetParameters(modelAEpf, labels = labLower("VD",nv), free = FALSE, values = 0)
    fitAEpf <- mxRun(modelAEpf, intervals = T)
    compnames <- c(compnames, fitAEpf)
    list <- c(list, modelAEpf, fitAEpf)
    listnames <- c(listnames, "modelAEpf", "fitAEpf")
  }
  
  # Run AE omitting bpd c/d and item c, keeping item d | AE vs ADE
  if('AEfp' %in% omit) {}
  else {
    modelAEfp <- mxModel(fit, name="twoAEvcfp")
    modelAEfp <- omxSetParameters(modelAEfp, labels = labLower("VC",nv), free = FALSE, values = 0)
    modelAEfp <- omxSetParameters(modelAEfp, labels = c("VD21", "VD22"), free = FALSE, values = 0)
    fitAEfp <- mxRun(modelAEfp, intervals = T )
    compnames <- c(compnames, fitAEfp)
    list <- c(list, modelAEfp, fitAEfp)
    listnames <- c(listnames, "modelAEfp", "fitAEfp")
  }
  
  # Run AE omitting both bpd and item c/d | AE v AE
  if('AEff' %in% omit) {}
  else {
    modelAEff <- mxModel(fit, name="twoAEvcff")
    modelAEff <- omxSetParameters(modelAEff, labels = labLower("VC",nv), free = FALSE, values = 0)
    modelAEff <- omxSetParameters(modelAEff, labels = labLower("VD",nv), free = FALSE, values = 0)
    fitAEff <- mxRun(modelAEff, intervals = T)
    compnames <- c(compnames, fitAEff)
    list <- c(list, modelAEff, fitAEff)
    listnames <- c(listnames, "modelAEff", "fitAEff")
  }
  compare <- mxCompare(fit, compnames)

  list <- c(list, list(compare))
  names(list) <- c(listnames, "compare")
  return(list)
}

bivariateACDEpath <- function(vars, data) {
  # Select Variables for Analysis
  nv <- 2                         # number of variables
  ntv <- nv*2                      # number of total variables
  varnames <- paste(vars, c(rep("1",nv), rep("2",nv)), sep="")
  
  # Select Data for Analysis
  mzdat <- subset(data, zyg1 == 1, varnames)
  dzdat <- subset(data, zyg1 == 2, varnames)
  adopdat <- subset(data, zyg1 == 0, varnames)
  
  # Generate Descriptive Statistics
  colMeans(mzdat, na.rm = TRUE)
  colMeans(dzdat, na.rm = TRUE)
  mzcov <- cov(mzdat, use = "complete")
  dzcov <- cov(dzdat, use = "complete")
  mzcor <- cor(mzdat, use = "complete")
  dzcor <- cor(dzdat, use = "complete")
  adoptcor <- cor(adopdat, use = "complete")
  
  # Set Starting Values 
  svMe <- c(.01, .1)                # start value for means
  svpA <- .3
  svPaD <- vech(diag(svpA,nv,nv)) 
  svpC <- .2
  svPcD <- vech(diag(svpC,nv,nv)) 
  svpD <- .2
  svPdD <- vech(diag(svpD,nv,nv)) 
  svpE <- .3
  svPeD <- vech(diag(svpE,nv,nv)) 
  
  # create mean matrices
  meanG <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                    values = svMe, labels = labVars("mean", vars), name = "meanG")
  
  # create path coef matrices
  coefA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svPaD, labels = labLower("a",nv), name = "a")
  coefC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svPcD, labels = labLower("c",nv), name = "c")
  coefD <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svPdD, labels = labLower("d",nv), name = "d")
  coefE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, free = TRUE, 
                    values = svPeD, labels = labLower("e",nv), name = "e")
  # create cov matrices
  covA <- mxAlgebra(expression = a %*% t(a), name = "A")
  covC <- mxAlgebra(expression = c %*% t(c), name = "C")
  covD <- mxAlgebra(expression = d %*% t(d), name = "D")
  covE <- mxAlgebra(expression = e %*% t(e), name = "E")
  
  covP <- mxAlgebra(expression = A + C + D + E, name = "V")
  covMZ <- mxAlgebra(expression = A + C + D, name = "cMZ")
  covDZ <- mxAlgebra(expression = 0.5%x%A + C + 0.25%x%D, name = "cDZ")
  covAdopt <- mxAlgebra(expression = C, name = "cAdopt")
  
  # Create Algebra for Standardization
  matI <- mxMatrix(type="Iden", nrow = nv, ncol = nv, name = "I")
  invSD <- mxAlgebra(expression = solve(sqrt(I*V)), name = "iSD")
  
  # Calculate genetic and environmental correlations
  corA <- mxAlgebra(expression = solve(sqrt(I*A))%&%A, name ="rA") #cov2cor()
  corC <- mxAlgebra(expression = solve(sqrt(I*C))%&%C, name ="rC")
  corD <- mxAlgebra(expression = solve(sqrt(I*D))%&%D, name ="rD")
  corE <- mxAlgebra(expression = solve(sqrt(I*E))%&%E, name ="rE")
  
  # expectation cov matrices
  expCovMZ <- mxAlgebra(expression = rbind(cbind(V, cMZ), cbind(cMZ, V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(expression = rbind(cbind(V, cDZ), cbind(cDZ, V)), name = "expCovDZ")
  expCovAdopt <- mxAlgebra(expression = rbind(cbind(V, cAdopt), cbind(cAdopt, V)), name = "expCovAdopt")
  
  # data matrices
  dataMZ <- mxData(observed = mzdat, type = "raw")
  dataDZ <- mxData(observed = dzdat, type = "raw")
  dataAdopt <- mxData(observed = adopdat, type = "raw")
  
  #expectation matrices for each group
  expMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "meanG", dimnames = varnames)
  expDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "meanG", dimnames = varnames)
  expAdopt <- mxExpectationNormal(covariance = "expCovAdopt", means = "meanG", dimnames = varnames)
  funML <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  pars <- list(meanG, matI, invSD, coefA, coefC, coefD, coefE, 
               covA, covC, covD, covE, covP, corA, corC, corD, corE)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name = "DZ")
  modelAdopt <- mxModel(pars, covAdopt, expCovAdopt, dataAdopt, expAdopt, funML, name = "Adopt")
  multi <- mxFitFunctionMultigroup(c("MZ","DZ", "Adopt"))
  
  # Create Algebra for Variance Components
  colVC <- vars
  rowVC <- rep(c('A','C','D','E','SA','SC','SD','SE'),each=nv)
  estVC <- mxAlgebra( expression=rbind(A,C,D,E,A/V,C/V,D/V,E/V), name="VC", dimnames=list(rowVC,colVC))
  
  # Create Confidence Interval Objects
  ciVC <- mxCI("VC")
  
  # Build Model with Confidence Intervals
  modelACDE <- mxModel("biACDEc", pars, modelMZ, modelDZ, modelAdopt, multi, estVC, ciVC)
  
  # run mdatr
  fitACDE <- mxRun(modelACDE, intervals=T)
  cors <- list(mz = mzcor, dz = dzcor, adopt = adoptcor)
  list <- list(model = modelACDE, fit = fitACDE, cors = cors)
  return(list)
}

biACDEpathsubs <- function(fit, omit) {
  nv <- 2
  # Run ADE omitting bpd and item c | ACE v ADE
  compnames <- c()
  if('ADEf' %in% omit) {}
  else {
    modelADEf <- mxModel(fit, name = "twoADEcf")
    modelADEf <- omxSetParameters(modelADEf, labels = labLower("c",nv), free = FALSE, values = 0)
    fitADEf <- mxRun(modelADEf, intervals = T)
    compnames <- c(fitADEf)
    list <- c(modelADEf, fitADEf)
    listnames <- c("modelADEf", "fitADEf")
  }
  
  # Run ADE omitting bpd c, keep item c | ADE v ACDE
  if('ADEp' %in% omit) {}
  else {
    modelADEp <- mxModel(fit, name="twoADEcp")
    modelADEp <- omxSetParameters(modelADEp, labels = c("c21", "c22"), free = FALSE, values = 0)
    fitADEp <- mxRun(modelADEp, intervals = T)
    compnames <- c(compnames, fitADEp)
    list <- c(list, modelADEp, fitADEp)
    listnames <- c(listnames, "modelADEp", "fitADEp")
  }
  
  # Run ACE omitting bpd d and item d | ACE v ACE
  if('ACEf' %in% omit) {}
  else {
    modelACEf <- mxModel(fit, name="twoACEcf")
    modelACEf <- omxSetParameters(modelACEf, labels = labLower("d",nv), free = FALSE, values = 0)
    fitACEf <- mxRun(modelACEf, intervals = T)
    compnames <- c(compnames, fitACEf)
    list <- c(list, modelACEf, fitACEf)
    listnames <- c(listnames, "modelACEf", "fitACEf")
  }
  
  # Run ACE omitting bpd d, keep item d | ACE v ACDE
  if('ACEp' %in% omit) {}
  else {
    modelACEp <- mxModel(fit, name="twoACEcp")
    modelACEp <- omxSetParameters(modelACEp, labels = c("d21", "d22"), free = FALSE, values = 0)
    fitACEp <- mxRun(modelACEp, intervals = T)
    compnames <- c(compnames, fitACEp)
    list <- c(list, modelACEp, fitACEp)
    listnames <- c(listnames, "modelACEp", "fitACEp")
  }
  
  # Run AE omitting bpd c/d, keeping item c/d | AE v ACDE
  if('AEpp' %in% omit) {}
  else {
    modelAEpp <- mxModel(fit, name="twoAEcpp")
    modelAEpp <- omxSetParameters(modelAEpp, labels = c("c21", "c22"), free = FALSE, values = 0)
    modelAEpp <- omxSetParameters(modelAEpp, labels = c("d21", "d22"), free = FALSE, values = 0)
    fitAEpp <- mxRun(modelAEpp, intervals = T)
    compnames <- c(compnames, fitAEpp)
    list <- c(list, modelAEpp, fitAEpp)
    listnames <- c(listnames, "modelAEpp", "fitAEpp")
  }
  
  # Run AE omitting bpd c/d and item d, keeping item c | AE vs ACE
  if('AEpf' %in% omit) {}
  else {
    modelAEpf <- mxModel(fit, name="twoAEcpf")
    modelAEpf <- omxSetParameters(modelAEpf, labels = c("c21", "c22"), free = FALSE, values = 0)
    modelAEpf <- omxSetParameters(modelAEpf, labels = labLower("d",nv), free = FALSE, values = 0)
    fitAEpf <- mxRun(modelAEpf, intervals = T)
    compnames <- c(compnames, fitAEpf)
    list <- c(list, modelAEpf, fitAEpf)
    listnames <- c(listnames, "modelAEpf", "fitAEpf")
  }
  
  # Run AE omitting bpd c/d and item c, keeping item d | AE vs ADE
  if('AEfp' %in% omit) {}
  else {
    modelAEfp <- mxModel(fit, name="twoAEcfp")
    modelAEfp <- omxSetParameters(modelAEfp, labels = labLower("c",nv), free = FALSE, values = 0)
    modelAEfp <- omxSetParameters(modelAEfp, labels = c("d21", "d22"), free = FALSE, values = 0)
    fitAEfp <- mxRun(modelAEfp, intervals = T )
    compnames <- c(compnames, fitAEfp)
    list <- c(list, modelAEfp, fitAEfp)
    listnames <- c(listnames, "modelAEfp", "fitAEfp")
  }
  
  # Run AE omitting both bpd and item c/d | AE v AE
  if('AEff' %in% omit) {}
  else {
    modelAEff <- mxModel(fit, name="twoAEcff")
    modelAEff <- omxSetParameters(modelAEff, labels = labLower("c",nv), free = FALSE, values = 0)
    modelAEff <- omxSetParameters(modelAEff, labels = labLower("d",nv), free = FALSE, values = 0)
    fitAEff <- mxRun(modelAEff, intervals = T)
    compnames <- c(compnames, fitAEff)
    list <- c(list, modelAEff, fitAEff)
    listnames <- c(listnames, "modelAEff", "fitAEff")
  }
  compare <- mxCompare(fit, compnames)
  
  list <- c(list, list(compare))
  names(list) <- c(listnames, "compare")
  return(list)
}
