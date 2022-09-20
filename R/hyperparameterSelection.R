validateHyperparams <- function(obsTrain, obsVali, hyperParmsSet, method) {
  validateFun <- getValidateFun(method)
  prepareValidationMemory(method, nrow(hyperParmsSet))
  vali <- sapply(
    seq_len(nrow(hyperParmsSet)),
    \(i) validateFun(obsTrain, obsVali, hyperParmsSet[i, ]))
  vali[is.na(vali)] <- Inf
  dplyr::bind_cols(hyperParmsSet, validationErr = vali)
}


selectHyperparams <- function(obs, hyperParmsSet, method) {
  n <- getCount(obs)
  iVali <- 1:floor(n/5) * 5
  iTrain <- setdiff(1:n, iVali)
  obsVali <- obs[iVali,]
  obsTrain <- obs[iTrain,]

  pt <- proc.time()
  validation <- validateHyperparams(obsTrain, obsVali, hyperParmsSet, method=method)
  message(as.vector((proc.time()-pt)["elapsed"]), "s")

  return(validation[which.min(validation$validationErr),])
}


estimateWithHyperparameterSelection <- function(obs, hyperParmsSet, method, outTimes) {
  optiHyperParms <- selectHyperparams(obs, hyperParmsSet, method = method)
  res <- getParmsAndIntitialState(obs, optiHyperParms, method = method)
  trajFinal <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(optiHyperParms$derivFun),
    times = outTimes,
    parms = res$parms)
  return(trajFinal)
}

getParmsAndIntitialState <- function(obs, hyperParms, method) {
  if (method == "Colloc") {
    smoothed <- estimateParmsColloc(obs, hyperParms$bwTime, hyperParms$kernelTime)
    parms <- as.list(smoothed)
    parms$bw <- hyperParms$bwState
    parms$kernel <- getKernel(hyperParms$kernelState)
    initialState <- getInitialState(smoothed)
  } else if (method == "Altopi") {
    preTraj <- getPreviousAltopiTraj(obs, hyperParms)
    parms <- oneAltopiStep(
      preTraj,
      obs,
      gamma = hyperParms$gamma,
      fitDeriv = fitLocalConst,
      fitDerivOpts = list(bw = hyperParms$bw, kernel = getKernel(hyperParms$kernel)),
      fitTraj = updateAltopiTraj)
    initialState <- getInitialState(parms)
  } else {
    stop("Unknown method ", method)
  }
  return(list(parms = parms, initialState = initialState))
}

