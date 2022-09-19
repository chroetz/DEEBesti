validateHyperparams <- function(obsTrain, obsVali, hyperParmsSet, validateFun) {
  clearMemory()
  vali <- sapply(
    seq_len(nrow(hyperParmsSet)),
    \(i) validateFun(obsTrain, obsVali, hyperParmsSet[i, ]))
  vali[is.na(vali)] <- Inf
  dplyr::bind_cols(hyperParmsSet, validationErr = vali)
}


selectHyperparams <- function(obs, hyperParmsSet, validateFun) {
  n <- getCount(obs)
  iVali <- 1:floor(n/5) * 5
  iTrain <- setdiff(1:n, iVali)
  obsVali <- obs[iVali,]
  obsTrain <- obs[iTrain,]

  pt <- proc.time()
  validation <- validateHyperparams(obsTrain, obsVali, hyperParmsSet, validateFun=validateFun)
  message(as.vector((proc.time()-pt)["elapsed"]), "s")

  return(validation[which.min(validation$validationErr),])
}


estimateWithHyperparameterSelection <- function(obs, hyperParmsSet, method) {
  optiHyperParms <- selectHyperparams(obs, hyperParmsSet, validateFun=getValidateFun(method))
  res <- getParmsAndIntitialState(obs, optiHyperParms, method = method)
  tMax <- max(obs$time)*2
  trajFinal <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(optiHyperParms$derivFun),
    tStep = tMax / 1e3,
    tMax = tMax,
    parms = res$parms)
  return(trajFinal)
}

getParmsAndIntitialState <- function(obs, hyperParms, method) {
  if (method == "Colloc") {
    smoothed <- estimateParmsColloc(obs, hyperParms$bwTime, hyperParms$kernelTime)
    parms <- as.list(smoothed)
    parms$bw <- hyperParms$bwState
    parms$kernel <- buildKernel(hyperParms$kernelState)
    initialState <- getInitialState(smoothed)
  } else if (method == "Altopi") {
    trajsInit <- initAltopi(obs, interSteps = hyperParms$interSteps)
    parms <- optimizeTrajs(trajsInit, obs, gamma=hyperParms$gamma, S=hyperParms$S, bw=hyperParms$bw, kernel=hyperParms$kernel)
    initialState <- getInitialState(parms)
  } else {
    stop("Unknown method ", method)
  }
  return(list(parms = parms, initialState = initialState))
}

