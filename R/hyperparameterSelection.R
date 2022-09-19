validateHyperparams <- function(obsTrain, obsVali, hyperParms, validateFun) {
  clearMemory()
  vali <- sapply(
    seq_len(nrow(hyperParms)),
    \(i) validateFun(obsTrain, obsVali, hyperParms[i, ]))
  vali[is.na(vali)] <- Inf
  dplyr::bind_cols(hyperParms, validationErr = vali)
}


selectHyperparams <- function(obs, hyperParms, validateFun) {
  n <- getCount(obs)
  iVali <- 1:floor(n/5) * 5
  iTrain <- setdiff(1:n, iVali)
  obsVali <- obs[iVali,]
  obsTrain <- obs[iTrain,]

  pt <- proc.time()
  validation <- validateHyperparams(obsTrain, obsVali, hyperParms, validateFun=validateFun)
  message(as.vector((proc.time()-pt)["elapsed"]), "s")

  return(validation[which.min(validation$validationErr),])
}


estimateWithHyperparameterSelectionAltopi <- function(obs, hyperParms) {
  optiHyperParms <- selectHyperparams(obs, hyperParms, validateFun=validateAltopi)

  trajsInit <- initAltopi(obs, interSteps = optiHyperParms$interSteps)
  trajsOpti <- optimizeTrajs(trajsInit, obs, gamma=optiHyperParms$gamma, S=optiHyperParms$S, bw=optiHyperParms$bw, optiHyperParms$kernel)

  tMax <- max(obs$time)*2
  trajFinal <- solveOde(
    u0 = getInitialState(trajsOpti),
    fun = buildDerivFun(optiHyperParms$derivFun),
    tStep = tMax / 1e3,
    tMax = tMax,
    parms = trajsOpti)

  return(trajFinal)
}

estimateWithHyperparameterSelectionColloc <- function(obs, hyperParms) {
  optiHyperParms <- selectHyperparams(obs, hyperParms, validateFun=validateColloc)

  smoothed <- estimateParmsColloc(obs, optiHyperParms$bwTime, optiHyperParms$kernelTime)
  tMax <- max(obs$time)*2
  trajFinal <- estimateTrajsColloc(
    smoothed, optiHyperParms$bwState, optiHyperParms$kernelState,
    tMax = tMax,
    tStep = tMax / 1e3)

  return(trajFinal)
}

