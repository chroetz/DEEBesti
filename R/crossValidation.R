validate <- function(
    obsTrain,
    obsVali,
    hyperParms,
    memoize,
    opts
  ) {

  parms <- getParms(
    obsTrain,
    hyperParms,
    memoize = memoize)

  esti <- solveOde(
    u0 = getInitialState(parms$trajs),
    fun = buildDerivFun(hyperParms$derivFun),
    timeRange = c(0, max(obsVali$time)),
    opts = opts$odeSolver,
    parms = parms)

  obsVali <- parms$normalization$normalize(obsVali)

  cleanUpParms(parms)

  err <- l2err(esti, obsVali) # TODO: make error type an option
  return(err)
}


l2err <- function(trajs, obs) {
  # TODO: apply per trajId
  trajsObs <- interpolateTrajs(trajs, obs$time)
  mean((obs$state - trajsObs$state)^2)
}


validateHyperparams <- function(
    obsTrain,
    obsVali,
    hyperParmsSet,
    opts
  ) {
  prepareMemory()
  vali <- sapply(
    seq_len(length(hyperParmsSet)),
    \(i) validate(
      obsTrain,
      obsVali,
      hyperParmsSet[[i]],
      memoize = TRUE,
      opts))
  vali[is.na(vali)] <- Inf
  return(vali)
}


