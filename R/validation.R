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

  if (
    is.null(parms$trajs) ||
    nrow(parms$trajs) == 0 ||
    !all(is.finite(parms$trajs$state)) ||
    !all(is.finite(parms$trajs$deriv))
  ) return(Inf)

  esti <- solveOde(
    u0 = getInitialState(parms$trajs),
    fun = buildDerivFun(hyperParms$derivFun),
    timeRange = c(0, max(obsVali$time)),
    opts = opts$odeSolver,
    parms = parms)

  obsVali <- parms$normalization$normalize(obsVali)

  cleanUpParms(parms)

  errVali <- lpErr(esti, obsVali, opts$errorPower)
  errTrain <- lpErr(esti, obsTrain, opts$errorPower)
  err <- (1-opts$trainErrorShare) * errVali + opts$trainErrorShare * errTrain
  return(err)
}


lpErr <- function(trajs, obs, p) {
  if (
    is.null(trajs) ||
    nrow(trajs) == 0 ||
    !all(is.finite(trajs$state))
  ) return(Inf)
  # TODO: apply per trajId
  trajsObs <- interpolateTrajs(trajs, obs$time)
  mean(abs(obs$state - trajsObs$state)^p)
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


