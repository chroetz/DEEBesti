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

  startTime <- opts$startTime * max(obsTrain$time)

  esti <- solveOde(
    u0 = getInitialState(parms$trajs, startTime),
    fun = buildDerivFun(hyperParms$derivFun),
    timeRange = c(startTime, max(obsVali$time)),
    opts = opts$odeSolver,
    parms = parms)

  obsVali <- parms$normalization$normalize(obsVali)

  cleanUpParms(parms)

  errTraj <- lpErr(parms$trajs, obsTrain, opts$errorPower)
  errVali <- lpErr(esti, obsVali, opts$errorPower)
  errTrain <- lpErr(esti, obsTrain[obsTrain$time>=startTime,], opts$errorPower)
  err <-
    opts$valiErrorWeight * errVali +
    opts$trainErrorWeight * errTrain +
    opts$trajsErrorWeight * errTraj
  return(err)
}


lpErr <- function(trajs, target, p) {
  if (
    is.null(trajs) ||
    nrow(trajs) == 0 ||
    !all(is.finite(trajs$state))
  ) return(Inf)
  trajsObs <- interpolateTrajs(trajs, target$time)
  errs <- apply2TrajId(
    trajsObs,
    target,
    \(x, y, p) mean(abs(y$state - x$state)^p),
    p = p,
    simplify = TRUE)
  return(mean(errs))
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


