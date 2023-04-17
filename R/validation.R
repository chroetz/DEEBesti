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

  esti <- estimateTrajs(
    getInitialState(parms$trajs, startTime),
    c(startTime, max(obsVali$time)),
    parms,
    hyperParms,
    opts)

  obsVali <- parms$normalization$normalize(obsVali)

  cleanUpParms(parms)

  # TODO: use other error measures such as followTime, validTime
  # TODO: re-organize EstiOpts
  errTraj <- lpErr(parms$trajs, obsTrain, opts$errorPower)
  errVali <- lpErr(esti, obsVali, opts$errorPower)
  errTrain <- lpErr(esti, obsTrain[obsTrain$time>=startTime,], opts$errorPower)
  err <-
    opts$valiErrorWeight * errVali +
    opts$trainErrorWeight * errTrain +
    opts$trajsErrorWeight * errTraj
  return(err)
}


lpErr <- function(trajs, targets, p) {
  if (
    is.null(trajs) ||
    nrow(trajs) == 0 ||
    !all(is.finite(trajs$state))
  ) return(Inf)
  errs <- apply2TrajId(
    trajs,
    targets,
    \(traj, target, p) {
      trajInterp <- interpolateTrajs(traj, target$time)
      mean(abs(target$state - trajInterp$state)^p)
    },
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


