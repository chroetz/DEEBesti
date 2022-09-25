initAltopi <- function(obs, interSteps) {
  trajIds <- getTrajIds(obs)
  trajsList <- lapply(trajIds, \(trajId) {
    trj <- getTrajsWithId(obs, trajId)
    n <- getCount(trj)
    steps <- (n-1) * interSteps + 1
    time <- seq(min(trj$time), max(trj$time), length.out = steps)
    trj <- interpolateTrajs(trj, time)
    trj <- setDeriv(trj)
    trj
  })
  return(bindTrajs(trajsList))
}


oneAltopiStep <- function(
    trajs, obs, gamma,
    fitDeriv, fitDerivOpts=NULL,
    fitTraj, fitTrajOpts=NULL
  ) {
  newTrajs <- trajs
  newTrajs$deriv <- do.call(fitDeriv, c(list(x=trajs$state, y=trajs$deriv), fitDerivOpts))
  newTrajs$state <- do.call(fitTraj, c(list(traj=newTrajs, obs=obs, gamma=gamma), fitTrajOpts))
  return(newTrajs)
}


getAltopiTraj <- function(obs, hyperParms, opts, memoize = FALSE) {
  opts <- asOpts(opts, c("Altopi", "Method"))
  hyperParms <- asOpts(hyperParms, c("Altopi", "HyperParms"))
  if (hyperParms$steps <= 0) return(initAltopi(obs, opts$interSteps))
  if (memoize) {
    traj <- getFromAltopiMemory(hyperParms)
    if (!is.null(traj)) return(traj)
  }
  preHyperParms <- getHyperParmPredecessorAltopi(hyperParms)
  preTraj <- getAltopiTraj(obs, preHyperParms, opts, memoize)
  nextTraj <- oneAltopiStep(
    preTraj,
    obs,
    gamma = hyperParms$gamma,
    fitDeriv = buildFitter(hyperParms$fitter),
    fitTraj = updateAltopiTraj)
  if (memoize) {
    addToAltopiMemory(hyperParms, nextTraj)
  }
  return(nextTraj)
}


altopiMemory <- utils::hashtab()

getHyperParmPredecessorAltopi <- function(hyperParms) {
  hyperParms$steps <- hyperParms$steps - 1
  return(hyperParms)
}

addToAltopiMemory <- function(hyperParms, trajs) {
  utils::sethash(altopiMemory, hyperParms, trajs)
}

getFromAltopiMemory <- function(hyperParms) {
  utils::gethash(altopiMemory, hyperParms, nomatch = NULL)
}

