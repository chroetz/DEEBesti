initAltopi <- function(obs, interSteps) {
  mapTrajs2Trajs(obs, \(trj) {
    n <- getCount(trj)
    steps <- (n-1) * interSteps + 1
    time <- seq(min(trj$time), max(trj$time), length.out = steps)
    trj <- interpolateTrajs(trj, time)
    trj <- setDeriv(trj, "center")
    trj
  })
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


getAltopiTraj <- function(obs, hyperParms, memoize = FALSE) {
  hyperParms <- asOpts(hyperParms, c("Altopi", "HyperParms"))
  if (hyperParms$steps <= 0) return(initAltopi(obs, hyperParms$interSteps))
  if (memoize) {
    traj <- getFromMemory(hyperParms)
    if (!is.null(traj)) return(traj)
  }
  preHyperParms <- getHyperParmPredecessorAltopi(hyperParms)
  preTraj <- getAltopiTraj(obs, preHyperParms, memoize)
  nextTraj <- oneAltopiStep(
    preTraj,
    obs,
    gamma = hyperParms$gamma,
    fitDeriv = buildFitter(hyperParms$fitter),
    fitTraj = updateAltopiTraj)
  if (memoize) {
    addToMemory(hyperParms, nextTraj)
  }
  return(nextTraj)
}


getHyperParmPredecessorAltopi <- function(hyperParms) {
  hyperParms$steps <- hyperParms$steps - 1
  return(hyperParms)
}


