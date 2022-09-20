initAltopi <- function(obs, hyperPrams) {
  n <- nrow(obs$state)
  steps <- (n-1)*hyperPrams$interSteps + 1
  time <- seq(min(obs$time), max(obs$time), length.out = steps)
  traj <- interpolateTrajs(obs, time)
  traj <- setDeriv(traj)
  return(traj)
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
  if (hyperParms$S <= 0) return(initAltopi(obs, hyperParms))
  if (memoize) {
    traj <- getFromAltopiMemory(hyperParms)
    if (!is.null(traj)) return(traj)
  }
  preHyperParms <- getHyperParmPredecessorAltopi(hyperParms)
  preTraj <- getAltopiTraj(obs, preHyperParms, memoize)
  nextTraj <- oneAltopiStep(
    preTraj,
    obs,
    gamma = hyperParms$gamma,
    fitDeriv = getFitter(hyperParms$fitter),
    fitDerivOpts = list(bw = hyperParms$fitterBw, kernel = getKernel(hyperParms$fitterKernel)),
    fitTraj = updateAltopiTraj)
  if (memoize) {
    addToAltopiMemory(hyperParms, nextTraj)
  }
  return(nextTraj)
}


altopiMemory <- NULL

getHyperParmPredecessorAltopi <- function(hyperParms) {
  hyperParms$S <- hyperParms$S - 1
  return(hyperParms)
}

addToAltopiMemory <- function(hyperParms, trajs) {
  utils::sethash(altopiMemory, hyperParms, trajs)
}

getFromAltopiMemory <- function(hyperParms) {
  utils::gethash(altopiMemory, hyperParms, nomatch = NULL)
}

