#' @export
initAltopi <- function(obs, hyperPrams) {
  n <- nrow(obs$state)
  steps <- (n-1)*hyperPrams$interSteps + 1
  time <- seq(min(obs$time), max(obs$time), length.out = steps)
  traj <- interpolateTrajs(obs, time)
  traj <- setDeriv(traj)
  return(traj)
}

#' @export
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


getPreviousAltopiTraj <- function(obs, hyperParms) {
  hyperParms$S <- hyperParms$S - 1
  if (hyperParms$S <= 0) return(initAltopi(obs, hyperParms))
  preTraj <- getPreviousAltopiTraj(obs, hyperParms)
  traj <- oneAltopiStep(
    preTraj,
    obs,
    gamma = hyperParms$gamma,
    fitDeriv = fitLocalConst,
    fitDerivOpts = list(bw = hyperParms$bw, kernel = getKernel(hyperParms$kernel)),
    fitTraj = updateAltopiTraj)
  return(traj)
}
