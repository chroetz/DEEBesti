#' @export
stepOptimization <- function(
    trajs, obs, gamma,
    fitDeriv, fitDerivOpts=NULL,
    fitTraj, fitTrajOpts=NULL
  ) {
  newTrajs <- trajs
  newTrajs$deriv <- do.call(fitDeriv, c(list(x=trajs$state, y=trajs$deriv), fitDerivOpts))
  newTrajs$state <- do.call(fitTraj, c(list(traj=newTrajs, obs=obs, gamma=gamma), fitTrajOpts))
  return(newTrajs)
}
