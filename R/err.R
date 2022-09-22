l2err <- function(trajs, obs) {
  # TODO: apply per trajId
  trajsObs <- interpolateTrajs(trajs, obs$time)
  mean((obs$state - trajsObs$state)^2)
}
