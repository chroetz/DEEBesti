#' @export
initAltopi <- function(obs, interSteps) {
  n <- nrow(obs$state)
  steps <- (n-1)*interSteps + 1
  time <- seq(min(obs$time), max(obs$time), length.out = steps)
  traj <- interpolateTrajs(obs, time)
  traj <- setDeriv(traj)
  return(traj)
}
