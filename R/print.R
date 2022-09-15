print_err <- function(traj, truth, maxObsTime) {
  iestitru <- apply(abs(outer(truth$time, traj$time, `-`)), 1, which.min)
  err_traj <- sqrt(rowSums((traj$state[iestitru,] - truth$state)^2))

  message("obs err: ", mean(err_traj[truth$time <= maxObsTime]))
  message("extrap err: ", mean(err_traj[truth$time > maxObsTime]))
}
