vali_err <- function(z, fun) {
  traj <- solveOde(z$state[1,], fun, 0.1, max(c(obs_vali$time, obs_train$time)), parms=z)
  i <- apply(abs(outer(traj$time, obs_vali$time, `-`)), 2, which.min)
  #mean((obs_vali$u - traj$state[i, ])^2 * traj$time[i]) # emphasize later points
  mean((obs_vali$state - traj$state[i, ])^2)
}

validate <- function(z, obs, bw, gamma, Smax, kernel=stats::dnorm, derivFun = derivFunNearestNeighbor) {
  errs <- double(Smax)
  for(i in seq_len(Smax)) {
    z <- step_optimization(
      z, obs, gamma,
      fit_deriv = fitLocalConst, fit_deriv_opts = list(bw = bw, kernel = kernel),
      fit_traj = update_trajectory)
    errs[i] <- vali_err(z, derivFun)
  }
  return(errs)
}
