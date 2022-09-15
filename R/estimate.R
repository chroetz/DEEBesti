
optimizeTraj <- function(z, obs, gamma, S, bw, kernel=dnorm) {
  for(i in seq_len(S)) {
    z <- step_optimization(
      z, obs, gamma,
      fit_deriv = fit_lc, fit_deriv_opts = list(bw = bw, kernel = kernel),
      fit_traj = update_trajectory)
  }
  return(z)
}
