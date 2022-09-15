#' @export
step_optimization <- function(z, obs, gamma, fit_deriv, fit_deriv_opts=NULL, fit_traj, fit_traj_opts=NULL) {
  z_new <- z

  x <- z$state[1:(z$steps-1),]
  y <- (z$state[2:z$steps,] - z$state[1:(z$steps-1),]) / z$step_size
  z_new$deriv <- do.call(fit_deriv, c(list(x=x, y=y), fit_deriv_opts))

  z_new$state <- do.call(fit_traj, c(list(z=z_new, obs=obs, gamma=gamma), fit_traj_opts))

  return(z_new)
}
