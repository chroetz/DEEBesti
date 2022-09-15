solveOde <- function(u0, f, delta_t, t_max, parms) {
  times <- seq(0, t_max, delta_t)
  res <- deSolve::ode(u0, times=times, func=f, parms=parms, method="rk4")
  return(tibble::tibble(time = res[,1], state = res[,-1]))
}
