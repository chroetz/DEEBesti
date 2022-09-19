solveOde <- function(fun, u0, tMax, tStep, opts = list(method = "rk4"), parms = NULL) {
  tm <- seq(0, tMax, by = tStep)
  suppressWarnings(suppressMessages(utils::capture.output( # make silent
    u <- do.call(
      deSolve::ode,
      c(list(y = u0, times = tm, func = fun, parms = parms), opts))
  )))
  colnames(u) <- c("time", paste0("state", seq_len(ncol(u)-1)))
  return(asTrajs(u))
}
