solveOde <- function(fun, u0, times, opts, parms = NULL) {
  suppressWarnings(suppressMessages(utils::capture.output( # make silent
    u <- do.call(
      deSolve::ode,
      c(list(y = u0, times = times, func = fun, parms = parms), opts))
  )))
  colnames(u) <- c("time", paste0("state", seq_len(ncol(u)-1)))
  return(asTrajs(u))
}
