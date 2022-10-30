solveOde <- function(fun, u0, timeRange, opts, parms = NULL) {
  opts <- asOpts(opts, "OdeSolver")
  if (!isTrajs(u0)) { # time column of u0 is ignored
    if (is.matrix(u0)) {
      n <- nrow(u0)
    } else if (is.vector(u0)) {
      n <- 1
    } else {
      stop("u0 must be Trajs, matrix or vector")
    }
    u0 <- makeTrajs(time = 0, state = u0, trajId = seq_len(n))
  }
  times <- seq(timeRange[1], timeRange[2], by = opts$timeStep)
  mapTrajs2Trajs(u0, function(init) {
    u <- do.call(
      deSolve::ode,
      c(
        list(y = init$state, times = times, func = fun, parms = parms),
        opts$additionalArgs
      )
    )
    colnames(u) <- c("time", paste0("state", seq_len(ncol(u)-1)))
    asTrajs(u)
  })
}
