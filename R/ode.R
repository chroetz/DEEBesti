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
  if (is.null(opts$timeStep)) {
    timeStep <- parms$obsTimeStep / opts$nInterTimeStepObs
  } else {
    timeStep <- opts$timeStep
  }
  times <- seq(
    timeRange[1],
    timeRange[2]+2*timeStep,
    by = timeStep)
  times <- times[1:which(times >= timeRange[2])[1]]
  mapTrajs2Trajs(u0, function(init) {
    u <- do.call(
      deSolve::ode,
      c(
        list(
          y = init$state,
          times = times,
          func = fun,
          parms = parms,
          method = opts$method),
        opts$additionalArgs
      )
    )
    colnames(u) <- c("time", paste0("state", seq_len(ncol(u)-1)))
    asTrajs(u)
  })
}
