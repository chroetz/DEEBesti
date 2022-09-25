solveOde <- function(fun, u0, times, opts, parms = NULL) {
  opts <- asOpts(opts, "OdeSolver")
  if (is.null(nrow(u0))) {
    u0 <- matrix(u0, nrow=1)
  }
  trajIds <- rownames(u0)
  trajList <- lapply(seq_len(nrow(u0)), \(i) {
    suppressWarnings(suppressMessages(utils::capture.output( # make silent
      u <- do.call(
        deSolve::ode,
        c(list(y = u0[i, ], times = times, func = fun, parms = parms), opts))
    )))
    colnames(u) <- c("time", paste0("state", seq_len(ncol(u)-1)))
    trajs <- asTrajs(u)
    if (length(trajIds) > 0)
      trajs <- setTrajId(trajs, trajIds[i])
    return(trajs)
  })

  return(bindTrajs(trajList))
}
