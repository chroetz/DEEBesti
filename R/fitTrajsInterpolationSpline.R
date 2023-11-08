fitTrajsInterpolationSpline <- function(obs, opts) {
  opts <- asOpts(opts, c("InterpolationSpline", "FitTrajs"))
  outTime <- interpolateTimeAdaptive(obs, opts$interSteps)
  applyToEachDimAndId(obs, fitTrajInterpolationSpline1D, outTime, target=opts$target)
}


fitTrajInterpolationSpline1D <- function(time, state, outTime, target) {
  outDeriv = NULL
  if (target == "point") {
    sf <- stats::splinefun(time, state)
    outState <- sf(outTime, 0)
    outDeriv <- sf(outTime, 1)
  } else if (target == "line") {
    outState <- stats::spline(time, state, xout = outTime)$y
  } else {
    stop("Unknown target ", target)
  }
  list(state = outState, deriv = outDeriv)
}
