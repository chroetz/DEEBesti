fitTrajsLocalPolynomial <- function(obs, opts) {
  opts <- asOpts(opts, c("LocalPolynomial", "FitTrajs"))
  outTime <- interpolateTime(obs$time, opts$interSteps)
  mapTrajs2Trajs(
    obs,
    fitTrajLocalPolynomial,
    outTime = outTime,
    bandwidth = opts$bandwidth,
    kernel = getKernel(opts$kernel))
}


fitTrajLocalPolynomial <- function(traj, outTime, bandwidth, kernel) {
  outState <- fitLocalPolynomial(
    traj$time,
    traj$state,
    xout = outTime,
    bandwidth = bandwidth,
    kernel = kernel)
  makeTrajs(
    time = outTime,
    state = outState)
}


fitLocalPolynomial <- function(x, y, xout, bandwidth, kernel) {
  if (!is.matrix(y)) y <- matrix(y, ncol=1)
  stopifnot(nrow(y) == length(x))

  # TODO
  stop("fitLocalPolynomial() is not implemented yet")
}

