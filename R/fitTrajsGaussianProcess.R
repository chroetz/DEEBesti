fitTrajsGaussianProcess <- function(obs, opts) {
  opts <- asOpts(opts, c("GaussianProcess", "FitTrajs"))
  outTime <- interpolateTime(obs$time, opts$interSteps)
  mapTrajs2Trajs(
    obs,
    fitTrajGaussianProcess,
    outTime = outTime,
    opts = opts)
}


fitTrajGaussianProcess <- function(traj, outTime, opts) {
  outState <- fitGaussianProcess(
    traj$time,
    traj$state,
    xout = outTime,
    bandwidth = opts$bandwidth,
    regulation = opts$regulation,
    neighborsHalf = opts$neighborsHalf)
  makeTrajs(
    time = outTime,
    state = outState)
}


fitGaussianProcess <- function(x, y, xout, bandwidth, regulation, neighborsHalf) {
  if (!is.matrix(y)) y <- matrix(y, ncol=1)
  stopifnot(nrow(y) == length(x))

  if (length(x) <=  2 * neighborsHalf || neighborsHalf == 0) {
    return(fitGaussianProcessAll(x, y, xout, bandwidth, regulation))
  }

  # Sliding Window of 2*neighborsHalf; at boundaries, also 2*neighborsHalf points are used.
  yout <- matrix(NA_real_, nrow = length(xout), ncol = ncol(y))
  for(i in (neighborsHalf+1):(length(x) - neighborsHalf)) {
    sel <- xout >= (x[i]+x[i-1])/2 & xout < (x[i]+x[i+1])/2
    yout[sel,] <- fitGaussianProcessAll(
      x[(i-neighborsHalf):(i+neighborsHalf)],
      y[(i-neighborsHalf):(i+neighborsHalf),],
      xout[sel],
      bandwidth=bandwidth, regulation=regulation)
  }
  sel <- xout < (x[neighborsHalf+1]+x[neighborsHalf])/2
  yout[sel,] <- fitGaussianProcessAll(
      x[1:(2*neighborsHalf)],
      y[1:(2*neighborsHalf),],
      xout[sel],
      bandwidth=bandwidth, regulation=regulation)
  sel <- xout >= (x[length(x) - neighborsHalf]+x[length(x) - neighborsHalf + 1])/2
  yout[sel,] <- fitGaussianProcessAll(
      x[(length(x) - (2*neighborsHalf) + 1):length(x)],
      y[(length(x) - (2*neighborsHalf) + 1):length(x),],
      xout[sel],
      bandwidth=bandwidth, regulation=regulation)
  return(yout)
}


fitGaussianProcessAll <- function(x, y, xout, bandwidth, regulation) {
  # TODO: Can I directly get the derivative?
  kernelMatrix <- expKernelMatrix1D(x, bandwidth, regulation)
  kernelVectors <- expKernelVectors1D(x, xout, bandwidth)
  crossprod(kernelVectors, solve.default(kernelMatrix, y))
}

