interpolate <- function(obs, hyperParms) {
  hyperParms <- asOpts(hyperParms, c("Interp", "HyperParms"))
  mapTrajs2Trajs(obs, .interpolate, hyperParms)
}

.interpolate <- function(obs, hyperParms) {
  d <- getDim(obs)
  interSteps <- hyperParms$interSteps
  time <- c( # TODO: code duplication with DEEBtrajs
    vapply(seq_len(length(obs$time)-1), \(i) {
      seq(obs$time[i], obs$time[i+1], length.out = interSteps+1)[-(interSteps+1)]
    }, double(interSteps)),
    obs$time[length(obs$time)]
  )
  stateAndDeriv <- lapply(seq_len(d), \(j) interpolate1D(obs$time, obs$state[,j], time, hyperParms))
  state <- sapply(stateAndDeriv, \(x) x$state)
  deriv <- sapply(stateAndDeriv, \(x) x$deriv)
  makeTrajs(
    time = time,
    state = state,
    # interpolation method may or may not calculate derivative
    deriv = if (is.matrix(deriv)) deriv else NULL,
    trajId = obs$trajId[1])
}

interpolate1D <- function(x, y, xout, opts) {
  deriv = NULL
  # TODO: allow to use interpolate1DViaGpKnn()
  if (otps$target == "point") {
    sf <- stats::splinefun(x, y)
    state <- sf(xout, 0)
    deriv <- sf(xout, 1)
  } else if (otps$target == "line") {
    state <- stats::spline(x, y, xout = xout)$y
  } else {
    stop("Unknown target ", otps$target)
  }
  list(state = state, deriv = deriv)
}

interpolate1DViaGp <- function(x, y, xout, bandwidth, regulation) {
  # TODO: Can I directly get the derivative?
  kernelMatrix <- expKernelMatrix1D(x, bandwidth, regulation)
  kernelVectors <- expKernelVectors1D(x, xout, bandwidth)
  crossprod(kernelVectors, solve.default(kernelMatrix, y))
}

interpolate1DViaGpKnn <- function(x, y, xout, bandwidth, regulation, neighborsHalf) {
  if (length(x) <=  2 * neighborsHalf) {
    return(interpolate1DViaGp(x, y, xout, bandwidth, regulation))
  }
  yout <- double(length(xout))
  for(i in (neighborsHalf+1):(length(x) - neighborsHalf)) {
    sel <- xout >= (x[i]+x[i-1])/2 & xout < (x[i]+x[i+1])/2
    yout[sel] <- interpolate1DViaGp(
      x[(i-neighborsHalf):(i+neighborsHalf)],
      y[(i-neighborsHalf):(i+neighborsHalf)],
      xout[sel],
      bandwidth=bandwidth, regulation=regulation)
  }
  sel <- xout < (x[neighborsHalf+1]+x[neighborsHalf])/2
  yout[sel] <- interpolate1DViaGp(
      x[1:neighborsHalf],
      y[1:neighborsHalf],
      xout[sel],
      bandwidth=bandwidth, regulation=regulation)
  sel <- xout >= (x[length(x) - neighborsHalf]+x[length(x) - neighborsHalf + 1])/2
  yout[sel] <- interpolate1DViaGp(
      x[(length(x) - neighborsHalf + 1):length(x)],
      y[(length(x) - neighborsHalf + 1):length(x)],
      xout[sel],
      bandwidth=bandwidth, regulation=regulation)
  return(yout)
}



