interpolate <- function(obs, hyperParms, opts) {
  opts <- asOpts(opts, c("Interp", "Method"))
  hyperParms <- asOpts(hyperParms, c("Interp", "HyperParms"))
  mapTrajs2Trajs(obs, .interpolate, hyperParms, opts)
}

.interpolate <- function(obs, hyperParms, opts) {
  d <- getDim(obs)
  interSteps <- opts$interSteps
  time <- c( # TODO: code duplication with DEEBtrajs
    vapply(seq_len(length(obs$time)-1), \(i) {
      seq(obs$time[i], obs$time[i+1], length.out = interSteps+1)[-(interSteps+1)]
    }, double(interSteps)),
    obs$time[length(obs$time)]
  )
  stateAndDeriv <- lapply(seq_len(d), \(j) interpolate1D(obs$time, obs$state[,j], time, opts$target))
  state <- sapply(stateAndDeriv, \(x) x$state)
  deriv <- sapply(stateAndDeriv, \(x) x$deriv)
  makeTrajs(
    time = time,
    state = state,
    deriv = if (is.matrix(deriv)) deriv else NULL,
    trajId = obs$trajId[1])
}

interpolate1D <- function(x, y, xout, target) {
  deriv = NULL
  if (target == "point") {
    sf <- splinefun(x, y)
    state <- sf(xout, 0)
    deriv <- sf(xout, 1)
  } else if (target == "line") {
    state <- spline(x, y, xout = xout)$y
  } else {
    stop("Unknown target ", target)
  }
  list(state = state, deriv = deriv)
}

