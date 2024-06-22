getPropagatorRegressionOut <- function(obs, targetType) {
  if (targetType == "deriv") {
    regressionOut <- do.call(
      rbind,
      applyTrajId(
        obs,
        \(traj) {
          (traj$state[-1,, drop=FALSE] - traj$state[-nrow(traj$state),, drop=FALSE]) / diff(traj$time)
        }
      )
    )
  } else if (targetType == "state") {
    regressionOut <- do.call(rbind, applyTrajId(obs, \(traj) traj$state[-1,, drop=FALSE]))
  } else {
    stop("Unknown targetType", targetType)
  }
  return(regressionOut)
}


getPropagatorNextState <- function(prevState, timeStep, prediction, targetType) {
  if (targetType == "deriv") {
    nextState <- as.vector(prevState) + timeStep * as.vector(prediction)
  } else if (targetType == "state") {
    nextState <- as.vector(prediction)
  } else {
    stop("Unknown targetType", targetType)
  }
  return(nextState)
}



polynomialInterpolation <- function(x, y) {
  stopifnot(length(x) == nrow(y))
  p <- length(x) - 1
  X <- outer(x, (0:p), `^`)
  coeff <- DEEButil::saveSolve(crossprod(X), crossprod(X, y))
  return(coeff)
}


estimateTrajsPropagator <- function(initState, timeRange, parms, opts) {

  opts <- asOpts(opts, c("Propagator", "HyperParms"))

  esti <- mapTrajs2Trajs(initState, \(startTraj) {
    predictPropagator(
      parms,
      opts,
      startTraj$state,
      timeRange = timeRange)
  })

  return(esti)
}


predictPropagatorDeriv <- function(parms, opts, states, derivOrder) {
  opts <- asOpts(opts, c("Propagator", "HyperParms"))
  if (opts$targetType == "deriv" || opts$targetType == "state") { # TODO: for deriv unnecessary computations but valid result
    deriv <- t(apply(states, 1, \(s) {
      predictedStates <- predictPropagator(parms, opts, s, len = derivOrder)$state
      polyInterpCoeffs <- polynomialInterpolation(parms$timeStep * 0:derivOrder, predictedStates)
      polyInterpCoeffs[2,] # derivative at 0 of polynomial is linear coefficient (second coeff)
    }))
  } else {
    stop("Unknown targetType", opts$targetType)
  }
  return(deriv)

}


predictPropagator <- function(parms, opts, startState, len = NULL, startTime = 0, timeRange = NULL) {

  opts <- asOpts(opts, c("Propagator", "HyperParms"))

  if (is.null(timeRange)) {
    stopifnot(length(len) == 1, len >= 0)
    time <- startTime + (0:len)*parms$timeStep
  } else {
    stopifnot(length(timeRange) == 2)
    time <- seq(timeRange[1], timeRange[2], by = parms$timeStep)
    if (time[length(time)] < timeRange[2]) {
      time <- c(time, time[length(time)] + parms$timeStep)
    }
    len <- length(time) - 1
  }

  name <- getClassAt(opts, 3)
  predictedStates <- switch(
    name,
    Esn = predictEsn(parms, opts, startState, len),
    Linear = predictLinear(parms, opts, startState, len),
    Transformer = predictTransformer(parms, opts, startState, len),
    Regression = predictRegression(parms, opts, startState, len),
    stop("Unknown propagator type", name))

  makeTrajs(
    time = time,
    state = predictedStates)
}
