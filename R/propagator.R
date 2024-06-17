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
    nextState <- prevState + timeStep * prediction
  } else if (targetType == "state") {
    nextState <- prediction
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


estimateTrajsPropagator <- function(initState, timeRange, parms, hyperParms) {

  esti <- mapTrajs2Trajs(initState, \(startTraj) {
    predictPropagator(
      parms$propagator,
      startTraj$state,
      timeRange = timeRange)
  })

  return(esti)
}


predictPropagatorDeriv <- function(propagator, states, derivOrder) {
  targetType <- propagator$targetType
  if (targetType == "deriv" || targetType == "state") { # TODO: for deriv unnecessary computations but valid result
    deriv <- t(apply(states, 1, \(s) {
      predictedStates <- predictPropagator(propagator, s, len = derivOrder)$state
      polyInterpCoeffs <- polynomialInterpolation(propagator$timeStep * 0:derivOrder, predictedStates)
      polyInterpCoeffs[2,] # derivative at 0 of polynomial is linear coefficient (second coeff)
    }))
  } else {
    stop("Unknown targetType", targetType)
  }
  return(deriv)

}




predictPropagator <- function(propagator, startState, len = NULL, startTime = 0, timeRange = NULL) {

  stopifnot(hasValue(propagator$propagatorType))

  switch(
    propagator$propagatorType,
    Esn = predictEsn(propagator, startState, len, startTime, timeRange),
    Linear = predictLinear(propagator, startState, len, startTime, timeRange),
    Tansformer = predictTransformer(propagator, startState, len, startTime, timeRange),
    stop("Unknown propagator type", propagator$propagatorType))

}
