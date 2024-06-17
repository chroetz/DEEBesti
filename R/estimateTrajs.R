estimateTrajs <- function(initState, timeRange, parms, hyperParms) {

  hyperParms <- asOpts(hyperParms, "HyperParms")

  name <- getClassAt(hyperParms, 2)
  esti <- switch(
    name,
    Trajs = estimateTrajsTrajs(initState, timeRange, parms, hyperParms),
    Esn = estimateTrajsPropagator(initState, timeRange, parms, hyperParms),
    Linear = estimateTrajsPropagator(initState, timeRange, parms, hyperParms),
    Transformer = estimateTrajsPropagator(initState, timeRange, parms, hyperParms),
    NeuralOde = estimateTrajsNeuralOde(initState, timeRange, parms, hyperParms),
    Direct = estimateTrajsDirect(initState, timeRange, parms, hyperParms),
    stop("Unknown HyperParms subclass"))

  return(esti)
}


estimateTrajsTrajs <- function(initState, timeRange, parms, hyperParms) {

  hyperParms <- asOpts(hyperParms, c("Trajs", "HyperParms"))

  if (hasValue(hyperParms$nInterTimeStepObs)) {
    if (hyperParms$odeSolver$timeStep != 0) {
      warning(
        "Overwriting odeSolver$timeStep ",
        hyperParms$odeSolver$timeStep,
        " by parms$obsTimeStep / hyperParms$nInterTimeStepObs")
    }
    hyperParms$odeSolver$timeStep <- parms$obsTimeStep / hyperParms$nInterTimeStepObs
  }
  solveOde(
    u0 = initState,
    fun = buildDerivFun(hyperParms$derivFun),
    timeRange = timeRange,
    opts = hyperParms$odeSolver,
    parms = parms)
}


estimateTrajsNeuralOde <- function(initState, timeRange, parms, hyperParms) {

  hyperParms <- asOpts(hyperParms, c("NeuralOde", "HyperParms"))

  if (hasValue(hyperParms$outTimeStep)) {
    timeStep <- hyperParms$outTimeStep
  } else if (hasValue(hyperParms$nOutSteps)) {
    timeStep <- diff(timeRange) / hyperParms$nOutSteps
  } else {
    stop("Specifiy outTimeStep or nOutSteps")
  }

  esti <- mapTrajs2Trajs(initState, \(startTraj) {
    predictNeuralOde(
      parms$neuralOde,
      startTraj$state,
      timeRange = timeRange,
      timeStep = timeStep)
  })

  return(esti)
}


estimateTrajsDirect  <- function(initState, timeRange, parms, hyperParms) {

  hyperParms <- asOpts(hyperParms, c("Direct", "HyperParms"))

  esti <- mapTrajs2Trajs(initState, \(startTraj) {
    extendAnalogue(
      makeTrajs(startTraj$time, matrix(startTraj$state, nrow=1)),
      parms,
      requireTime = max(timeRange))
  })

  return(esti)
}

