estimateTrajs <- function(initState, timeRange, parms, hyperParms) {

  hyperParms <- asOpts(hyperParms, "HyperParms")

  name <- getClassAt(hyperParms, 2)
  esti <- switch(
    name,
    Trajs = solveOde(
      u0 = initState,
      fun = buildDerivFun(hyperParms$derivFun),
      timeRange = timeRange,
      opts = hyperParms$odeSolver,
      parms = parms),
    Esn = estimateTrajsEsn(initState, timeRange, parms, hyperParms),
    Transformer = estimateTrajsTransformer(initState, timeRange, parms, hyperParms),
    Direct = estimateTrajsDirect(initState, timeRange, parms, hyperParms),
    stop("Unknown HyperParms subclass"))

  return(esti)
}


estimateTrajsEsn <- function(initState, timeRange, parms, hyperParms) {

  hyperParms <- asOpts(hyperParms, c("Esn", "HyperParms"))

  esti <- mapTrajs2Trajs(initState, \(startTraj) {
    predictEsn(
      parms$esn,
      startTraj$state,
      timeRange = timeRange)
  })

  return(esti)
}



estimateTrajsTransformer <- function(initState, timeRange, parms, hyperParms) {

  hyperParms <- asOpts(hyperParms, c("Transformer", "HyperParms"))

  esti <- mapTrajs2Trajs(initState, \(startTraj) {
    predictTransformer(
      parms$transformer,
      startTraj$state,
      timeRange = timeRange)
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

