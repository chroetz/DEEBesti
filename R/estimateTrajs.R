estimateTrajs <- function(initState, timeRange, parms, hyperParms, opts) {

  hyperParms <- asOpts(hyperParms, "HyperParms")

  name <- getClassAt(hyperParms, 2)
  esti <- switch(
      name,
      Trajs = solveOde(
        u0 = initState,
        fun = buildDerivFun(hyperParms$derivFun),
        timeRange = timeRange,
        opts = opts$odeSolver,
        parms = parms),
      Esn = estimateTrajsEsn(initState, timeRange, parms, hyperParms, opts))

  return(esti)
}


estimateTrajsEsn <- function(initState, timeRange, parms, hyperParms, opts) {

  hyperParms <- asOpts(hyperParms, c("Esn", "HyperParms"))

  startIdx <- DEEButil::whichMinDist(
    target = parms$esn$obs$state,
    query = initState$state)

  timeStep <- getTimeStepTrajs(parms$esn$obs)
  time <- seq(timeRange[1], timeRange[2], timeStep)

  estiState <- predictEsn(
    parms$esn,
    parms$esn$reservoirSeries[startIdx, ],
    length(time))

  esti <- makeTrajs(
    time = time,
    state = estiState)

  return(esti)
}
