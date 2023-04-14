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

  esti <- mapTrajs2Trajs(initState, \(startTraj) {
    predictEsn(
      parms$esn,
      startTraj$state,
      timeRange = timeRange)
  })

  return(esti)
}
