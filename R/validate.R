validate <- function(
    obsTrain,
    obsVali,
    hyperParms,
    method,
    memoize,
    odeSteps,
    odeSolverOpts
  ) {
  res <- getParmsAndIntitialState(
    obsTrain,
    hyperParms,
    method,
    memoize = memoize)
  esti <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(hyperParms$derivFun),
    times = seq(0, max(obsVali$time), length.out = odeSteps),
    opts = odeSolverOpts,
    parms = res$parms)
  err <- l2err(esti, obsVali)
  return(err)
}

