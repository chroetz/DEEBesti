validate <- function(obsTrain, obsVali, hyperParms, method, memoize) {
  res <- getParmsAndIntitialState(
    obsTrain,
    hyperParms,
    method,
    memoize = memoize)
  esti <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(hyperParms$derivFun),
    times = seq(0, max(obsVali$time), length.out = 1e3), # TODO put magic number into opts
    parms = res$parms)
  err <- l2err(esti, obsVali)
  return(err)
}

