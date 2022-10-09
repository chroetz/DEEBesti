validate <- function(
    obsTrain,
    obsVali,
    hyperParms,
    methodOpts,
    memoize,
    opts
  ) {
  res <- getParmsAndIntitialState(
    obsTrain,
    hyperParms,
    methodOpts,
    memoize = memoize)
  esti <- solveOde(
    u0 = res$initial,
    fun = buildDerivFun(hyperParms$derivFun),
    times = seq(0, max(obsVali$time), length.out = opts$odeSteps),
    opts = opts$odeSolver,
    parms = res$parms)
  err <- l2err(esti, obsVali)
  return(err)
}


validateHyperparams <- function(
    obsTrain,
    obsVali,
    hyperParmsSet,
    methodOpts,
    opts
  ) {
  prepareMemory(methodOpts, length(hyperParmsSet))
  vali <- sapply(
    seq_len(length(hyperParmsSet)),
    \(i) validate(
      obsTrain,
      obsVali,
      hyperParmsSet[[i]],
      methodOpts,
      memoize = TRUE,
      opts))
  vali[is.na(vali)] <- Inf
  return(vali)
}


