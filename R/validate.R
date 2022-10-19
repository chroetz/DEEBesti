validate <- function(
    obsTrain,
    obsVali,
    hyperParms,
    memoize,
    opts
  ) {
  res <- getParmsAndIntitialState(
    obsTrain,
    hyperParms,
    memoize = memoize)
  esti <- solveOde(
    u0 = res$initial,
    fun = buildDerivFun(hyperParms$derivFun),
    times = seq(0, max(obsVali$time), length.out = opts$odeSteps),
    opts = opts$odeSolver,
    parms = res$parms)
  cleanUpParms(res$parms)
  err <- l2err(esti, obsVali) # TODO: make error type an option
  return(err)
}


validateHyperparams <- function(
    obsTrain,
    obsVali,
    hyperParmsSet,
    opts
  ) {
  prepareMemory()
  vali <- sapply(
    seq_len(length(hyperParmsSet)),
    \(i) validate(
      obsTrain,
      obsVali,
      hyperParmsSet[[i]],
      memoize = TRUE,
      opts))
  vali[is.na(vali)] <- Inf
  return(vali)
}


