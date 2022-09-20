getValidateFun <- function(method) {
  switch(
    method,
    Altopi = validateAltopi,
    \(obsTrain, obsVali, hyperParms) validateDefault(obsTrain, obsVali, hyperParms, method = method))
}

validationMemory <- NULL

validateAltopi <- function(
    obsTrain,
    obsVali,
    hyperParms
  ) {
  err <- getErrFromValidationMemory(hyperParms)
  if (!is.null(err)) return(err)
  preHyperParms <- getHyperParmaPredecessorAltopi(hyperParms)
  trajs <- getTrajsFromValidationMemory(preHyperParms)
  if (is.null(trajs)) {
    trajs <- initAltopi(obsTrain, hyperParms)
  }
  newTrajs <- oneAltopiStep(
    trajs,
    obsTrain,
    hyperParms$gamma,
    fitDeriv = fitLocalConst,
    fitDerivOpts = list(bw = hyperParms$bw, kernel = getKernel(hyperParms$kernel)),
    fitTraj = updateAltopiTraj)
  esti <- solveOde(
    buildDerivFun(hyperParms$derivFun),
    getInitialState(newTrajs),
    times = seq(0, max(obsVali$time), length.out = 1e3),
    parms = newTrajs)
  err <- l2err(esti, obsVali)
  addToValidationMemory(hyperParms, err, newTrajs)
  return(err)
}

getHyperParmaPredecessorAltopi <- function(hyperParms) {
  hyperParms$S <- hyperParms$S - 1
  return(hyperParms)
}

addToValidationMemory <- function(hyperParms, err, trajs) {
  utils::sethash(validationMemory, hyperParms, list(err = err, trajs = trajs))
}

getErrFromValidationMemory <- function(hyperParms) {
  utils::gethash(validationMemory, hyperParms, nomatch = NULL)$err
}

getTrajsFromValidationMemory <- function(hyperParms) {
  utils::gethash(validationMemory, hyperParms, nomatch = NULL)$trajs
}

prepareValidationMemory <- function(method, n) {
  validationMemory <<- utils::hashtab("identical", n)
}

validateDefault <- function(obsTrain, obsVali, hyperParms, method) {
  res <- getParmsAndIntitialState(obsTrain, hyperParms, method = method)
  esti <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(hyperParms$derivFun),
    times = seq(0, max(obsVali$time), length.out = 1e3),
    parms = res$parms)
  err <- l2err(esti, obsVali)
  return(err)
}

