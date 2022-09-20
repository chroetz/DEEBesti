validateHyperparams <- function(
    obsTrain,
    obsVali,
    hyperParmsSet,
    method,
    odeSteps,
    odeSolverOpts
  ) {
  prepareMemory(method, length(hyperParmsSet))
  vali <- sapply(
    seq_len(length(hyperParmsSet)),
    \(i) validate(
      obsTrain,
      obsVali,
      hyperParmsSet[[i]],
      method,
      memoize = TRUE,
      odeSteps = odeSteps,
      odeSolverOpts = odeSolverOpts))
  vali[is.na(vali)] <- Inf
  return(vali)
}


selectHyperparams <- function(obs, hyperParmsSet, method, opts) {
  # TODO: implement cross validation
  splitedObs <- splitIntoTrainAndValidation(obs, ratio = opts$ratio)
  pt <- proc.time()
  validationErr <- validateHyperparams(
    splitedObs$train,
    splitedObs$vali,
    hyperParmsSet,
    method = method,
    odeSteps = opts$odeSteps,
    odeSolverOpts = opts$odeSolver)
  message(
    as.vector((proc.time()-pt)["elapsed"]), "s. ",
    "err:", min(validationErr))
  minRowIdx <- which.min(validationErr)
  return(hyperParmsSet[[minRowIdx]])
}

splitIntoTrainAndValidation <- function(trajs, ratio) {
  n <- sum(getCount(trajs))
  iVali <- floor(1:(n*ratio) / ratio)
  iTrain <- setdiff(1:n, iVali)
  vali <- trajs[iVali,]
  train <- trajs[iTrain,]
  return(list(train = train, vali = vali))
}


#' @export
estimateWithHyperparameterSelection <- function(
    obs,
    hyperParmsSet,
    method,
    outTimes,
    opts,
    verbose = FALSE
  ) {
  optiHyperParms <- selectHyperparams(obs, hyperParmsSet, method, opts$crossValidation)
  if (verbose) {
    print(optiHyperParms)
  }
  res <- getParmsAndIntitialState(obs, optiHyperParms, method)
  trajFinal <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(optiHyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$parms)
  return(trajFinal)
}
