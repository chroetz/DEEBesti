validateHyperparams <- function(
    obsTrain,
    obsVali,
    hyperParmsSet,
    method,
    odeSteps,
    odeSolverOpts
  ) {
  prepareMemory(method, nrow(hyperParmsSet))
  vali <- sapply(
    seq_len(nrow(hyperParmsSet)),
    \(i) validate(
      obsTrain,
      obsVali,
      hyperParmsSet[i, ],
      method,
      memoize = TRUE,
      odeSteps = odeSteps,
      odeSolverOpts = odeSolverOpts))
  vali[is.na(vali)] <- Inf
  dplyr::bind_cols(hyperParmsSet, validationErr = vali)
}


selectHyperparams <- function(obs, hyperParmsSet, method, opts) {
  # TODO: implement cross validation
  splitedObs <- splitIntoTrainAndValidation(obs, ratio = opts$ratio)
  pt <- proc.time()
  validatedHyperParms <- validateHyperparams(
    splitedObs$train,
    splitedObs$vali,
    hyperParmsSet,
    method = method,
    odeSteps = opts$odeSteps,
    odeSolverOpts = opts$odeSolver)
  message(as.vector((proc.time()-pt)["elapsed"]), "s")
  minRowIdx <- which.min(validatedHyperParms$validationErr)
  return(validatedHyperParms[minRowIdx,])
}

splitIntoTrainAndValidation <- function(trajs, ratio) {
  n <- getCount(trajs)
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
    opts
  ) {
  optiHyperParms <- selectHyperparams(obs, hyperParmsSet, method, opts$crossValidation)
  res <- getParmsAndIntitialState(obs, optiHyperParms, method)
  trajFinal <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(optiHyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$parms)
  return(trajFinal)
}
