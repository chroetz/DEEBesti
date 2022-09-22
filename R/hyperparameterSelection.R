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
  pt <- proc.time()
  validationFoldErrors <- sapply(
    seq_len(opts$folds),
    function(k) {
      splitedObs <- splitCrossValidation(obs, fold = k, maxFolds = opts$folds)
      validateHyperparams(
        splitedObs$train,
        splitedObs$vali,
        hyperParmsSet,
        method = method,
        odeSteps = opts$odeSteps,
        odeSolverOpts = opts$odeSolver)
    })
  validationError <- rowMeans(validationFoldErrors)
  message(
    as.vector((proc.time()-pt)["elapsed"]), "s. ",
    "err:", min(validationError))
  minRowIdx <- which.min(validationError)
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

splitCrossValidation <- function(trajs, fold, maxFolds) {
  n <- sum(getCount(trajs))
  iVali <- 1:(n/maxFolds) * maxFolds - (maxFolds - fold)
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
  if (verbose) printHyperParms(optiHyperParms, method)
  res <- getParmsAndIntitialState(obs, optiHyperParms, method)
  trajFinal <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(optiHyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$parms)
  return(list(trajs = trajFinal, hyperParms = optiHyperParms))
}

printHyperParms <- function(hyperParms, method) {
  cat(
    "[", method, "] ",
    paste0(names(hyperParms), ": ", unlist(hyperParms), collapse=", "),
    "\n",
    sep="")
}
