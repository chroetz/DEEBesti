validateHyperparams <- function(obsTrain, obsVali, hyperParmsSet, method) {
  prepareMemory(method, nrow(hyperParmsSet))
  vali <- sapply(
    seq_len(nrow(hyperParmsSet)),
    \(i) validate(obsTrain, obsVali, hyperParmsSet[i, ], method, memoize = TRUE))
  vali[is.na(vali)] <- Inf
  dplyr::bind_cols(hyperParmsSet, validationErr = vali)
}


selectHyperparams <- function(obs, hyperParmsSet, method, opts) {
  splitedObs <- splitIntoTrainAndValidation(obs, ratio = opts$ratio)

  pt <- proc.time()
  validatedHyperParms <- validateHyperparams(
    splitedObs$train,
    splitedObs$vali,
    hyperParmsSet,
    method = method)
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
  optiHyperParms <- selectHyperparams(obs, hyperParmsSet, method, opts$hyper)
  res <- getParmsAndIntitialState(obs, optiHyperParms, method)
  trajFinal <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(optiHyperParms$derivFun),
    times = outTimes,
    parms = res$parms)
  return(trajFinal)
}
