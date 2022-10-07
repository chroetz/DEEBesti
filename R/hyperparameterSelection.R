selectHyperparams <- function(obs, hyperParmsListOpts, methodOpts, opts) {
  hyperParmsListOpts <- asOpts(hyperParmsListOpts, c("HyperParms", "List"))
  hyperParmsList <- hyperParmsListOpts$list
  len <- length(hyperParmsList)
  if (len == 0) return(NULL)
  if (len == 1) return(hyperParmsList[[1]])
  pt <- proc.time()
  validationFoldErrors <- vapply(
    seq_len(opts$folds),
    function(k) {
      splitedObs <- splitCrossValidation(obs, fold = k, maxFolds = opts$folds)
      validateHyperparams(
        splitedObs$train,
        splitedObs$vali,
        hyperParmsList,
        methodOpts = methodOpts,
        opts = opts)
    },
    FUN.VALUE = double(len))
  validationFoldErrors <- matrix(validationFoldErrors, nrow = len)
  validationError <- rowMeans(validationFoldErrors)
  message(
    as.vector((proc.time()-pt)["elapsed"]), "s. ",
    "err:", min(validationError))
  minRowIdx <- which.min(validationError)
  return(hyperParmsList[[minRowIdx]])
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
    hyperParmsListOpts,
    opts,
    verbose = FALSE
  ) {
  normalization <- calculateNormalization(obs)
  obsNormed <- normalization$normalize(obs)
  if (
    !inheritsOptsClass(hyperParmsListOpts, "List") &&
    inheritsOptsClass(hyperParmsListOpts, "HyperParms")
  ) {
    optiHyperParms <- asOpts(hyperParmsListOpts, "HyperParms")
  } else {
    optiHyperParms <- selectHyperparams(
      obsNormed, hyperParmsListOpts, opts$method, opts$hyperParmsSelection)
  }
  if (verbose) printHyperParms(optiHyperParms)
  res <- getParmsAndIntitialState(obsNormed, optiHyperParms, opts$method)
  outTimes <- seq(0, max(obs$time), length.out = opts$odeSteps)
  trajFinal <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(optiHyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$parms)
  trajDenormed <- normalization$denormalize(trajFinal)
  return(c(list(trajs = trajDenormed, hyperParms = optiHyperParms)))
}

printHyperParms <- function(hyperParms) {
  method <- getClassAt(hyperParms, 2)
  # TODO write print / format method for Opts
  cat(
    "[", method, "] ",
    paste0(
      names(hyperParms), ": ",
      sapply(
        hyperParms,
        \(x) if (isOpts(x)) paste(oldClass(x)[-length(oldClass(x))], collapse="_") else format(x)),
      collapse=", "),
    "\n",
    sep="")
}
