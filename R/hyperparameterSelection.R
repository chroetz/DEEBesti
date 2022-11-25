selectHyperparams <- function(obs, hyperParmsListOpts, opts) {
  opts <- asOpts(opts, "HyperParmsSelection")
  name <- getClassAt(opts, 2)
  switch(
    name,
    None = selectHyperparamsNone(hyperParmsListOpts),
    CrossValidation = selectHyperparamsCrossValidation(obs, hyperParmsListOpts, opts),
    EndValidation = selectHyperparamsEndValidation(obs, hyperParmsListOpts, opts))
}


selectHyperparamsNone <- function(hyperParmsListOpts) {
  if (
    !inheritsOptsClass(hyperParmsListOpts, "List") &&
    inheritsOptsClass(hyperParmsListOpts, "HyperParms")
  ) {
    return(list(hyperParmshyperParmsListOpts))
  } else {
    hyperParmsList <- hyperParmsListOpts$list
    len <- length(hyperParmsList)
    if (len == 1) return(list(hyperParmshyperParmsList[[1]]))
    stop("There are ", len, " hyperparms to select from but selection method is `None`.")
  }
}


selectHyperparamsCrossValidation <- function(obs, hyperParmsListOpts, opts) {
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
        opts = opts)
    },
    FUN.VALUE = double(len))
  validationFoldErrors <- matrix(validationFoldErrors, nrow = len)
  validationErrors <- rowMeans(validationFoldErrors)
  message(
    as.vector((proc.time()-pt)["elapsed"]), "s. ",
    "err:", min(validationErrors))
  minRowIdx <- which.min(validationErrors)
  return(list(
    hyperParmshyperParmsList[[minRowIdx]],
    validationErrors = validationErrors))
}

splitCrossValidation <- function(trajs, fold, maxFolds) { # TODO: options
  n <- sum(getCount(trajs))
  iVali <- 1:(n/maxFolds) * maxFolds - (maxFolds - fold)
  iTrain <- setdiff(1:n, iVali)
  vali <- trajs[iVali,]
  train <- trajs[iTrain,]
  return(list(train = train, vali = vali))
}


selectHyperparamsEndValidation <- function(obs, hyperParmsListOpts, opts) {
  hyperParmsListOpts <- asOpts(hyperParmsListOpts, c("HyperParms", "List"))
  hyperParmsList <- hyperParmsListOpts$list
  len <- length(hyperParmsList)
  if (len == 0) return(NULL)
  if (len == 1) return(hyperParmsList[[1]])
  pt <- proc.time()
  splitedObs <- splitEndValidation(obs, opts$ratio)
  validationErrors <- validateHyperparams(
    splitedObs$train,
    splitedObs$vali,
    hyperParmsList,
    opts = opts)
  message(
    as.vector((proc.time()-pt)["elapsed"]), "s. ",
    "err:", min(validationErrors))
  minRowIdx <- which.min(validationErrors)
  return(list(
    hyperParmshyperParmsList[[minRowIdx]],
    validationErrors = validationErrors))
}

splitEndValidation <- function(trajs, ratio) {
  train <- mapTrajs2Trajs(trajs, \(traj) traj[1:floor(nrow(traj) * (1-ratio)), ])
  vali <- mapTrajs2Trajs(trajs, \(traj) traj[(floor(nrow(traj) * (1-ratio)) + 1):nrow(traj), ])
  return(list(train = train, vali = vali))
}



#' @export
estimateWithHyperparameterSelection <- function(
    obs,
    hyperParmsListOpts,
    opts,
    verbose = FALSE
  ) {
  selection <- selectHyperparams(
      obs,
      hyperParmsListOpts,
      opts$hyperParmsSelection)
  if (verbose) printHyperParms(selection)
  parms <- getParms(obs, selection$hyperParms)
  return(c(list(parms = parms), selection))
}

printHyperParms <- function(selection) {
  if ("validationErrors" %in% selection)
    cat(paste0(sprintf("%.1f", selection$validationErrors), collapse=", "), "\n")
}
