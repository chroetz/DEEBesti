getParms <- function(obs, hyperParms, memoize = FALSE) {

  hyperParms <- asOpts(hyperParms, "HyperParms")

  normalization <- calculateNormalization(obs, hyperParms$normalize)
  obs <- normalization$normalize(obs)

  name <- getClassAt(hyperParms, 2)
  parms <- c(
    list(normalization = normalization),
    switch(
      name,
      Trajs = getParmsTrajs(obs, hyperParms, memoize),
      Esn = getParmsEsn(obs, hyperParms, memoize)))

  return(parms)
}


getParmsEsn <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Esn", "HyperParms"))

  esn <- createEsn(
    size = hyperParms$size,
    inDim = getDim(obs),
    degree = hyperParms$degree,
    spectralRadius = hyperParms$spectralRadius,
    inWeightScale = hyperParms$inWeightScale,
    bias = hyperParms$bias)

  esn <- trainEsn(
    esn, obs,
    l2Penalty = hyperParms$l2Penalty,
    warmUpLen = hyperParms$warmUpLen,
    initReservoirScale = hyperParms$initReservoirScale)

  return(list(esn = esn))
}


getParmsTrajs <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Trajs", "HyperParms"))

  parms <- list()

  # Fit the trajectory.
  method <- getClassAt(hyperParms$fitTrajs, 2)
  trajs <- switch(
    method,
    "Identity" = obs,
    "Const" = makeTrajsStateConst(obs, mean),
    "InterpolationSpline" = fitTrajsInterpolationSpline(obs, hyperParms$fitTraj),
    "GaussianProcess" = fitTrajsGaussianProcess(obs, hyperParms$fitTraj),
    "LocalPolynomial" = fitTrajsLocalPolynomial(obs, hyperParms$fitTraj),
    "AltOpt" = fitTrajsAltOpt(obs, hyperParms$fitTraj, memoize = memoize),
    "TrajOpt" = fitTrajsTrajOpt(obs, hyperParms$fitTraj, memoize = memoize),
    stop("Unknown method ", method)
  )
  if (is.null(trajs)) return(NULL)
  if (!hasDeriv(trajs)) { # TODO: check meaningful setting of derivMethod...
    trajs <- setDeriv(trajs, hyperParms$derivMethod)
  }
  trajs <- applyDenoisers(trajs, hyperParms$denoiserState, hyperParms$denoiserDeriv)
  parms$trajs <- trajs

  # Pre-calculations for derivFuns.
  if (hyperParms$derivFun$neighbors >= 1) {
    parms$knnFun <- FastKNN::buildKnnFunction(
      trajs$state,
      hyperParms$derivFun$neighbors)
  }
  derivFunName <- getClassAt(hyperParms$derivFun, 2)
  parms <- switch(
    derivFunName,
    GlobalLm = prepareParmsGlobalLm(parms, hyperParms$derivFun),
    Glmnet = prepareParmsGlmnet(parms, hyperParms$derivFun),
    ThreshLm = prepareParmsThreshLm(parms, hyperParms$derivFun),
    parms
  )

  # Pre-calculations for modifiers.
  for (modifierOpts in hyperParms$derivFun$modifierList$list) {
    name <- getClassAt(modifierOpts, 2)
    if (name == "Attraction") {
      parms$attactionKnnFun <- FastKNN::buildKnnFunction(
        trajs$state,
        modifierOpts$neighbors)
    }
  }

  return(parms)
}


cleanUpParms <- function(parms) {
  if ("knnFun" %in% names(parms)) {
    FastKNN::deleteQueryFunction(parms$knnFun)
  }
  if ("attactionKnnFun" %in% names(parms)) {
    FastKNN::deleteQueryFunction(parms$attactionKnnFun)
  }
}

