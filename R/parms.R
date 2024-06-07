getParms <- function(obs, hyperParms, memoize = FALSE) {

  hyperParms <- asOpts(hyperParms, "HyperParms")

  normalization <- calculateNormalization(obs, hyperParms$normalize)

  if (!length(hyperParms$obsTimeSpanTarget) == 0) {
    avgTimeSpan <- mean(unlist(applyTrajId(obs, \(traj) diff(range(traj$time)))))
    timeScaling <- hyperParms$obsTimeSpanTarget / avgTimeSpan
  } else {
    timeScaling <- 1
  }

  normalizeParms <- list(
    normalization = normalization,
    timeScaling = timeScaling)

  obs <- normalize(obs, normalizeParms)

  name <- getClassAt(hyperParms, 2)
  parms <- c(
    normalizeParms,
    switch(
      name,
      Trajs = getParmsTrajs(obs, hyperParms, memoize),
      Esn = getParmsEsn(obs, hyperParms, memoize),
      Linear = getParmsLinear(obs, hyperParms, memoize),
      Transformer = getParmsTransformer(obs, hyperParms, memoize),
      Direct = getParmsDirect(obs, hyperParms, memoize),
      stop("Unknown HyperParms subclass")
    )
  )

  return(parms)
}


# Propagator map Estimation: Linear (Next Generation Reservoir Computing)
getParmsLinear <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Linear", "HyperParms"))

  linear <- createLinear(
    obs,
    timeStepAsInput = hyperParms$timeStepAsInput,
    pastSteps = hyperParms$pastSteps,
    skip = hyperParms$skip,
    polyDeg = hyperParms$polyDeg,
    l2Penalty = hyperParms$l2Penalty)

  return(list(linear = linear))
}


# Propagator map Estimation: Echo State Network and Random Features
getParmsEsn <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Esn", "HyperParms"))

  inDim <- getDim(obs)
  if (hyperParms$timeStepAsInput) {
    inDim <- inDim + 1
  }

  esn <- createEsn(
    size = hyperParms$size,
    inDim = inDim,
    degree = hyperParms$degree,
    spectralRadius = hyperParms$spectralRadius,
    inWeightScale = hyperParms$inWeightScale,
    bias = hyperParms$bias,
    seed = hyperParms$seed)

  esn <- trainEsn(
    esn, obs,
    l2Penalty = hyperParms$l2Penalty,
    warmUpLen = hyperParms$warmUpLen,
    initReservoirScale = hyperParms$initReservoirScale,
    timeStepAsInput = hyperParms$timeStepAsInput,
    skip = hyperParms$skip)

  return(list(esn = esn))
}



# Propagator map Estimation: Transformer
getParmsTransformer <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Transformer", "HyperParms"))

  inDim <- getDim(obs)
  if (hyperParms$timeStepAsInput) {
    inDim <- inDim + 1
  }

  transformer <- createTransformer(hyperParms, stateDim = inDim)

  transformer <- trainTransformer(transformer, obs, hyperParms)

  return(list(transformer = transformer))
}


# No training required
getParmsDirect <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Direct", "HyperParms"))

  obs <-
    obs |>
    dplyr::arrange(trajId, time)
  store <-
    obs |>
    dplyr::group_by(trajId) |>
    tibble::rowid_to_column("obsIdx") |>
    dplyr::filter(
      order(time, decreasing=TRUE) > hyperParms$requiredFutures)

  parms <- list()

  parms$knnIdxToStoreIdx <- store$obsIdx
  parms$stored <- obs
  parms$knnFun <- FastKNN::buildKnnFunction(store$state, 1)

  return(parms)
}


# Two-Step Procedure:
# 1. estimate observed trajectory
# 2. estimate model function
getParmsTrajs <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Trajs", "HyperParms"))

  parms <- list()

  parms$obsTimeStep <- getTimeStepTrajs(obs, requireConst=FALSE)

  # Fit the trajectory.
  method <- getClassAt(hyperParms$fitTrajs, 2)
  trajs <- switch(
    method,
    "Identity" = obs,
    "ConstMean" = makeTrajsStateConst(obs, mean),
    "ConstLast" = makeTrajsStateConst(obs, dplyr::last),
    "InterpolationSpline" = fitTrajsInterpolationSpline(obs, hyperParms$fitTraj),
    "GaussianProcess" = fitTrajsGaussianProcess(obs, hyperParms$fitTraj),
    "LocalPolynomial" = fitTrajsLocalPolynomial(obs, hyperParms$fitTraj),
    "AltOpt" = fitTrajsAltOpt(obs, hyperParms$fitTraj, memoize = memoize),
    "TrajOpt" = fitTrajsTrajOpt(obs, hyperParms$fitTraj, memoize = memoize),
    stop("Unknown method ", method)
  )
  if (is.null(trajs)) return(NULL)
  if (nrow(trajs) <= 1) stop("Something is wrong: Fitted Trajs has 0 or 1 rows.")
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

