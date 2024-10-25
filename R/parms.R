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
      Propagator = getParmsPropagator(obs, hyperParms, memoize),
      NeuralOde = getParmsNeuralOde(obs, hyperParms, memoize),
      Direct = getParmsDirect(obs, hyperParms, memoize),
      Const = getParmsConst(obs, hyperParms, memoize),
      stop("Unknown HyperParms subclass")
    )
  )

  return(parms)
}


getParmsPropagator <- function(obs, hyperParms, memoize) {
  hyperParms <- asOpts(hyperParms, c("Propagator", "HyperParms"))
  name <- getClassAt(hyperParms, 3)
  parms <- switch(
    name,
    Esn = getParmsEsn(obs, hyperParms, memoize),
    Linear = getParmsLinear(obs, hyperParms, memoize),
    Transformer = getParmsTransformer(obs, hyperParms, memoize),
    Regression = getParmsRegression(obs, hyperParms, memoize),
    RecurrentNet = getParmsRecurrentNet(obs, hyperParms, memoize),
    stop("Unknown Propagator subclass")
  )
  return(parms)
}


# Propagator map Estimation: Linear (Next Generation Reservoir Computing)
getParmsLinear <- function(obs, hyperParms, memoize) {
  hyperParms <- asOpts(hyperParms, c("Linear", "Propagator", "HyperParms"))
  parms <- createLinear(obs, hyperParms)
  return(parms)
}


getParmsRegression <- function(obs, hyperParms, memoize) {
  hyperParms <- asOpts(hyperParms, c("Regression", "Propagator", "HyperParms"))
  parms <- createRegression(obs, hyperParms)
  return(parms)
}


# Propagator map Estimation: Echo State Network and Random Features
getParmsEsn <- function(obs, hyperParms, memoize) {
  hyperParms <- asOpts(hyperParms, c("Esn", "Propagator", "HyperParms"))
  inDim <- getDim(obs)
  if (hyperParms$timeStepAsInput) {
    inDim <- inDim + 1
  }
  parms <- createEsn(hyperParms, inDim = inDim)
  parms <- trainEsn(parms, obs, hyperParms)
  return(parms)
}



# Propagator map Estimation: RNN and LSTM
getParmsRecurrentNet <- function(obs, hyperParms, memoize) {
  hyperParms <- asOpts(hyperParms, c("RecurrentNet", "Propagator", "HyperParms"))
  parms <- createRecurrentNet(obs, hyperParms)
  return(parms)
}



# Propagator map Estimation: Transformer
getParmsTransformer <- function(obs, hyperParms, memoize) {
  hyperParms <- asOpts(hyperParms, c("Transformer", "Propagator", "HyperParms"))
  inDim <- getDim(obs)
  if (hyperParms$timeStepAsInput) {
    inDim <- inDim + 1
  }
  parms <- createTransformer(obs, hyperParms)
  return(parms)
}



getParmsNeuralOde <- function(obs, hyperParms, memoize) {
  hyperParms <- asOpts(hyperParms, c("NeuralOde", "HyperParms"))
  neuralOde <- createAndTrainNeuralOde(hyperParms, obs)
  return(list(neuralOde = neuralOde))
}


# No training required
getParmsDirect <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Direct", "HyperParms"))

  obs <-
    obs |>
    dplyr::arrange(.data$trajId, .data$time)
  store <-
    obs |>
    dplyr::group_by(.data$trajId) |>
    tibble::rowid_to_column("obsIdx") |>
    dplyr::filter(
      order(.data$time, decreasing=TRUE) > hyperParms$requiredFutures)

  parms <- list()

  parms$knnIdxToStoreIdx <- store$obsIdx
  parms$stored <- obs
  parms$knnFun <- FastKNN::buildKnnFunction(store$state, 1)

  return(parms)
}


# Constant Output
getParmsConst <- function(obs, hyperParms, memoize) {

  hyperParms <- asOpts(hyperParms, c("Const", "HyperParms"))

  constState <- switch(
    hyperParms$method,
    Mean = colMeans(obs$state),
    Last = obs$state[nrow(obs$state), ],
    stop("Unknown method for extracting constant state")
  )

  parms <- list(state = constState)

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
    "InterpolationSpline" = fitTrajsInterpolationSpline(obs, hyperParms$fitTrajs),
    "GaussianProcess" = fitTrajsGaussianProcess(obs, hyperParms$fitTrajs),
    "LocalPolynomial" = fitTrajsLocalPolynomial(obs, hyperParms$fitTrajs),
    "AltOpt" = fitTrajsAltOpt(obs, hyperParms$fitTrajs, memoize = memoize),
    "TrajOpt" = fitTrajsTrajOpt(obs, hyperParms$fitTrajs, memoize = memoize),
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

