getParms <- function(obs, hyperParms, memoize = FALSE) {

  hyperParms <- asOpts(hyperParms, "HyperParms")
  parms <- list()

  # Preprocessing: Normalization.
  parms$normalization <- calculateNormalization(obs, hyperParms$normalize)
  obs <- parms$normalization$normalize(obs)

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
    stop("Unknown method ", method)
  )
  if (!hasDeriv(trajs)) {
    trajs <- setDeriv(trajs, hyperParms$derivMethod)
  }
  parms$trajs <- trajs

  # Pre-calculations for derivFuns.
  if (hyperParms$derivFun$neighbors >= 1) {
    parms$knnFun <- FastKNN::buildKnnFunction(
      trajs$state,
      hyperParms$derivFun$neighbors)
  }
  if (getClassAt(hyperParms$derivFun, 2) == "GlobalLm") {
    lmFuns <- buildLmFuns(hyperParms$derivFun)
    z <- lmFuns$matrix$transform(trajs$state, trajs$deriv)
    coef <- lapply(seq_len(ncol(z)), \(j) {
      X <- lmFuns$matrix$features(trajs$state, j)
      solve.default(crossprod(X), crossprod(X, z[,j]))
    })
    parms$lmFuns <- lmFuns
    parms$coef <- coef
  }

  return(parms)
}


cleanUpParms <- function(parms) {
  if ("knnFun" %in% names(parms)) {
    FastKNN::deleteQueryFunction(parms$knnFun)
  }
}

