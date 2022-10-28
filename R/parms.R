getParms <- function(obs, hyperParms, memoize = FALSE) {

  hyperParms <- asOpts(hyperParms, "HyperParms")
  parms <- list()

  # Preprocessing: Normalization.
  parms$normalization <- calculateNormalization(obs, hyperParms$normalize)
  obs <- parms$normalization$normalize(obs)

  # Fit the trajectory.
  method <- getClassAt(hyperParms, 2)
  trajs <- switch(
    method,
    "Colloc" = estimateParmsColloc(obs, hyperParms),
    "Altopi" = getAltopiTraj(obs, hyperParms, memoize = memoize),
    "Trivial" = obs,
    "Interp" = interpolate(obs, hyperParms),
    "Const" = makeTrajsStateConst(obs, mean),
    stop("Unknown method ", method)
  )
  if (!hasDeriv(trajs) && "derivMethod" %in% names(hyperParms)) {
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

