
getParmsAndIntitialState <- function(obs, hyperParms, memoize = FALSE) {
  parms <- getParms(obs, hyperParms, memoize = FALSE)
  list(
    parms = parms,
    initial = getInitialState(parms))
}

getParms <- function(obs, hyperParms, memoize = FALSE) {
  hyperParms <- asOpts(hyperParms, "HyperParms")
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
  result <- list(trajs = trajs)
  if ("derivFun" %in% names(hyperParms)) {
    derivFunClass <- getClassAt(hyperParms$derivFun, 2)
    if (derivFunClass == "GaussianProcess") {
      result$look <- buildLookUpGrid(trajs, distance=hyperParms$derivFun$maxDist)
    }
  }
  return(result)
}


