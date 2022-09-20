
getParmsAndIntitialState <- function(obs, hyperParms, method, memoize = FALSE) {
  if (method == "Colloc") {
    trajs <- estimateParmsColloc(obs, hyperParms)
  } else if (method == "Altopi") {
    trajs <- getAltopiTraj(obs, hyperParms, memoize = memoize)
  } else if (method == "Trivial") {
    trajs <- setDeriv(obs)
  } else {
    stop("Unknown method ", method)
  }
  initialState <- getInitialState(trajs)
  derivFunParms <- list()
  if (hyperParms$derivFun == "NearestNeighbor") {
  } else if (hyperParms$derivFun == "WeightedKNN") {
    derivFunParms$k <- hyperParms$derivFunK
    derivFunParms$p <- hyperParms$derivFunP
  } else if (hyperParms$derivFun %in% c("LocalConst", "LocalLinear")) {
    derivFunParms$bw <- hyperParms$derivFunBw
    derivFunParms$kernel <- getKernel(hyperParms$derivFunKernel)
  } else {
    stop("Unknown derivFun name ", hyperParms$derivFun)
  }
  list(
    parms = list(
      trajs = trajs,
      derivFun = derivFunParms),
    initialState = initialState)
}

