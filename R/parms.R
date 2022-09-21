
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
  derivFunParms <- selectDerivFunHyperParms(hyperParms)
  list(
    parms = list(
      trajs = trajs,
      derivFun = derivFunParms),
    initialState = initialState)
}

