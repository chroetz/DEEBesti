
getParmsAndIntitialState <- function(obs, hyperParms, method, memoize = FALSE) {
  if (method == "Colloc") {
    smoothed <- estimateParmsColloc(obs, hyperParms$bwTime, hyperParms$kernelTime)
    parms <- as.list(smoothed)
    parms$bw <- hyperParms$bwState
    parms$kernel <- getKernel(hyperParms$kernelState)
    initialState <- getInitialState(smoothed)
  } else if (method == "Altopi") {
    parms <- getAltopiTraj(obs, hyperParms, memoize = memoize)
    initialState <- getInitialState(parms)
  } else {
    stop("Unknown method ", method)
  }
  return(list(parms = parms, initialState = initialState))
}

