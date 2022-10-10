
getParmsAndIntitialState <- function(obs, hyperParms, opts, memoize = FALSE) {
  parms <- getParms(obs, hyperParms, opts, memoize = FALSE)
  list(
    parms = parms,
    initial = getInitialState(parms))
}

getParms <- function(obs, hyperParms, opts, memoize = FALSE) {
  opts <- asOpts(opts, "Method")
  method <- getClassAt(opts, 2)
  if (method == "Colloc") {
    trajs <- estimateParmsColloc(obs, hyperParms, opts)
  } else if (method == "Altopi") {
    trajs <- getAltopiTraj(obs, hyperParms, opts, memoize = memoize)
  } else if (method == "Trivial") {
    trajs <- obs
  } else if (method == "Const") {
    trajs <- makeTrajsStateConst(obs, mean)
  } else {
    stop("Unknown method ", method)
  }
  # TODO: check where it makes sense to set the derivative
  if (!hasDeriv(trajs)) trajs <- setDeriv(trajs, opts$derivMethod)
  return(trajs)
}


