
getParmsAndIntitialState <- function(obs, hyperParms, opts, memoize = FALSE) {
  opts <- asOpts(opts, "Method")
  method <- getClassAt(opts, 2)
  if (method == "Colloc") {
    trajs <- estimateParmsColloc(obs, hyperParms, opts)
  } else if (method == "Altopi") {
    trajs <- getAltopiTraj(obs, hyperParms, opts, memoize = memoize)
  } else if (method == "Trivial") {
    trajs <- setDeriv(obs)
  } else if (method == "Const") {
    trajs <- makeTrajsStateConst(obs, mean)
  } else {
    stop("Unknown method ", method)
  }
  initialState <- getInitialState(trajs)
  list(
    parms = trajs,
    initialState = initialState)
}

