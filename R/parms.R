
getParmsAndIntitialState <- function(obs, hyperParms, opts, memoize = FALSE) {
  parms <- getParms(obs, hyperParms, opts, memoize = FALSE)
  list(
    parms = parms,
    initial = getInitialState(parms))
}

getParms <- function(obs, hyperParms, opts, memoize = FALSE) {
  opts <- asOpts(opts, "Method")
  method <- getClassAt(opts, 2)
  trajs <- switch(
    method,
    "Colloc" = estimateParmsColloc(obs, hyperParms, opts),
    "Altopi" = getAltopiTraj(obs, hyperParms, opts, memoize = memoize),
    "Trivial" = obs,
    "Interp" = interpolate(obs, hyperParms, opts),
    "Const" = makeTrajsStateConst(obs, mean),
    stop("Unknown method ", method)
  )
  # TODO: check where it makes sense to set the derivative
  if (!hasDeriv(trajs)) trajs <- setDeriv(trajs, opts$derivMethod)
  return(trajs)
}


