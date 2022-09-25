estimateParmsColloc <- function(obs, hyperParms, opts) {
  opts <- asOpts(opts, c("Colloc", "Method"))
  switch(
    hyperParms$smoothingMethod,
    predict = estimateParmsCollocPredict(obs, hyperParms, opts),
    fit = estimateParmsCollocFit(obs, hyperParms, opts),
    stop("Unknown smoothing method: ", hyperParms$smoothingMethod))
}

estimateParmsCollocFit <- function(obs, hyperParms, opts) {
  fitter <- buildFitter(hyperParms$fitter)
  mapTrajs2Trajs(obs, function(trj) {
    smoothed <- makeTrajs(
      time = trj$time,
      state = fitter(trj$time, trj$state))
    setDeriv(smoothed, hyperParms$derivFitMethod)
  })
}

estimateParmsCollocPredict <- function(obs, hyperParms, opts) {
  stop("estimateParmsCollocPredict() is not working! need sub method of Colloc?")
  mapTrajs2Trajs(obs, function(trj) {
    time <- seq(
      min(obs$time),
      max(obs$time),
      length.out = length(obs$time) * opts$interSteps)
    smoothed <- makeTrajs(
      time = time,
      state = sapply(
        seq_len(getDim(trj)),
        \(j) predictLoess(trj$time, trj$state[, j], time, hyperParms$fitterSpan, hyperParms$fitterDegree)))
    setDeriv(smoothed, hyperParms$derivFitMethod)
  })
}

