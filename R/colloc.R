estimateParmsColloc <- function(obs, hyperParms) {
  hyperParms <- asOpts(hyperParms, c("Colloc", "HyperParms"))
  switch(
    hyperParms$smoothingMethod,
    predict = estimateParmsCollocPredict(obs, hyperParms),
    fit = estimateParmsCollocFit(obs, hyperParms),
    stop("Unknown smoothing method: ", hyperParms$smoothingMethod))
}

estimateParmsCollocFit <- function(obs, hyperParms) {
  fitter <- buildFitter(hyperParms$fitter)
  mapTrajs2Trajs(obs, function(trj) {
    smoothed <- makeTrajs(
      time = trj$time,
      state = fitter(trj$time, trj$state))
    setDeriv(smoothed, hyperParms$derivFitMethod)
  })
}

estimateParmsCollocPredict <- function(obs, hyperParms) {
  stop("estimateParmsCollocPredict() is not working! need sub method of Colloc?")
  mapTrajs2Trajs(obs, function(trj) {
    time <- seq(
      min(obs$time),
      max(obs$time),
      length.out = length(obs$time) * hyperParms$interSteps)
    smoothed <- makeTrajs(
      time = time,
      state = sapply(
        seq_len(getDim(trj)),
        \(j) predictLoess(trj$time, trj$state[, j], time, hyperParms$fitterSpan, hyperParms$fitterDegree)))
    setDeriv(smoothed, hyperParms$derivFitMethod)
  })
}

