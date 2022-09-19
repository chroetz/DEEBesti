
optimizeTrajs <- function(trajs, obs, gamma, S, bw, kernel="Normal") {
  for(i in seq_len(S)) {
    trajs <- stepOptimization(
      trajs, obs, gamma,
      fitDeriv = fitLocalConst, fitDerivOpts = list(bw = bw, kernel = buildKernel(kernel)),
      fitTraj = updateTrajectory)
  }
  return(trajs)
}


estimateParmsColloc <- function(obs, bwTime, kernelTime="Normal") {
  smoothed <- makeTrajs(
    time = obs$time,
    trajId = obs$trajId,
    state = fitLocalLinear(obs$time, obs$state, bw=bwTime, kernel=buildKernel(kernelTime)))
  smoothed <- setDeriv(smoothed)
  return(smoothed)
}

estimateTrajsColloc <- function(smoothed, bwState, kernelState, tMax, tStep, derivFun = "LocalConst") {
  parms <- as.list(smoothed)
  parms$bw <- bwState
  parms$kernel <- buildKernel(kernelState)

  mat <- solveOde(buildDerivFun(derivFun), getInitialState(smoothed), tMax=tMax, tStep=tStep, parms=parms)
  esti <- asTrajs(mat)
  esti <- setTrajId(esti, smoothed$trajId[1])

  return(esti)
}
