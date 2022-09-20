estimateParmsColloc <- function(obs, bwTime, kernelTime="Normal") {
  smoothed <- makeTrajs(
    time = obs$time,
    trajId = obs$trajId,
    state = fitLocalLinear(obs$time, obs$state, bw=bwTime, kernel=getKernel(kernelTime)))
  smoothed <- setDeriv(smoothed)
  return(smoothed)
}

