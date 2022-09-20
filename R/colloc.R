estimateParmsColloc <- function(obs, hyperPrams) {
  fitter <- getFitter(hyperPrams$fitter)
  kernel <- getKernel(hyperPrams$fitterKernel)
  smoothed <- makeTrajs(
    time = obs$time,
    trajId = obs$trajId,
    state = fitter(
      obs$time,
      obs$state,
      bw = hyperPrams$fitterBw,
      kernel = kernel))
  smoothed <- setDeriv(smoothed)
  return(smoothed)
}

