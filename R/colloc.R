estimateParmsColloc <- function(obs, hyperPrams) {
  fitter <- getFitter(hyperPrams$fitter)
  kernel <- getKernel(hyperPrams$fitterKernel)
  mapTrajs2Trajs(obs, function(trj) {
    smoothed <- makeTrajs(
      time = trj$time,
      state = fitter(
        trj$time,
        trj$state,
        bw = hyperPrams$fitterBw,
        kernel = kernel))
    setDeriv(smoothed)
  })
}

