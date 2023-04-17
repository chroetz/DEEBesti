applyDenoisers <- function(trajs, optsState, optsDeriv) {
  mapTrajs2Trajs(
    trajs,
    applyDenoisersOne,
    optsState=optsState,
    optsDeriv=optsDeriv)
}


applyDenoisersOne <- function(traj, optsState, optsDeriv) {

  name <- getClassAt(optsState, 2)
  traj$state <- switch(
    name,
    TotalVariation = denoiseTotalVariation(traj$state, optsState),
    None = traj$state)

  if ("deriv" %in% names(traj)) {
    name <- getClassAt(optsDeriv, 2)
    traj$deriv <- switch(
      name,
      TotalVariation = denoiseTotalVariation(traj$deriv, optsDeriv),
      None = traj$deriv)
  }

  return(traj)
}


denoiseTotalVariation <- function(mat, opts) {
  opts <- asOpts(opts, c("TotalVariation", "Denoiser"))
  if (opts$lambda <= 0) return(mat)
  apply(
    mat,
    2,
    tvR::denoise1,
    lambda = opts$lambda,
    niter = opts$niter,
    method = opts$method)
}
