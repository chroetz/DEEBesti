validateAltopiMemory <- tibble::tibble(
  bw = double(0),
  S = integer(0),
  gamma = double(0),
  kernel = character(0),
  derivFun = character(0),
  interSteps = integer(0),
  err = double(0),
  trajs = list()
)

validateAltopi <- function(
    obsTrain,
    obsVali,
    parms
  ) {
  if (nrow(validateAltopiMemory) > 0) {
    selMem <-
      validateAltopiMemory$bw == parms$bw &
      validateAltopiMemory$gamma == parms$gamma &
      validateAltopiMemory$kernel == parms$kernel &
      validateAltopiMemory$derivFun == parms$derivFun &
      validateAltopiMemory$interSteps == parms$interSteps
    selMemThis <- selMem & validateAltopiMemory$S == parms$S
    if (any(selMemThis)) return(validateAltopiMemory$err[which(selMemThis)[1]])
    selMemPrev <- selMem & validateAltopiMemory$S == parms$S-1
  } else {
    selMemPrev <- FALSE
  }
  if (any(selMemPrev)) {
    trajs <- validateAltopiMemory$trajs[[which(selMemPrev)[1]]]
    sStart <- parms$S
  } else {
    trajs <- initAltopi(obsTrain, interSteps = parms$interSteps)
    sStart <- 1
  }
  for(i in sStart:parms$S) {
    trajs <- stepOptimization(
      trajs, obsTrain, parms$gamma,
      fitDeriv = fitLocalConst, fitDerivOpts = list(bw = parms$bw, kernel = buildKernel(parms$kernel)),
      fitTraj = updateTrajectory)
  }
  tMax <- max(obsVali$time)
  esti <- solveOde(
    buildDerivFun(parms$derivFun),
    getInitialState(trajs),
    tStep = tMax / 1e3, tMax = tMax,
    parms = trajs)
  err <- l2err(esti, obsVali)
  validateAltopiMemory <<- dplyr::bind_rows(
    validateAltopiMemory,
    tibble::tibble(
      S = parms$S,
      bw = parms$bw,
      gamma = parms$gamma,
      kernel = parms$kernel,
      derivFun = parms$derivFun,
      interSteps = parms$interSteps,
      err = err,
      trajs = list(trajs)))
  return(err)
}

clearMemory <- function() {
  validateAltopiMemory <<- validateAltopiMemory[0,]
}

validateColloc <- function(obsTrain, obsVali, parms) {
  smoothed <- estimateParmsColloc(obsTrain, parms$bwTime, kernelTime=parms$kernelTime)
  tMax <- max(obsVali$time)
  esti <- estimateTrajsColloc(smoothed, parms$bwState, parms$kernelState, tMax, tMax / 1e3)
  err <- l2err(esti, obsVali)
  return(err)
}
