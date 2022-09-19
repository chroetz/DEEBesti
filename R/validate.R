getValidateFun <- function(method) {
  switch(
    method,
    Altopi = validateAltopi,
    \(obsTrain, obsVali, hyperParms) validateDefault(obsTrain, obsVali, hyperParms, method = method))
}

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
    hyperParms
  ) {
  if (nrow(validateAltopiMemory) > 0) {
    selMem <-
      validateAltopiMemory$bw == hyperParms$bw &
      validateAltopiMemory$gamma == hyperParms$gamma &
      validateAltopiMemory$kernel == hyperParms$kernel &
      validateAltopiMemory$derivFun == hyperParms$derivFun &
      validateAltopiMemory$interSteps == hyperParms$interSteps
    selMemThis <- selMem & validateAltopiMemory$S == hyperParms$S
    if (any(selMemThis)) return(validateAltopiMemory$err[which(selMemThis)[1]])
    selMemPrev <- selMem & validateAltopiMemory$S == hyperParms$S-1
  } else {
    selMemPrev <- FALSE
  }
  if (any(selMemPrev)) {
    trajs <- validateAltopiMemory$trajs[[which(selMemPrev)[1]]]
    sStart <- hyperParms$S
  } else {
    trajs <- initAltopi(obsTrain, interSteps = hyperParms$interSteps)
    sStart <- 1
  }
  for(i in sStart:hyperParms$S) {
    trajs <- stepOptimization(
      trajs, obsTrain, hyperParms$gamma,
      fitDeriv = fitLocalConst, fitDerivOpts = list(bw = hyperParms$bw, kernel = buildKernel(hyperParms$kernel)),
      fitTraj = updateTrajectory)
  }
  tMax <- max(obsVali$time)
  esti <- solveOde(
    buildDerivFun(hyperParms$derivFun),
    getInitialState(trajs),
    tStep = tMax / 1e3, tMax = tMax,
    parms = trajs)
  err <- l2err(esti, obsVali)
  validateAltopiMemory <<- dplyr::bind_rows(
    validateAltopiMemory,
    tibble::tibble(
      S = hyperParms$S,
      bw = hyperParms$bw,
      gamma = hyperParms$gamma,
      kernel = hyperParms$kernel,
      derivFun = hyperParms$derivFun,
      interSteps = hyperParms$interSteps,
      err = err,
      trajs = list(trajs)))
  return(err)
}

clearMemory <- function() {
  validateAltopiMemory <<- validateAltopiMemory[0,]
}

validateDefault <- function(obsTrain, obsVali, hyperParms, method) {
  res <- getParmsAndIntitialState(obsTrain, hyperParms, method = method)
  tMax <- max(obsVali$time)
  esti <- solveOde(
    u0 = res$initialState,
    fun = buildDerivFun(hyperParms$derivFun),
    tStep = tMax / 1e3,
    tMax = tMax,
    parms = res$parms)
  err <- l2err(esti, obsVali)
  return(err)
}

