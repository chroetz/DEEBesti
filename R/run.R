applyMethodToModel <- function(
    opts,
    observationPath,
    taskPath = observationPath,
    submissionPath = observationPath,
    hyperParmsList = NULL,
    obsNrFilter = NULL,
    truthNrFilter = NULL
) {

  opts <- asOpts(opts, "Estimation")
  method <- getClassAt(opts$method, 2)
  if (is.null(hyperParmsList)) {
    hyperParmsList <- makeOpts(c(method, "HyperParms"))
  }

  outDir <- file.path(submissionPath, method)
  if (!file.exists(outDir)) dir.create(outDir)

  writeOpts(hyperParmsList, dir = outDir)
  writeOpts(opts, dir = outDir)

  taskMeta <- DEEBpath::getMetaGeneric(taskPath, tagsFilter = "task")
  meta <- DEEBpath::getMetaGeneric(
    observationPath,
    tagsFilter = c("truth", "obs"),
    nrFilters = list(obs = obsNrFilter, truth = truthNrFilter))

  for (i in seq_len(nrow(meta))) {
    obs <- readTrajs(meta$obsPath[i])
    res <- estimateWithHyperparameterSelection(
      obs,
      hyperParmsList,
      opts,
      verbose = TRUE)
    hpPath <- file.path(outDir, sprintf(
        "truth%04dobs%04dhyperParms", meta$truthNr[i], meta$obsNr[i]))
    writeOpts(res$hyperParms, hpPath)

    # TODO: check where it makes sense to set the derivative
    if (!hasDeriv(res$trajs)) res$trajs <- setDeriv(res$trajs)

    for (j in seq_len(nrow(taskMeta))) {
      info <- c(as.list(meta[i,]), as.list(taskMeta[j,]), list(outDir = outDir))
      writeTaskResult(res, opts, info, obs)
    }
  }
}

writeTaskResult <- function(res, opts, info, obs) {
  info$task <- ConfigOpts::readOptsBare(info$taskPath)
  taskClass <- getClassAt(info$task, 2)
  switch(
    taskClass,
    estiObsTrajs = writeTaskResultEstiObsTrajs(res, opts, info, obs),
    newTrajs = writeTaskResultNewTrajs(res, opts, info),
    velocity = writeTaskResultVelocity(res, opts, info),
    stop("Unknown task class ", taskClass)
  )
}

writeTaskResultNewTrajs <- function(res, opts, info) {
  outTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  init <- makeTrajs(
    time = 0,
    trajId = seq_len(nrow(info$task$initialState)),
    state = info$task$initialState)
  result <- solveOde(
    u0 = init,
    fun = buildDerivFun(res$hyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$trajs)
  writeTrajs(
    result,
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

writeTaskResultVelocity <- function(res, opts, info) {
  gridSides <- lapply(seq_along(info$task$gridSteps), \(i) seq(
    info$task$gridRanges[i,1],
    info$task$gridRanges[i,2],
    info$task$gridSteps[i]
    ))
  states <- as.matrix(expand.grid(gridSides))
  derivFun <- buildDerivFun(res$hyperParms$derivFun)
  derivs <- t(apply(states, 1, \(s) derivFun(0, s, res$trajs)[[1]]))
  result <- makeDerivTrajs(state = states, deriv = derivs)
  writeDerivTrajs(
    result,
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

writeTaskResultEstiObsTrajs <- function(res, opts, info, obs) {
  targetTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  init <- estimateInitialStateAndTime(res, opts, info$task$predictionTime[1], info$task$timeStep, obs)
  outTimes <- seq(
    init$time[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  result <- solveOde(
    u0 = init$initial,
    fun = buildDerivFun(res$hyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$trajs)
  writeTrajs(
    interpolateTrajs(result, targetTimes),
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

estimateInitialStateAndTime <- function(res, opts, startTime, timeStep, obs) {
  name <- getClassAt(opts$initialState, 2)
  if (name == "FromTrajs" || (name == "Choose" && startTime == res$trajs$time[1])) {
    return(list(time = startTime, initial = getInitialState(res$trajs, startTime)))
  } else if (name == "FromObs" || (name == "Choose" && startTime > res$trajs$time[1])) {
    startObs <- getClosestInTime(obs, startTime)
    return(list(time = mean(startObs$time), initial = startObs))
  } else {
    stop("Unkown initOpts ", name)
  }
}

