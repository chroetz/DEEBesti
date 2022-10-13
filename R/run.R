applyMethodToModel <- function(
    opts,
    hyperParmsList = NULL,
    observationPath,
    taskPath = observationPath,
    submissionPath = observationPath,
    obsNrFilter = NULL,
    truthNrFilter = NULL
) {

  opts <- asOpts(opts, "Estimation")

  outDir <- file.path(submissionPath, hyperParmsList$name)
  if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)

  writeOpts(hyperParmsList, dir = outDir)
  writeOpts(opts, dir = outDir)

  taskMeta <- DEEBpath::getMetaGeneric(taskPath, tagsFilter = "task")
  meta <- DEEBpath::getMetaGeneric(
    observationPath,
    tagsFilter = c("truth", "obs"),
    nrFilters = list(obsNr = obsNrFilter, truthNr = truthNrFilter))

  for (i in seq_len(nrow(meta))) {
    info <- meta[i,]
    obs <- readTrajs(info$obsPath)
    res <- estimateWithHyperparameterSelection(
      obs,
      hyperParmsList,
      opts,
      verbose = TRUE)
    hpPath <- file.path(outDir, DEEBpath::hyperParmsFile(info))
    writeOpts(res$hyperParms, hpPath)

    for (j in seq_len(nrow(taskMeta))) {
      allInfo <- c(as.list(info), as.list(taskMeta[j,]), list(outDir = outDir))
      writeTaskResult(res, opts, allInfo)
    }
  }
}

writeTaskResult <- function(res, opts, info) {
  info$task <- ConfigOpts::readOptsBare(info$taskPath)
  taskClass <- getClassAt(info$task, 2)
  switch(
    taskClass,
    estiObsTrajs = writeTaskResultEstiObsTrajs(res, opts, info),
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
    parms = res$parms)
  resultDenormed <- res$normalization$denormalize(result)
  writeTrajs(
    resultDenormed,
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

writeTaskResultVelocity <- function(res, opts, info) {
  gridSides <- lapply(seq_along(info$task$gridSteps), \(i) seq(
    info$task$gridRanges[i,1],
    info$task$gridRanges[i,2],
    info$task$gridSteps[i]
    ))
  grid <- makeDerivTrajs(state = as.matrix(expand.grid(gridSides)))
  gridNormed <- res$normalization$normalize(grid)
  derivFun <- buildDerivFun(res$hyperParms$derivFun)
  derivs <- t(apply(gridNormed$state, 1, \(s) derivFun(0, s, res$parms)[[1]]))
  resultNormed <- makeDerivTrajs(state = gridNormed$state, deriv = derivs)
  result <- res$normalization$denormalize(resultNormed)
  writeDerivTrajs(
    result,
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

writeTaskResultEstiObsTrajs <- function(res, opts, info) {
  targetTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  init <- estimateInitialStateAndTime(res, opts, info$task$predictionTime[1], info$task$timeStep)
  outTimes <- seq(
    init$time[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  result <- solveOde(
    u0 = init$initial,
    fun = buildDerivFun(res$hyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$parms)
  resultDenormed <- res$normalization$denormalize(result)
  writeTrajs(
    interpolateTrajs(resultDenormed, targetTimes),
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

estimateInitialStateAndTime <- function(res, opts, startTime, timeStep) {
  name <- getClassAt(opts$initialState, 2)
  if (
    name == "FromTrajs" ||
    (name == "Choose" && isTrajs(res$parms) && startTime == res$parms$time[1])) {
    return(list(time = startTime, initial = getInitialState(res$parms, startTime)))
  } else if (
    name == "FromObs" ||
    (name == "Choose" && (!isTrajs(res$parms) || startTime > res$parms$time[1]))) {
    startObs <- getClosestInTime(res$obsNormed, startTime)
    return(list(time = mean(startObs$time), initial = startObs))
  } else {
    stop("Unkown initOpts ", name)
  }
}

