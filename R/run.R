applyMethodToModel <- function(
    opts,
    hyperParmsList = NULL,
    observationPath,
    taskPath = observationPath,
    submissionPath = observationPath,
    obsNrFilter = NULL,
    truthNrFilter = NULL,
    verbose = TRUE
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
      verbose = verbose)
    hpPath <- file.path(outDir, DEEBpath::hyperParmsFile(info))
    writeOpts(res$hyperParms, hpPath)

    for (j in seq_len(nrow(taskMeta))) {
      allInfo <- c(as.list(info), as.list(taskMeta[j,]), list(outDir = outDir))
      writeTaskResult(res$parms, res$hyperParms, obs, opts, allInfo)
    }

    cleanUpParms(res$parms)
  }
}

writeTaskResult <- function(parms, hyperParms, obs, opts, info) {
  info$task <- ConfigOpts::readOptsBare(info$taskPath)
  taskClass <- getClassAt(info$task, 2)
  switch(
    taskClass,
    estiObsTrajs = writeTaskResultEstiObsTrajs(parms, hyperParms, obs, opts, info),
    newTrajs = writeTaskResultNewTrajs(parms, hyperParms, opts, info),
    velocity = writeTaskResultVelocity(parms, hyperParms, opts, info),
    stop("Unknown task class ", taskClass)
  )
}

writeTaskResultNewTrajs <- function(parms, hyperParms, opts, info) {
  odeTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    length.out = opts$odeSteps)
  targetTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  init <- makeTrajs(
    time = 0,
    trajId = seq_len(nrow(info$task$initialState)),
    state = info$task$initialState)
  init <- parms$normalization$normalize(init)
  result <- solveOde(
    u0 = init,
    fun = buildDerivFun(hyperParms$derivFun),
    times = odeTimes,
    opts = opts$odeSolver,
    parms = parms)
  result <- parms$normalization$denormalize(result)
  writeTrajs(
    interpolateTrajs(result, targetTimes),
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

writeTaskResultVelocity <- function(parms, hyperParms, opts, info) {
  gridSides <- lapply(seq_along(info$task$gridSteps), \(i) seq(
    info$task$gridRanges[i,1],
    info$task$gridRanges[i,2],
    info$task$gridSteps[i]
    ))
  grid <- makeDerivTrajs(state = as.matrix(expand.grid(gridSides)))
  gridNormed <- parms$normalization$normalize(grid)
  derivFun <- buildDerivFun(hyperParms$derivFun)
  derivs <- t(apply(gridNormed$state, 1, \(s) derivFun(0, s, parms)[[1]]))
  resultNormed <- makeDerivTrajs(state = gridNormed$state, deriv = derivs)
  result <- parms$normalization$denormalize(resultNormed)
  writeDerivTrajs(
    result,
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

writeTaskResultEstiObsTrajs <- function(parms, hyperParms, obs, opts, info) {
  init <- estimateInitialStateAndTime(parms, hyperParms, obs, info$task$predictionTime[1], info$task$timeStep)
  odeTimes <- seq(
    init$time[1],
    info$task$predictionTime[2],
    length.out = opts$odeSteps)
  targetTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  result <- solveOde(
    u0 = init$initial,
    fun = buildDerivFun(hyperParms$derivFun),
    times = odeTimes,
    opts = opts$odeSolver,
    parms = parms)
  resultDenormed <- parms$normalization$denormalize(result)
  writeTrajs(
    interpolateTrajs(resultDenormed, targetTimes),
    file.path(info$outDir, DEEBpath::estiFile(info)))
}

estimateInitialStateAndTime <- function(parms, hyperParms, obs, startTime, timeStep) {
  name <- getClassAt(hyperParms$initialState, 2)
  if (
    name == "FromTrajs" ||
    (name == "Choose" && "trajs" %in% names(parms) && startTime == parms$trajs$time[1])) {
    return(list(time = startTime, initial = getInitialState(parms$trajs, startTime)))
  } else if (
    name == "FromObs" ||
    (name == "Choose" && (!"trajs" %in% names(parms) || startTime > parms$trajs$time[1]))) {
    startObs <- getClosestInTime(obs, startTime)
    startObs <- parms$normalization$normalize(startObs)
    return(list(time = mean(startObs$time), initial = startObs))
  } else {
    stop("Unkown initOpts ", name)
  }
}

