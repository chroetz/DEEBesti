#' @export
run <- function(
    dbPath,
    methodOptsDir,
    truthNrFilter = NULL,
    obsNrFilter = NULL,
    modelPattern = NULL,
    methodPattern = NULL,
    copyTruth = FALSE
) {

  if (copyTruth) {
    message("Copying Truth...")
    copyTruth(
      dbPath,
      modelPattern = modelPattern,
      obsNrFilter = obsNrFilter,
      truthNrFilter = truthNrFilter)
  }

  models <- DEEBpath::getModels(dbPath, modelPattern)

  for (model in models) {

    message("MODEL: ", model)
    obsNrs <- DEEBpath::getObservationNrs(dbPath, model, obsNrFilter)

    for (obsNr in obsNrs) {

      message("obsNr: ", obsNr)

      hyperParmsPaths <- DEEBpath::getMethodPaths(methodOptsDir)
      paths <- DEEBpath::getPaths(dbPath, model)

      for (hyperParmsPath in hyperParmsPaths) {
        cat(hyperParmsPath)
        hyperParmsList <- ConfigOpts::readOpts(hyperParmsPath)
        if (nchar(hyperParmsList$name) == 0) {
          x <- basename(hyperParmsPath)
          hyperParmsList$name <- substring(x, 1, x-5)
        }
        pt <- proc.time()
        applyMethodToModel(
          hyperParmsList,
          obsNrFilter = obsNr,
          truthNrFilter = truthNrFilter,
          observationPath = paths$obs,
          submissionPath = paths$esti,
          taskPath = paths$task,
          verbose = FALSE)
        cat(" took ", format((proc.time()-pt)[3]), "s\n", sep="")
      }

    }
  }
}


#' @export
runOne <- function(
    dbPath,
    truthNrFilter = NULL,
    obsNr,
    model,
    method,
    expansionNr = NULL,
    warningAsError = FALSE
) {

  if (warningAsError) options(warn=2)

  hyperParms <- loadHyperParms(
    dbPath,
    method,
    expansionNr)

  pt <- proc.time()
  paths <- DEEBpath::getPaths(dbPath, model)
  applyMethodToModel(
    hyperParms,
    obsNrFilter = obsNr,
    truthNrFilter = truthNrFilter,
    observationPath = paths$obs,
    submissionPath = paths$esti,
    taskPath = paths$task,
    verbose = TRUE)
  cat(" took ", format((proc.time()-pt)[3]), "s\n", sep="")
}


loadHyperParms <- function(
    dbPath,
    methodFile,
    expansionNr = NULL
) {
  cat("loadHyperParms:", methodFile)
  hyperParmsList <- loadAsHyperParmsList(dbPath, methodFile)
  if (is.null(expansionNr)) {
    stopifnot(length(loadHyperParms$list) == 1)
    expansionNr <- 1
  }
  cat(",", expansionNr)
  hyperParms <- hyperParmsList$list[[expansionNr]]
  cat(",", hyperParms$name, "\n")
  return(hyperParms)
}


#' @export
loadAsHyperParmsList <- function(
    dbPath,
    methodFile
) {
  hyperParmsPath <- DEEBpath::getMethodFile(dbPath, methodFile)
  hyperParmsList <- ConfigOpts::readOptsBare(hyperParmsPath)
  if (!hasValue(hyperParmsList$name)) {
    hyperParmsList$name <- basename(methodFile)
  }
  if (ConfigOpts::getClassAt(hyperParmsList, 1) == "List") {
    hyperParmsList <- ConfigOpts::expandList(hyperParmsList)
  } else {
    hyperParmsList <- ConfigOpts::makeOpts(
      c("HyperParms", "List"),
      name = hyperParmsList$name,
      list = list(hyperParmsList))
  }
  for (i in seq_along(hyperParmsList$list)) {
    hyperParmsList$list[[i]]$name <- DEEBpath::nameWithHash(hyperParmsList$name, hyperParmsList$list[[i]])
  }
  return(hyperParmsList)
}


applyMethodToModel <- function(
    hyperParms = NULL,
    observationPath,
    taskPath = observationPath,
    submissionPath = observationPath,
    obsNrFilter = NULL,
    truthNrFilter = NULL,
    verbose = TRUE,
    saveParms = FALSE
) {

  outDir <- file.path(submissionPath, hyperParms$name)
  if (verbose) cat("outDir:", outDir, "\n")
  dir.create(outDir, showWarnings=FALSE, recursive=TRUE)

  writeOpts(hyperParms, dir = outDir, warn = FALSE)

  taskMeta <- DEEBpath::getMetaGeneric(taskPath, tagsFilter = "task")
  meta <- DEEBpath::getMetaGeneric(
    observationPath,
    tagsFilter = c("truth", "obs"),
    nrFilters = list(obsNr = obsNrFilter, truthNr = truthNrFilter))

  for (i in seq_len(nrow(meta))) {
    info <- meta[i,]
    if (verbose) cat(paste0("truth: ", info$truthNr, ", obs: ", info$obsNr, ".\n"))
    estiStartProcTime <- proc.time()
    obs <- readTrajs(info$obsPath)
    pt <- proc.time()
    parms <- getParms(obs, hyperParms)
    normalizedObs <- normalize(obs, parms)
    parmsElapsedTime <- (proc.time() - pt)[3]
    if (saveParms) writeParms(parms, info, outDir)
    taskElapsedTimes <- c()
    for (j in seq_len(nrow(taskMeta))) {
      allInfo <- c(as.list(info), as.list(taskMeta[j,]), list(outDir = outDir))
      pt <- proc.time()
      writeTaskResult(parms, hyperParms, normalizedObs, allInfo)
      taskElapsedTimes[[j]] <- (proc.time() - pt)[3]
    }
    estiElapsedTime <- (proc.time() - estiStartProcTime)[3]
    writeEstiInfo(
      lst(
        estiElapsedTime,
        parmsElapsedTime,
        taskElapsedTimes),
      parms,
      info,
      outDir)
    cleanUpParms(parms)
  }
}


writeEstiInfo <- function(data, parms, info, outDir) {
  estiInfoPath <- file.path(outDir, DEEBpath::estiInfoFile(info))
  info <- as.list(info)
  data <- as.list(data)
  jsonlite::write_json(c(data, info), path=estiInfoPath, pretty=TRUE, digits=I(8), auto_unbox=TRUE)
}


writeParms <- function(parms, info, outDir) {
  if (hasValue(parms$trajs)) {
    smoothPath <- file.path(outDir, DEEBpath::smoothFile(info))
    writeTrajs(unnormalize(parms$trajs, parms), smoothPath)
  }
}


writeTaskResult <- function(parms, hyperParms, normalizedObs, info) {
  info$task <- ConfigOpts::readOptsBare(info$taskPath)
  taskClass <- getClassAt(info$task, 2)
  switch(
    taskClass,
    estiObsTrajs = writeTaskResultEstiObsTrajs(parms, hyperParms, normalizedObs, info),
    newTrajs = writeTaskResultNewTrajs(parms, hyperParms, info),
    velocity = writeTaskResultVelocity(parms, hyperParms, info),
    stop("Unknown task class ", taskClass)
  )
}


writeTaskResultNewTrajs <- function(parms, hyperParms, info) {
  normalizedPredictionTime <- normalizeTime(info$task$predictionTime, parms)
  unnormalizedTargetTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  normalizedTargetTimes <- normalizeTime(unnormalizedTargetTimes, parms)
  init <- makeTrajs(
    time = 0,
    trajId = seq_len(nrow(info$task$initialState)),
    state = info$task$initialState)
  init <- normalize(init, parms)
  normalizedResult <- estimateTrajs(
    init,
    normalizedPredictionTime,
    parms = parms,
    hyperParms)
  if (is.null(normalizedResult)) {
    normalizedResult <- makeTrajs(
      state = matrix(nrow = length(normalizedTargetTimes), ncol = ncol(info$task$initialState)),
      time = normalizedTargetTimes)
  }
  result <- unnormalize(normalizedResult, parms)
  writeTrajs(
    interpolateTrajs(result, unnormalizedTargetTimes),
    file.path(info$outDir, DEEBpath::estiFile(info)))
}


writeTaskResultVelocity <- function(parms, hyperParms, info) {
  gridSides <- lapply(seq_along(info$task$gridSteps), \(i) seq(
    info$task$gridRanges[i,1],
    info$task$gridRanges[i,2],
    info$task$gridSteps[i]
    ))
  grid <- makeDerivTrajs(state = as.matrix(expand.grid(gridSides)))
  gridNormed <- normalize(grid, parms)

  # TODO: move to its own function
  hyperParamsClass <- getClassAt(hyperParms, 2)
  derivs <- switch(
    hyperParamsClass,
    Trajs = {
      derivFun <- buildDerivFun(hyperParms$derivFun)
      t(apply(gridNormed$state, 1, \(s) derivFun(0, s, parms)[[1]]))
    },
    Esn = predictPropagatorDeriv(parms$propagator, hyperParms, gridNormed$state, hyperParms$derivOrder),
    Linear = predictPropagatorDeriv(parms$propagator, hyperParms, gridNormed$state, hyperParms$derivOrder),
    Transformer = predictPropagatorDeriv(parms$propagator, hyperParms, gridNormed$state, hyperParms$derivOrder),
    Direct = predictDirectDeriv(gridNormed$state, parms, hyperParms),
    NeuralOde = predictNeuralOdeDeriv(parms$neuralOde, gridNormed$state),
    stop("Unknown HyperParms subclass")
  )

  resultNormed <- makeDerivTrajs(state = gridNormed$state, deriv = derivs)
  result <- unnormalize(resultNormed, parms) # TODO: this does not work with time normalization
  stopifnot(max(abs(result$state - grid$state)) < sqrt(.Machine$double.eps))
  result$state <- grid$state # more robust against machine imprecision
  writeDerivTrajs(
    result,
    file.path(info$outDir, DEEBpath::estiFile(info)))
}


writeTaskResultEstiObsTrajs <- function(parms, hyperParms, normalizedObs, info) {
  normalizedPredictionTime <- normalizeTime(info$task$predictionTime, parms)
  normalizedTimeStep <- normalizeDuration(info$task$timeStep, parms)
  init <- estimateInitialStateAndTime(
    parms,
    hyperParms,
    normalizedObs,
    normalizedPredictionTime[1],
    normalizedTimeStep)
  unnormalizedTargetTime <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  normalizedTargetTime <- normalizeTime(unnormalizedTargetTime, parms)
  normalizedResult <- estimateTrajs(
    init$initial,
    c(init$time[1], normalizedPredictionTime[2]),
    parms,
    hyperParms)
  if (is.null(normalizedResult)) {
    normalizedResult <- makeTrajs(
      state = matrix(nrow = length(normalizedTargetTime), ncol = ncol(normalizedObs$state)),
      time = normalizedTargetTime)
  }
  unnormalizedResult <- unnormalize(normalizedResult, parms)
  writeTrajs(
    interpolateTrajs(unnormalizedResult, unnormalizedTargetTime),
    file.path(info$outDir, DEEBpath::estiFile(info)))
}


estimateInitialStateAndTime <- function(parms, hyperParms, normalizedObs, startTime, timeStep) {
  name <- getClassAt(hyperParms$initialState, 2)
  if (
    name == "FromTrajs" ||
    (name == "Choose" && "trajs" %in% names(parms) && startTime == parms$trajs$time[1])) {
    return(list(time = startTime, initial = getInitialState(parms$trajs, startTime)))
  } else if (
    name == "FromObs" ||
    (name == "Choose" && (!"trajs" %in% names(parms) || startTime > parms$trajs$time[1]))) {
    startObs <- getClosestInTime(normalizedObs, startTime)
    return(list(time = mean(startObs$time), initial = startObs))
  } else {
    stop("Unkown initOpts ", name)
  }
}

