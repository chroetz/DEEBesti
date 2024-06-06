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
    model,
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
    model,
    method,
    expansionNr
) {
  hyperParmsPath <- DEEBpath::getMethodFile(dbPath, method)
  paths <- DEEBpath::getPaths(dbPath, model)
  cat(hyperParmsPath)
  hyperParmsObject <- ConfigOpts::readOptsBare(hyperParmsPath)
  if (nchar(hyperParmsObject$name) == 0) hyperParmsObject$name <- basename(method)
  if (ConfigOpts::getClassAt(hyperParmsObject, 1) == "List") {
    if (is.null(expansionNr)) {
      stopifnot(length(hyperParmsObject$list) == 1)
      expansionNr <- 1
    }
    hyperParmsList <- ConfigOpts::expandList(hyperParmsObject)
    cat(",", expansionNr)
    hyperParms <- hyperParmsList$list[[expansionNr]]
    hyperParms$name <- DEEBpath::nameWithHash(hyperParmsList$name, hyperParms)
  } else {
    hyperParms <- hyperParmsObject
  }
  cat(",", hyperParms$name)
  return(hyperParms)
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
    obs <- readTrajs(info$obsPath)
    parms <- getParms(obs, hyperParms)
    if (saveParms) writeParms(parms, obs, info, outDir)
    for (j in seq_len(nrow(taskMeta))) {
      allInfo <- c(as.list(info), as.list(taskMeta[j,]), list(outDir = outDir))
      writeTaskResult(parms, hyperParms, obs, allInfo)
    }
    cleanUpParms(parms)
  }
}


writeParms <- function(parms, obs, info, outDir) {
  if (!is.null(parms$trajs)) {
    smoothPath <- file.path(outDir, DEEBpath::smoothFile(info))
    writeTrajs(unnormalize(parms$trajs, parms), smoothPath)
  }
}


writeTaskResult <- function(parms, hyperParms, obs, info) {
  info$task <- ConfigOpts::readOptsBare(info$taskPath)
  taskClass <- getClassAt(info$task, 2)
  switch(
    taskClass,
    estiObsTrajs = writeTaskResultEstiObsTrajs(parms, hyperParms, obs, info),
    newTrajs = writeTaskResultNewTrajs(parms, hyperParms, info),
    velocity = writeTaskResultVelocity(parms, hyperParms, info),
    stop("Unknown task class ", taskClass)
  )
}


writeTaskResultNewTrajs <- function(parms, hyperParms, info) {
  targetTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  init <- makeTrajs(
    time = 0,
    trajId = seq_len(nrow(info$task$initialState)),
    state = info$task$initialState)
  init <- normalize(init, parms)
  result <- estimateTrajs(
    init,
    info$task$predictionTime,
    parms = parms,
    hyperParms)
  if (is.null(result)) {
    writeTrajs(
      makeTrajs(
        state = matrix(nrow = length(targetTimes), ncol = ncol(info$task$initialState)),
        time = targetTimes),
      file.path(info$outDir, DEEBpath::estiFile(info)))
    return(invisible(NULL))
  }
  result <- unnormalize(result, parms)
  writeTrajs(
    interpolateTrajs(result, targetTimes),
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
    Esn = predictEsnDeriv(parms$esn, gridNormed$state, hyperParms$derivOrder),
    Linear = predictLinearDeriv(parms$linear, gridNormed$state, hyperParms$derivOrder),
    Direct = predictDirectDeriv(gridNormed$state, parms, hyperParms),
    Transformer = stop("Deriv not implemented for Transformer"),
    stop("Unknown HyperParms subclass")
  )

  resultNormed <- makeDerivTrajs(state = gridNormed$state, deriv = derivs)
  result <- unnormalize(resultNormed, parms)
  stopifnot(max(abs(result$state - grid$state)) < sqrt(.Machine$double.eps))
  result$state <- grid$state # more robust against machine imprecision
  writeDerivTrajs(
    result,
    file.path(info$outDir, DEEBpath::estiFile(info)))
}


writeTaskResultEstiObsTrajs <- function(parms, hyperParms, obs, info) {
  init <- estimateInitialStateAndTime(
    parms,
    hyperParms,
    obs,
    info$task$predictionTime[1],
    info$task$timeStep)
  targetTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  result <- estimateTrajs(
    init$initial,
    c(init$time[1], info$task$predictionTime[2]),
    parms,
    hyperParms)
  if (is.null(result)) {
    writeTrajs(
      makeTrajs(
        state = matrix(nrow = length(targetTimes), ncol = ncol(obs$state)),
        time = targetTimes),
      file.path(info$outDir, DEEBpath::estiFile(info)))
    return(invisible(NULL))
  }
  resultDenormed <- unnormalize(result, parms)
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
    startObs <- normalize(startObs, parms)
    return(list(time = mean(startObs$time), initial = startObs))
  } else {
    stop("Unkown initOpts ", name)
  }
}

