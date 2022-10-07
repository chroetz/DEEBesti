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
      writeTaskResult(res, opts, info)
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
  iniState <- info$task$initialState
  rownames(iniState) <- seq_len(nrow(iniState))
  result <- solveOde(
    u0 = iniState,
    fun = buildDerivFun(res$hyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$trajs)
  writeTrajs(
    result,
    file.path(info$outDir, sprintf(
      "truth%04dobs%04dtask%02desti.csv",
      info$truthNr, info$obsNr, info$taskNr)))
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
    file.path(info$outDir, sprintf(
      "truth%04dobs%04dtask%02desti.csv",
      info$truthNr, info$obsNr, info$taskNr)))
}

writeTaskResultEstiObsTrajs <- function(res, opts, info) {
  outTimes <- seq(
    info$task$predictionTime[1],
    info$task$predictionTime[2],
    by = info$task$timeStep)
  result <- solveOde(
    u0 = getInitialState(res$trajs, info$task$predictionTime[1]),
    fun = buildDerivFun(res$hyperParms$derivFun),
    times = outTimes,
    opts = opts$odeSolver,
    parms = res$trajs)
  writeTrajs(
    result,
    file.path(info$outDir, sprintf(
      "truth%04dobs%04dtask%02desti.csv",
      info$truthNr, info$obsNr, info$taskNr)))
}
