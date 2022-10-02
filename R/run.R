applyMethodToModel <- function(
    opts,
    hyperParmsList,
    observationPath = NULL,
    submissionPath = observationPath,
    obsNrFilter = NULL,
    truthNrFilter = NULL
) {

  opts <- asOpts(opts, "Estimation")
  hyperParmsList <- asOpts(hyperParmsList, c("HyperParms", "List"))

  method <- getClassAt(opts$method, 2)

  outPath <- file.path(submissionPath, method)
  if (!file.exists(outPath)) dir.create(outPath)

  predictionTimes <- loadTasksPredictionTime(observationPath)
  totalPredictionTime <- range(unlist(predictionTimes))
  opts <- overwriteOpts(
    opts,
    list(outTime = list(range = totalPredictionTime)))

  writeOpts(hyperParmsList, file.path(outPath, "Opts_List_HyperParms"))
  writeOpts(opts, file.path(outPath, "Opts_Estimation"))

  obsMeta <- getObsMeta(observationPath, obsNrFilter, truthNrFilter)

  for (i in seq_len(nrow(obsMeta))) {
    obs <- readTrajs(obsMeta$path[i])
    res <- estimateWithHyperparameterSelection(
      obs,
      hyperParmsList,
      opts,
      verbose = TRUE)
    baseOfPath <- file.path(outPath, sprintf(
        "truth%04dobs%04d", obsMeta$truthNr[i], obsMeta$obsNr[i]))
    writeOpts(res$hyperParms, paste0(baseOfPath, "hyperParms"))
    for (j in seq_along(predictionTimes)) {
      delta <- getTimeStep(res$trajs$time)
      writeTrajs(
        res$trajs |> dplyr::filter(
          .data$time+delta >= predictionTimes[[j]][1],
          .data$time-delta <= predictionTimes[[j]][2]),
        paste0(baseOfPath, sprintf("task%02desti.csv", j)))
    }
  }
}


getObsMeta <- function(observationPath, obsNrFilter, truthNrFilter) {
  obsFiles <-
    observationPath |>
    dir() |>
    stringr::str_subset("^truth\\d+obs\\d+\\.csv$")
  parts <- stringr::str_match(obsFiles, "^truth(\\d+)obs(\\d+)\\.csv$")
  obsMeta <- tibble::tibble(
    fileName = obsFiles,
    dir = normalizePath(observationPath),
    path = normalizePath(file.path(.data$dir, .data$fileName)),
    obsNr = as.integer(parts[,3]),
    truthNr = as.integer(parts[,2]))
  if (!is.null(obsNrFilter))
    obsMeta <- dplyr::filter(obsMeta, .data$obsNr %in% obsNrFilter)
  if (!is.null(truthNrFilter))
      obsMeta <- dplyr::filter(obsMeta, .data$truthNr %in% truthNrFilter)
  return(obsMeta)
}

loadTasksPredictionTime <- function(observationPath) {
  taskFiles <-
    observationPath |>
    dir() |>
    stringr::str_subset("task\\d+.json$")
  predictionTimes <-
    tibble::tibble(fileNames = taskFiles) |>
    dplyr::mutate(fullPath = file.path(observationPath, .data$fileNames)) |>
    dplyr::mutate(task = lapply(.data$fullPath, ConfigOpts::readOptsBare)) |>
    dplyr::mutate(task = lapply(.data$task, unclass)) |>
    tidyr::unnest_wider(.data$task) |>
    dplyr::pull(.data$predictionTime)
  return(predictionTimes)
}
