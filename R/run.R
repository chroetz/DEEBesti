applyMethodToModel <- function(
    dbPath,
    statModel,
    example,
    opts,
    hyperParmsList,
    predTimeFromTasks = TRUE,
    obsFileNrs = NULL
) {

  opts <- asOpts(opts, "Estimation")
  hyperParmsList <- asOpts(hyperParmsList, c("HyperParms", "List"))

  basePath <- normalizePath(file.path(dbPath, statModel), mustWork=TRUE)

  if (example) {
    examplePath <- file.path(basePath, "example")
    observationPath <- examplePath
    submissionPath <- examplePath
  } else {
    observationPath <- file.path(basePath, "observation")
    submissionPath <- file.path(basePath, "estimation")
  }

  method <- getClassAt(opts$method, 2)

  outPath <- file.path(submissionPath, method)
  if (!file.exists(outPath)) dir.create(outPath)

  if (predTimeFromTasks) {
    taskFiles <-
      observationPath |>
      dir() |>
      stringr::str_subset("task\\d+.json$")
    # TODO: these lines create too many depenencies
    predictionTime <-
      tibble::tibble(fileNames = taskFiles) |>
      dplyr::mutate(fullPath = file.path(observationPath, .data$fileNames)) |>
      dplyr::mutate(task = lapply(.data$fullPath, ConfigOpts::readOptsBare)) |>
      dplyr::mutate(task = lapply(.data$task, unclass)) |>
      tidyr::unnest_wider(.data$task) |>
      dplyr::pull(.data$predictionTime) |>
      unlist() |>
      range()
    opts <- overwriteOpts(opts, list(outTime = list(range = predictionTime)))
  }

  writeOpts(hyperParmsList, file.path(outPath, "Opts_List_HyperParms")) # TODO
  writeOpts(opts, file.path(outPath, "Opts_Estimation")) # TODO: set default name of out file for write function

  obsFiles <-
    observationPath |>
    dir() |>
    stringr::str_subset("^truth\\d+obs\\d+\\.csv$")

  if (is.null(obsFileNrs)) {
    obsFileSel <- seq_along(obsFiles)
  } else {
    obsFileSel <- intersect(seq_along(obsFiles), obsFileNrs)
  }

  for (obsFile in obsFiles[obsFileSel]) {
    obs <- readTrajs(file.path(observationPath, obsFile))
    res <- estimateWithHyperparameterSelection(
      obs,
      hyperParmsList,
      opts,
      verbose = TRUE)
    outFileBase <- substr(obsFile, 1, nchar(obsFile)-4)
    writeTrajs(res$trajs, file.path(outPath, paste0(outFileBase, "esti.csv")))
    writeOpts(res$hyperParms, file.path(outPath, paste0(outFileBase, "hyperParms")))
  }
}
