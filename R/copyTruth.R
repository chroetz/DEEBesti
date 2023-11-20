#' @export
copyTruth <- function(
    dbPath,
    modelPattern = NULL,
    obsNrFilter = NULL,
    truthNrFilter = NULL
) {
  cat("Start Copy Truth.\n")
  models <- DEEBpath::getModels(dbPath, modelPattern)
  for (model in models) {
    paths <- DEEBpath::getPaths(dbPath, model)
    cat("Model:", model, "\n")
    meta <- DEEBpath::getMetaGeneric(
      c(paths$obs, paths$truth),
      tagFileFilter = list(c("truth", "obs"), c("task", "truth")),
      nrFilters = list(obsNr = obsNrFilter, truthNr = truthNrFilter))
    outPath <- file.path(paths$esti, "Truth")
    dir.create(outPath, showWarnings=FALSE, recursive=TRUE)
    for (i in seq_len(nrow(meta))) {
      info <- as.list(meta[i,])
      cat("\ttruth:", info$truthNr, ", obs:", info$obsNr, ", task:", info$taskNr, "\n")
      file.copy(info$truthPath, file.path(outPath, DEEBpath::estiFile(info)))
    }
  }
  cat("End Copy Truth.\n")
}

