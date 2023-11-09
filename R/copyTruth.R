#' @export
copyTruth <- function(
    dbPath,
    modelPattern = NULL,
    obsNrFilter = NULL,
    truthNrFilter = NULL,
    example = FALSE
) {
  models <- DEEBpath::getModels(dbPath, modelPattern)
  for (model in models) {
    paths <- DEEBpath::getPaths(dbPath, model, example=example)
    cat(model, "\n")
    truthMeta <- DEEBpath::getMetaGeneric(
      paths$truth,
      c("task", "truth"),
      nrFilters = list(obsNr = obsNrFilter, truthNr = truthNrFilter))
    obsMeta <- DEEBpath::getMetaGeneric(
      paths$obs,
      c("truth", "obs"),
      nrFilters = list(obsNr = obsNrFilter, truthNr = truthNrFilter))
    meta <- dplyr::left_join(obsMeta, truthMeta, by = "truthNr", multiple = "all")
    outPath <- file.path(paths$esti, "Truth")
    dir.create(outPath, showWarnings=FALSE, recursive=TRUE)
    for (i in seq_len(nrow(meta))) {
      info <- as.list(meta[i,])
      file.copy(info$truthPath, file.path(outPath, DEEBpath::estiFile(info)))
    }
  }
}

