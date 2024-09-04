#options(warn=2) # warnings as error

dbPath <- "~/DeebDbDystsNoisefreeTune_debug"
model <- "Aizawa"
devtools::load_all("~/DEEBesti")

hyperParmsPath <- "~/DEEBextras/hyper/Rnn.json"
hyperParmsList <- ConfigOpts::readOptsBare(hyperParmsPath)
if (!hasValue(hyperParmsList$name)) {
  hyperParmsList$name <- tools::file_path_sans_ext(basename(hyperParmsPath))
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

paths <- DEEBpath::getPaths(dbPath, model)
outPath <- file.path(paths$esti, "debug")

cumMaxErrs <- double(length(hyperParmsList$list))

for (i in seq_along(hyperParmsList$list)) {
  pt <- proc.time()
  applyMethodToModel(
      hyperParms = hyperParmsList$list[[i]],
      observationPath = paths$obs,
      taskPath = paths$task,
      submissionPath = paths$esti,
      verbose = TRUE,
      saveParms = FALSE
  )
  cat(hyperParmsList$name, "took", format((proc.time()-pt)[3]), "s\n")

  truth <- readTrajs(file.path(paths$truth, "task01truth0001.csv"))
  esti <- readTrajs(file.path(outPath, "truth0001obs0001task01esti.csv"))
  taskOpts <- ConfigOpts::readOpts(file.path(paths$task, "task01.json"))
  cumMaxErrs[i] <- DEEBeval:::scoreCumMaxErr(esti, truth, taskOpts$scoreList$list[[1]])
  DEEBplots::plotTimeState(truth[1:20,], esti[1:20,]) |> plot()
}

dirs <- list.dirs(paths$esti)
res <- lapply(dirs, \(d) {
  estiPath <- file.path(d, "truth0001obs0001task01esti.csv")
  if (!file.exists(estiPath)) return(NULL)
  esti <- readTrajs(estiPath)
  taskOpts <- ConfigOpts::readOpts(file.path(paths$task, "task01.json"))
  DEEBeval:::scoreCumMaxErr(esti, truth, taskOpts$scoreList$list[[1]])
})
sort(unlist(res))
idx <- which.min(unlist(res)) + 1
res[[idx]]
dirs[[idx]]
