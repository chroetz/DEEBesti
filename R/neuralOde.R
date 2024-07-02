createAndTrainNeuralOde <- function(opts, obs) {

  workingDirRoot <- normalizePath(".DEEB_tmp", mustWork=FALSE)
  workingName <- DEEButil::getUniqueFileName(workingDirRoot, identifyingObject=list(opts, obs, Sys.time()))
  workingDir <- file.path(workingDirRoot, workingName)
  dir.create(workingDir, showWarnings=FALSE, recursive=TRUE)

  wd <- getwd()
  on.exit(setwd(wd))

  cat("cd to", workingDir, "\n")
  setwd(workingDir)

  DEEBtrajs::writeTrajs(obs, "obs.csv")

  projectPath <- getDeebJlPath(opts)
  scriptPath <- file.path(projectPath, "train.jl")
  schedulePath <- file.path(projectPath, paste0("schedule_", opts$learningRateSchedule, ".toml"))

  cmd <- paste0(
    'julia',
    ' --project="', projectPath, '"',
    ' "', scriptPath, '"',
    ' --rng-seed ', opts$seed,
    ' --dim ', ncol(obs$state),
    ' --data-path "obs.csv"',
    ' --hidden-layers ', opts$hiddenLayers,
    ' --hidden-width ', opts$hiddenWidth,
    ' --activation ', opts$activation,
    ' --steps ', opts$steps,
    ' --train-frac ', opts$trainFrac,
    ' --optimiser-rule AdamW',
    ' --optimiser-hyperparams "lambda=', opts$weightDecay,'"',
    ' --epochs ', opts$epochs,
    ' --schedule-file "', schedulePath, '"',
    ' --sensealg BacksolveAdjoint',
    ' --vjp ZygoteVJP')

  cat("Run command:\n", cmd, "\n")

  errorCode <- system(cmd)
  stopifnot(errorCode == 0)

  return(lst(workingDir, opts))
}


getDeebJlPath <- function(opts) {
  normalizePath(opts$deebNeuralOdeProjectPath, mustWork=FALSE)
}



predictNeuralOde <- function(neuralOde, startState, timeRange, timeStep) {

  stopifnot(length(timeStep) == 1)
  stopifnot(is.numeric(timeStep))
  stopifnot(is.finite(timeStep))
  stopifnot(timeStep > 0)

  opts <- neuralOde$opts
  stateDim <- ncol(startState)

  wd <- getwd()
  on.exit(setwd(wd))

  cat("cd to", neuralOde$workingDir, "\n")
  setwd(neuralOde$workingDir)

  startTrajs <- DEEBtrajs::makeTrajs(timeRange[1], matrix(startState, nrow=1))
  DEEBtrajs::writeTrajs(startTrajs, "start.csv")

  projectPath <- getDeebJlPath(opts)
  scriptPath <- file.path(projectPath, "infer.jl")

  cmd <- paste0(
    'julia',
    ' --project="', projectPath, '"',
    ' "', scriptPath, '"',
    ' --dim ', stateDim,
    ' --data-path "start.csv"',
    ' --hidden-layers ', opts$hiddenLayers,
    ' --hidden-width ', opts$hiddenWidth,
    ' --activation ', opts$activation,
    ' --epochs ', 0,
    ' --schedule-file "-"',
    ' --pred-traj',
    ' --pred-path "esti.csv"',
    ' --t0 ', timeRange[1],
    ' --t1 ', timeRange[2],
    ' --dt ', timeStep)

  cat("Run command:\n", cmd, "\n")

  errorCode <- system(cmd)
  stopifnot(errorCode == 0)

  esti <- DEEBtrajs::readTrajs("esti.csv")

  return(esti)
}






predictNeuralOdeDeriv <- function(neuralOde, grid) {

  opts <- neuralOde$opts

  wd <- getwd()
  on.exit(setwd(wd))

  cat("cd to", neuralOde$workingDir, "\n")
  setwd(neuralOde$workingDir)

  startTrajs <- DEEBtrajs::makeTrajs(rep(0, nrow(grid)), grid, seq_len(nrow(grid)))
  stateDim <- DEEBtrajs::getDim(startTrajs)
  DEEBtrajs::writeTrajs(startTrajs, "grid.csv")

  projectPath <- getDeebJlPath(opts)
  scriptPath <- file.path(projectPath, "infer.jl")

  cmd <- paste0(
    'julia',
    ' --project="', projectPath, '"',
    ' "', scriptPath, '"',
    ' --dim ', stateDim,
    ' --data-path "grid.csv"',
    ' --hidden-layers ', opts$hiddenLayers,
    ' --hidden-width ', opts$hiddenWidth,
    ' --activation ', opts$activation,
    ' --epochs ', 0,
    ' --schedule-file "-"',
    ' --pred-grid',
    ' --pred-path "estigrid.csv"')

  cat("Run command:\n", cmd, "\n")

  errorCode <- system(cmd)
  stopifnot(errorCode == 0)

  esti <- DEEBtrajs::readDerivTrajs("estigrid.csv")

  return(esti)
}
