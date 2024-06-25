unnormalize <- function(trajs, parms) {
  unnormed <- parms$normalization$denormalize(trajs)
  unnormed$time <- unnormalizeTime(unnormed$time, parms)
  return(unnormed)
}

normalize <- function(trajs, parms) {
  normed <- parms$normalization$normalize(trajs)
  normed$time <- normalizeTime(normed$time, parms)
  return(normed)
}

normalizeTime <- function(time, parms) {
  time * parms$timeScaling
}

unnormalizeTime <- function(time, parms) {
  time / parms$timeScaling
}

normalizeDuration <- function(time, parms) {
  time * parms$timeScaling
}

unnormalizeDuration <- function(time, parms) {
  time / parms$timeScaling
}
