unnormalize <- function(trajs, parms) {
  unnormedState <- parms$normalization$denormalize(trajs)
  unnormed <- unnormedState$time / parms$timeScaling
  return(unnormed)
}

normalize <- function(trajs, parms) {
  normedState <- parms$normalization$normalize(trajs)
  normed <- normedState$time * parms$timeScaling
  return(normed)
}
