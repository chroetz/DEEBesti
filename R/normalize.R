unnormalize <- function(trajs, parms) {
  unnormed <- parms$normalization$denormalize(trajs)
  unnormed$time <- unnormed$time / parms$timeScaling
  return(unnormed)
}

normalize <- function(trajs, parms) {
  normed <- parms$normalization$normalize(trajs)
  normed$time <- normed$time * parms$timeScaling
  return(normed)
}
