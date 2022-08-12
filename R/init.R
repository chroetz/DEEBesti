#' @export
init <- function(obs, inter_steps) {
  n <- nrow(obs$u)
  inter_steps <- inter_steps
  steps <- (n-1)*inter_steps + 1
  d <- ncol(obs$u)
  t <- seq(min(obs$t), max(obs$t), length.out = steps)
  step_size <- t[2] - t[1]
  u <- cbind(
    stats::approx(obs$t, obs$u[,1], t)$y,
    stats::approx(obs$t, obs$u[,2], t)$y)
  a <- cbind(
    (u[2:steps,1]-u[1:(steps-1),1]) / step_size,
    (u[2:steps,2]-u[1:(steps-1),2]) / step_size)
  has_obs <- 1:steps %% inter_steps == 1
  lst <- as.list(environment())
  return(lst[setdiff(names(lst), "obs")])
}
