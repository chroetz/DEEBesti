l2err <- function(trajs, obs) {
  i <- apply(abs(outer(trajs$time, obs$time, `-`)), 2, which.min)
  mean((obs$state - trajs$state[i, ])^2)
}
