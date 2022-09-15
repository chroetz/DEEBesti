plot_traj <- function(traj, truth, obs, obs_vali) {

  u_all <- rbind(traj$state, truth$state, obs$state, obs_vali$state)
  par(mar = c(0.1,0.1,0.1,0.1))
  plot(
    NA,
    xlim=range(u_all[,1]),
    ylim=range(u_all[,2]),
    xlab = NA, ylab = NA, axes = F, asp = 1)
  box()
  grid()
  lines(truth$state, col=1, lwd=2)

  lines(z_esti$state, col=5, lwd=2)
  lines(traj$state[traj$time<=max(obs$time),], col=2, lwd=2)

  points(obs$state,col=2)
  iobstru <- apply(abs(outer(truth$time, obs$time, `-`)), 2, which.min)
  segments(obs$state[,1], obs$state[,2], truth$state[iobstru,1], truth$state[iobstru,2], col=2)

  points(obs_vali$state,col=3)
  iobsvtru <- apply(abs(outer(truth$time, obs_vali$time, `-`)), 2, which.min)
  segments(obs_vali$state[,1], obs_vali$state[,2], truth$state[iobsvtru,1], truth$state[iobsvtru,2], col=3)

  tmaxobs <- max(c(obs$time, obs_vali$time))
  lines(traj$state[traj$time<=tmaxobs,], col=2, lwd=2)
  lines(traj$state[traj$time>tmaxobs,], col=4, lwd=2)

}
