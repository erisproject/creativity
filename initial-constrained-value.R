R <- matrix(c(-1,-1,10,-10,1,0,-1,0),ncol=2,byrow=T)
r <- matrix(c(-3,10,0.5,0.5))
rint <- r / R[,2]
z <- c(1.2, -0.5)
x <- z[1]
y <- z[2]
plot(x, y, xlim=c(-2,2), ylim=c(-2,4), asp=1)
points(0, 0, pch=4)
P <- list()
M <- list()
for (i in 1:nrow(R)) {
  int <- rint[i,1]
  slope <- -R[i,1]/R[i,2]
  if (slope == Inf || slope == -Inf) abline(v=r[i]/R[i,1])
  else abline(rint[i,1], -R[i,1]/R[i,2])
}

pnext <- z
restrs <- R %*% pnext - r
count <- 0
while (any(restrs > 0) && count < 100) {
  count <- count + 1
  # NB: When doing this for beta values, normalize
  # (although I think that's already being done by the Gibbs sampler!)
  d <- (R %*% pnext - r) / rowSums(R^2)
  dist <- (R %*% pnext - r) / sqrt(rowSums(R^2))
  next_i <- which.max(ifelse(restrs > 0, dist, NA))
  vec <- 1.5 * -d[next_i] * R[next_i,]
  lines(c(pnext[1], pnext[1]+vec[1]), c(pnext[2], pnext[2]+vec[2]), col="red")
  pnext <- pnext + vec
  restrs <- R %*% pnext - r
}
