R <- matrix(c(-1,-1,1,-1,1,0,-1,0,0,1),ncol=2,byrow=T)
r <- matrix(c(-3,-2.5,0.5,0.75,3.3))
# Intersections:
poly_x <- c(-0.3, 0.5, 0.5, 0.25)
poly_y <- c( 3.3, 3.3, 3.0, 2.75)
rint <- r / R[,2]
z <- c(1, 1)
x <- z[1]
y <- z[2]
plot(x, y, xlim=c(-2,2), ylim=c(0,4), asp=1)
polygon(poly_x, poly_y, lwd=2)
P <- list()
M <- list()
for (i in 1:nrow(R)) {
  int <- rint[i,1]
  slope <- -R[i,1]/R[i,2]
  if (slope == Inf || slope == -Inf) abline(v=r[i]/R[i,1], col="gray")
  else abline(rint[i,1], -R[i,1]/R[i,2], col="gray")
}

pnext <- z
restrs <- R %*% pnext - r
count <- 0

fixsample <- function(x, ...) x[sample.int(length(x), ...)]
while (any(restrs > 0)) {
  count <<- count + 1
  d <- (R %*% pnext - r) / rowSums(R^2)
  next_i <- fixsample(which(restrs > 0), 1)
  vec <- 1.5 * -d[next_i] * R[next_i,]
  arrows(pnext[1], pnext[2], pnext[1]+vec[1], pnext[2]+vec[2], col="red")
  pnext <<- pnext + vec
  restrs <<- R %*% pnext - r
}
