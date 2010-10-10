
LVmatrix <- function(t, n, parms) {
  with(parms, {
    dn <- r * n + n * (A %*% n)
    return(list(c(dn)))
  })
}
parms <- list(
  r = c(r1 = 0.1, r2 = 0.1, r3 = -0.1, r4 = -0.1),
  A = matrix(c(0.0, 0.0, -0.2, 0.01,      # prey 1
               0.0, 0.0, 0.02, -0.1,      # prey 2
               0.2, 0.02, 0.0, 0.0,       # predator 1; prefers prey 1
               0.01, 0.1, 0.0, 0.0),      # predator 2; prefers prey 2
               nrow = 4, ncol = 4, byrow=TRUE)
)
times <- seq(from = 0, to = 500, by = 0.1)
y     <- c(prey1 = 1, prey2 = 1, pred1 = 2, pred2 = 2)

out <- ode(y, times, LVmatrix, parms)

parms2 <- parms
parms2$r[1] <- 0.2
out2 <- ode(y, times, LVmatrix, parms2)

parms3 <- parms
parms3$r[2] <- 0.2
out3 <- ode(y, times, LVmatrix, parms3)




## Basic line plot
plot(out, type = "l")

## User-specified axis labels
plot(out, type = "l", ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", main = "Time Series")

## add second output 
plot(out, x2=out2, type = "l", ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", main = "Time Series", col = c("red", "blue"))

## add second output
plot(out, x2=list(out2, out3), type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1)
  

 
## ylim is vector
plot(out, x2=list(out2, out3), type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim = c(0,20))

  
## ylim is list
plot(out, x2=list(out2, out3), type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim=list(c(0,20)))
  
## individual ylim
plot(out, x2=list(out2, out3), type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim=list(c(0,15), c(0, 12), c(0,20), c(0, 25)))
  
## incomplete ylim list
plot(out, x2=list(out2, out3), type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim=list(c(0,15), c(0, 12), c(0,20)))
  
  
  
  
