library(deSolve)
library(ReacTran)

## conventional
transport <- function(t, Y, parms) {
  dY  <- - as.vector(v * diff(c(1e-6, Y)) / dx)
  list(c(dY))
}

## ReacTran

transport2 <- function(t, Y, parms) {
#  cat(t, ": ", timestep(TRUE), " -- ", timestep(FALSE), "\n")
  dY <- advection.1D(C=Y, C.up=0, v=v,dx=dx, adv.method="muscl")
  list(c(dY$dC))
}

dx  <- .1    # grid size
v   <- 1     # velocity
x   <- seq(dx/2, 100, by = dx)
N   <- length(x)
Y0  <- c(2, 5, 10, 5, 2, rep(0, N-5))


testfunc <- function (dt, func = transport2, dt2 = dt) {
  times <- seq(0, 50, by = dt2)
  out1 <- ode.1D(y = Y0, times, func, method="euler", hini=dt,
                parms = NULL, nspec = 1, ynames=FALSE, rtol=1e-5, atol=1e-5)
  out2 <- euler.1D(y = Y0, times, func,  parms = NULL, nspec = 1)
  plot(out1[nrow(out1),-1])
  lines(out2[nrow(out2),-1])
}
testfunc(0.5, transport2, 1.25)



times <- seq(0, 50, by = 1)
out1 <- ode.1D(y = Y0, times, transport2, method="euler", hini=2,
               parms = NULL, nspec = 1)
out2 <- euler.1D(y = Y0, times, transport2,  parms = NULL, nspec = 1)
plot(out1[nrow(out1),-1], type = "l", col = "red")
lines(out2[nrow(out2),-1], col="grey")

