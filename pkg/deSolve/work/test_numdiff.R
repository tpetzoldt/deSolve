library(deSolve)
library(ReacTran)

dx  <- .1   # grid size
v   <- 1     # velocity
x   <- seq(dx/2, 100, by = dx)
N   <- length(x)

#Ntimes <- 201   ## works
Ntimes <- 101   ## works
#Ntimes <- 51    ## works
#Ntimes <- 21    ## works
#Ntimes <- 401   ## FAIL: diffusion!
#Ntimes <- 122   ## FAIL: diffusion!

times <- seq(0, 100, length=Ntimes)

## hard test with irregular time steps, FAIL
#times <- sort(c(0, 100, runif(Ntimes-2, min=0, max=100))) 

Y0 <- c(2, 5, 10, 5, 2, rep(0, N-5))

## conventional
transport <- function(t, Y, parms) {
  dY  <- - as.vector(v * diff(c(1e-6, Y)) / dx)
  list(c(dY))
}

## ReacTran
#xgrid <- setup.grid.1D(x.up=0, L=dx*N, N=N)

transport <- function(t, Y, parms) {
  cat(t, ": ", timestep(TRUE), " -- ", timestep(FALSE), "\n")
  #dY <- tran.1D(C=Y, C.up=1e-6, v=v,dx=dx)
  dY <- advection.1D(C=Y, C.up=0, v=v,dx=dx, adv.method="muscl")
  list(c(dY$dC))
}


## hini: NULL, 1, 0.5 works; 0.25: diffusion; other: diffusion or unstable
system.time(
out   <- ode.1D(y = Y0, times, transport, method="euler", hini=0.5,
  parms = NULL, nspec = 1, ynames=FALSE, rtol=1e-4, atol=1e-4)
)

## euler.1D uses "special euler" code
system.time(
  out2   <- euler.1D(y = Y0, times, transport,  parms = NULL, nspec = 1)
)


image(out)
plot.1D(out, type="l", ylim=c(-10, 12), delay = 20)

image(out2)
plot.1D(out2, type="l", ylim=c(-10, 12), delay = 20)

  
  