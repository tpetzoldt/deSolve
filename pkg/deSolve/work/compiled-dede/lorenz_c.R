library(deSolve)
library(scatterplot3d)
source("lorenz_R.R")

system("R CMD SHLIB lorenzc.c")
dyn.load("lorenzc.dll")

## non chaotic version for testing
parameters <- c(a = 0.1, b = 0, c =  0)
state <- c(X = 1, Y = 0, Z = 0)

times <- seq(0, 10, 0.1)

out   <- dede(state, 0:1, Lorenz, parms = parameters)

out1 <- ode(state, times=times, func = "derivs", parms = parameters,
  dllname = "lorenzc", initfunc = "initmod", nout = 1, lags = list(mxhist = 1e4))

out2 <- dede(state, times=times, func = "derivs", parms = parameters,
  dllname = "lorenzc", initfunc = "initmod",
  nout = 1, control = list(mxhist = 1e4))
      
      
#out3 <- lsoda(state, times=1:10, func = "derivs", parms = parameters,
#  dllname = "lorenzc", initfunc = "initmod")

dyn.unload("lorenzc.dll")


plot(out2)
