library(deSolve)
library(scatterplot3d)
source("lorenz_R.R")

system("R CMD SHLIB lorenzc.c")
system("R CMD SHLIB lorenzf.f")

dyn.load("lorenzc.dll")
dyn.load("lorenzf.dll")

system.time(
  for(i in 1:10) out   <- ode(state, times, Lorenz, parms = parameters)
)

system.time(
  for(i in 1:10)
    out <- ode(state, times, func = "derivs", parms = parameters,
      dllname = "lorenzc", initfunc = "initmod")
)

system.time(
  for(i in 1:10)
    out <- ode(state, times, func = "derivs", parms = parameters,
      dllname = "lorenzf", initfunc = "initmod")
)

dyn.unload("lorenzc.dll")
dyn.unload("lorenzf.dll")


plot(out)

library(scatterplot3d)
scatterplot3d(out[,-1], type="l")