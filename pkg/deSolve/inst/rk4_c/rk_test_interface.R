library("deSolve")
source("helpfun.R")
source("lvmodel.R")


## compile and load DLL
ext <- if (Sys.info()["sysname"] == "Windows") "dll" else "so"

if (is.loaded("dlotka")) dyn.unload(paste("lotka", ext, sep="."))  # unload DLL if already loaded
system("R CMD SHLIB lotka.c")
lotkadll <- dyn.load(paste("lotka", ext, sep="."))


## lsoda, Solver C, Model R
system.time(
out <- as.data.frame(lsoda(xstart, times, lvmodel, parms,
  hmax=1, atol=1e-6, rtol=1e-6))
)


## rk, Solver R, Model R
system.time(
out <- as.data.frame(rk(xstart, times, lvmodel, parms,  method="ode23",
  hmax=1, atol=1e-6, rtol=1e-6, verbose=F))
)

## rk, Solver C, Model R
system.time(
out.c1 <- as.data.frame(rk_c(xstart, times, lvmodel, parms,  method="ode23",
  hmax=1, atol=1e-6, rtol=1e-6))
)

## Solver C, Model R, fixed step method
system.time(
  out.c <- as.data.frame(rk_c(xstart, times, lvmodel, parms, hmax=0.01, method="rk4"))
)

## Solver C, Model C
system.time(
out.c2 <- as.data.frame(rk_c(xstart, times, func="dlotka", parms, initfunc="ilotka",
  hmax=1, atol=1e-6, rtol=1e-6, nout=2, method="ode23", dllname="lotka"))
)

# two different models
par(mfrow=c(2,1))
matplot(out.c1[,1], out.c[2:4], type="l")

matplot(out.c2[,1], out.c2[2:4], type="l")
