library("deSolve")
source("helpfun.R")
source("lvmodel.R")


longtime  <- seq(0, 10000, length=100001)

## compile and load DLL
ext <- if (Sys.info()["sysname"] == "Windows") "dll" else "so"

if (is.loaded("dlotka")) dyn.unload(paste("lotka", ext, sep="."))  # unload DLL if already loaded
system("R CMD SHLIB lotka.c")
lotkadll <- dyn.load(paste("lotka", ext, sep="."))


## lsoda, Solver C, Model R
system.time(
out <- lsoda(xstart, times, lvmodel, parms,
  hmax=1, atol=1e-6, rtol=1e-6)
)

## lsoda, Solver C, Model C, long run
system.time(
out.l1 <- lsoda(xstart, longtime, func="dlotka", parms, initfunc="ilotka",
  hmax=1, atol=1e-12, rtol=1e-8, nout=2,dllname="lotka")
)


## rk, Solver R, Model R
system.time(
out <- rk_r(xstart, times, lvmodel, parms,  method="ode23",
  hmax=1, atol=1e-6, rtol=1e-6, verbose=FALSE)
)

## rk, Solver C, Model R
system.time(
out.c1 <- rk(xstart, times, lvmodel, parms,  method="ode23",
  hmax=1, atol=1e-6, rtol=1e-6)
)

## Solver C, Model R, fixed step method
system.time(
  out.c <- rk(xstart, times, lvmodel, parms, hmax=0.01, method="rk4")
)

## Solver C, Model C
system.time(
out.c2 <- rk(xstart, times, func="dlotka", parms, initfunc="ilotka",
  hmax=1, atol=1e-6, rtol=1e-6, nout=2, method="ode23", dllname="lotka")
)

## Solver C, Model C, long run
system.time(
out.l2 <- rk(xstart, longtime, func="dlotka", parms, initfunc="ilotka",
  hmax=1, atol=1e-12, rtol=1e-8, nout=2, method="ode45", dllname="lotka", maxsteps=5e5)
)

# two different models
out.c1 <- as.data.frame(out.c1)
out.c2 <- as.data.frame(out.c2)

par(mfrow=c(2,1))
matplot(out.c1[,1], out.c1[2:4], type="l", lwd=3)

matplot(out.c2[,1], out.c2[2:4], type="l", lwd=3)

diff <- out.l1 - out.l2
summary(diff)
