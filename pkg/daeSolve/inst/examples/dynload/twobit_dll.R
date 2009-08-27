## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Two bit adding unit
##        index 1 DAE of dimension 350
##
##     The implementation is as a DLL; code in "tube.f"
##
## =============================================================================
# before trying this code, the fortran programme has to be compiled
# this can be done in R:
# dyn.unload("twobit.dll");system("R CMD SHLIB twobit.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

require(daeSolve)
dyn.load("twobit.dll")

N <- 350

# Initial conditions in a fortran function
Init <- .Fortran("twobinit",N=as.integer(N),
                  T=as.double(0),Y=as.double(rep(0.,N)),
                  dY=as.double(rep(0.,N)))
yini   <- Init$Y
yprime <- Init$dY

# index of the system
ind  <- c(N,0,0)

DLLres(res="twobres",time=0.,y=yini,dy=yprime,parms= NULL,
  initfunc=NULL,dllname="twobit")

# run for tend= 320
times <- c(0,320)

print("DAE solved with mebdfi - using res, jac, DLL")
print(system.time(
DAE_dll <- mebdfi(y=yini,dy=yprime,times=times,res="twobres", nind=ind,
          dllname="twobit",  initfunc=NULL, parms=NULL,
          hini=1e-4,atol=1e-5,rtol=1e-5,maxsteps=100000)
))

# The "correct" solution, a fortran function
sol <- .Fortran("twobsoln",N=as.integer(N),
                  Y=as.double(rep(0.,N)))

# the maximal difference
max(abs(sol$Y-DAE_dll[nrow(DAE_dll),-1]))

# run dynamically
times <- seq(0,320,by=0.1)
print(system.time(
twobit <- mebdfi(y=yini,dy=yprime,times=times,res="twobres", nind=ind,
          dllname="twobit", initfunc=NULL, parms=NULL,
          hini=1e-4,atol=1e-4,rtol=1e-4,maxsteps=100000)
))

dyn.unload("twobit.dll")
#plot(twobit,type="l",lwd=2)
plot(twobit,which=c(49,130,148),type="l",lwd=2,mfrow=c(3,1))

