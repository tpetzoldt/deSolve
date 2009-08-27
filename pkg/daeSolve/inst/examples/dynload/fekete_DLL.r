## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Fekete problem (in stabilized index 2 formulation)
##        index 2 DAE of dimension 160
##
##     The implementation is as a DLL; code in "fekete.f"
##
## =============================================================================
# before trying this code, the fortran programme has to be compiled
# this can be done in R:
# dyn.unload("fekete.dll");system("R CMD SHLIB fekete.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
# This does not work WITH the analytic jacobian...
#-------------------------------------------------------------------------------

require(mebdfi)
dyn.load("fekete.dll")

# Initial conditions in a fortran function
Init <- .Fortran("fekinit",N=as.integer(160),
                  T=as.double(0),Y=as.double(rep(0.,160)),
                  dY=as.double(rep(0.,160)))
                  
yini   <- Init$Y
yprime <- Init$dY

# index of the system
nart <- 20
ind  <- c(6*20,2*20,0)
DLLres(res="fekres",time=0.,y=yini,dy=yprime,parms= NULL,
  initfunc=NULL,dllname="fekete")

# run for tend= 1000
times <- c(0,1000)

print("DAE solved with mebdfi - using res, jac, DLL")
print(system.time(
DAE_dll <- mebdfi(y=yini,dy=yprime,times=times,res="fekres", nind=ind,
          dllname="fekete",  initfunc=NULL, parms=NULL,
          hini=1e-8,atol=1e-8,rtol=1e-8,jactype="fullint",maxsteps=100000)
))

# The "correct" solution, a fortran function
sol <- .Fortran("feksoln",N=as.integer(160),
                  T=as.double(1000),Y=as.double(rep(0.,160)))

# the maximal difference
max(abs(sol$Y-DAE_dll[nrow(DAE_dll),-1]))

# run for shorter period, and make plot of first 6 values...
times <- seq(0,20,by=0.1)
print(system.time(
fekete <- mebdfi(y=yini,dy=yprime,times=times,res="fekres", nind=ind,
          dllname="fekete", initfunc=NULL, parms=NULL,
          hini=1e-8,atol=1e-8,rtol=1e-8,maxsteps=100000)
))


dyn.unload("fekete.dll")
plot(fekete,which=1:6,type="l",lwd=2)
