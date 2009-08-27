## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Water tube system
##        index 2 DAE of dimension 49
##
##     The implementation is as a DLL; code in "tube.f"
##
## =============================================================================
# before trying this code, the fortran programme has to be compiled
# this can be done in R:
# dyn.unload("tube.dll")
# system("R CMD SHLIB tube.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

require(daeSolve)
dyn.load("tube.dll")

# Initial conditions
yini <- rep(0,49)
yprime <- rep(0,49)
yini[19:36]<- 0.47519404529185289807e-1
yini[37:49] <-109800

pars<- c(nu = 1.31e-6, g = 9.8, rho = 1.0e3, rcrit = 2.3e3,
    length= 1.0e3, k = 2.0e-4, d = 1.0, b = 2.0e2)

# index of the system
ind  <- c(38,11,0)

DLLres(res="tuberes",time=0.,y=yini,dy=yprime,parms= pars,
  initfunc="tubeinit",dllname="tube")

# run for tend= 1000
times <- c(0,17.0*3600)
atol         <- rep(1e-12,49)
atol[37:49]  <- 10000
rtol<-rep(1e-12,49)
print("DAE solved with mebdfi - using res, jac, DLL")
print(system.time(
DAE_dll <- mebdfi(y=yini,dy=yprime,times=times,res="tuberes", nind=ind,
          dllname="tube",initfunc="tubeinit", parms=pars,
          hini=1e-6,atol=atol,rtol=rtol,maxsteps=100000)
))
DAE_dll <- mebdfi(y=yini,dy=yprime,times=times,res="tuberes", nind=ind,
          dllname="daeSolve",initfunc=NULL, parms=NULL,
          hini=1e-2,atol=atol,rtol=rtol,maxsteps=100000)

# The "correct" solution, a fortran function
sol <- .Fortran("tubesoln",N=as.integer(49),
                 Y=as.double(rep(0.,49)))

# the maximal difference
max(abs(sol$Y-DAE_dll[nrow(DAE_dll),-1]))

# small hini terminates R!
times <- seq(0,17.0*3600,by=10)
tuber <- mebdfi(y=yini,dy=yprime,times=times,res="tuberes", nind=ind,
          dllname="tube",initfunc=NULL, parms=NULL,
           atol=atol,rtol=rtol,jactype="fullint",maxsteps=100000)
dyn.unload("tube.dll")
plot(tuber,which=1:6,type="l",lwd=2)
plot(tuber,which=4,type="l",lwd=2,xlim=c(10000,60000),ylim=c(0.000145,0.000185))
