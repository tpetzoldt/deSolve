## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Wheelset in index-2 formulation
##        index 2 IDE of dimension 17
##
##     The implementation is as a DLL; code in "wheelset.f"
##
## =============================================================================
# before trying this code, the fortran programme has to be compiled
# this can be done in R:
# dyn.unload("wheelset.dll")
# system("R CMD SHLIB wheelset.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

require(daeSolve)
dyn.load("wheelset.dll")


# index of the system
ind  <- c(15,2,0)

# initial conditions
yini <- c( 0.14941e-02,0.40089e-06,0.11241e-05,-.28573e-03,
           0.26459e-03,0,0,0,0,0,0, -7.4122380357667139e-06,
           -0.1521364296121248,7.5634406395172940e-06,0.1490635714733819,
           -8.3593e-3,-7.4144e-3)
yprime <- c(0,0,0,0,0,-1.975258894011285,-1.0898297102811276e-03,
           7.8855083626142589e-02,-5.533362821731549,-0.3487021489546511,
           -2.132968724380927,0,0,0,0,0,0)

DLLres(res="wheelres",time=0.,y=yini,dy=yprime,parms= NULL,
  initfunc=NULL,dllname="wheelset")


# run for tend= 1000
times <- c(0,10)
atol         <- rep(1e-6,17)
atol[16:17]  <- 1e10
rtol<-rep(1e-6,17)
rtol[16:17]  <- 1e10

print("DAE solved with mebdfi - using res, jac, DLL")
print(system.time(
DAE_dll <- mebdfi(y=yini,dy=yprime,times=times,res="wheelres", nind=ind,
          dllname="wheelset",initfunc=NULL, parms=NULL,
          hini=1e-8,atol=atol,rtol=rtol,maxsteps=100000)
))

# The "correct" solution, a fortran function
sol <- .Fortran("wheelsoln",N=as.integer(17),t=as.double(0.),
                 Y=as.double(rep(0.,17)))

# the maximal difference
max(abs(sol$Y-DAE_dll[nrow(DAE_dll),-1]))

# small hini terminates R!
times <- seq(0,4,by=0.001)
wheel <- mebdfi(y=yini,dy=yprime,times=times,res="wheelres", nind=ind,
          dllname="wheelset",initfunc=NULL, parms=NULL,
          hini=1e-15,atol=atol,rtol=rtol,maxsteps=100000)

plot(wheel,which=1:6,type="l",lwd=2)
dyn.unload("wheelset.dll")
