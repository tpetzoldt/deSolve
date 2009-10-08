## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Car Axis problem (in index 3 formulation)
##        index 3 DAE of dimension 10
##
##     The implementation is as a DLL; code in "andrews.f"
##
## =============================================================================
#-------------------------------------------------------------------------------
#
# before trying this code, the fortran/C programme has to be compiled
# this can be done in R:
# dyn.unload("car.dll");system("R CMD SHLIB car.f")
# dyn.unload("carc.dll");system("R CMD SHLIB carc.c")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#
#-------------------------------------------------------------------------------
require(daeSolve)
dyn.load("carc.dll")       # contains c-code...
#dyn.load("car.dll")       # contains fortra-code...

# initial conditions: state variables
# parameters
pars <- c(eps = 1e-2, M = 10, L = 1, L0 = 0.5,
          r   = 0.1,  w = 10, g = 1)
#
yini <-  with (as.list(pars),
   c(xl=0, yl=L0, xr=L, yr=L0, xla=-L0/L,
     yla=0, xra=-L0/L, yra=0, lam1=0, lam2=0)
              )

# initial conditions: derivates
dyini <- rep(0,10)
# 10 extra output variables (nout)...

FF <- DLLres(res="carresc",time=0.,y=yini,dy=dyini,parms= pars,
             initfunc="carparc", dllname="carc", nout=10)$var
dyini[1:4] <- yini[5:8]
dyini[5:8] <- 2/pars["M"]/(pars["eps"])^2*FF[5:8]

# check consistency of initial condition: delt should be = 0.
DLLres(res="carresc",time=0., y=yini,dy=dyini,parms= pars,
             initfunc="carparc", dllname="carc", nout=10)

# running the model
times <- seq(0,3,by=0.01)
nind  <- c(4,4,2)   # index 1, 2 and 3 variables
# 10 extra output variables...
out   <- mebdfi(y=yini,dy=dyini,times,res="carresc",parms=pars, nind=nind,
                initfunc="carparc", dllname="carc", nout = 10,
                rtol=1e-10, atol=1e-10, outnames=c("f1","f2","f3","f4",
                 "f5","f6","f7","f8","f9","f10"))

par(mar=c(3,3,3,1))
plot(out,type="l",lwd=2,mfrow=c(4,5))

