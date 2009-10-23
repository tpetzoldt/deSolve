## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Slider crank, index 2 IDE, 24 eqns
##
##     The implementation is as a DLL; code in "crank.f"
##
## =============================================================================

#---------------------------------------------------------------------------
# before trying this code, the fortran programme has to be compiled
# this can be done in R:
# dyn.unload("crank.dll")
# system("R CMD SHLIB crank.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

require(daeSolve)
dyn.load("crank.dll")

# 2 possible initial conditions...

yini2 <- rep(0,24)
yini2[c(3,8,9,17,20,21,23)]<-c(0.45,150,-75,-3.789473684210526e3,
     1.924342105263158e2,1.273026315789474e3,2.863023157894737e2)

## First initial conditions
y      <- yini2
names(y) <- c("phi1","phi2","x3","q1","q2","q3","q4",
              "vphi1","vphi2","vx3","vq1","vq2","vq3","vq4",
              "aphi1","aphi2","ax3","aq1","aq2","aq3","aq4",
              "la1","la2","la3")
yprime <- c(y[8:21],rep(0,10))
ipar   <- c(0,0)
nind   <- c(14,10,0)

DLLres(res="crankres",time=0.,y=y,dy=yprime,parms= NULL,
  initfunc=NULL,dllname="crank")
Tol <- 1e-5
times <- seq(0,0.1,by=0.001)
crank1 <- mebdfi(y=y,dy=yprime,times=times,res="crankres", nind=nind,
          dllname="crank",initfunc=NULL, parms=NULL, ipar=ipar,
          hini=Tol*1e-2,atol=Tol,rtol=Tol,maxsteps=100000)
plot(crank1, type="l", lwd=2,
  which=c("phi2","x3","q1","q2","q3","q4","la1","la2","la3"))

## Second set of initial conditions
y <- c(0,0,0.450016933,0,0,0.103339863e-4,0.169327969e-4,
       0.150000000e3,-.749957670e2,-.268938672e-5,0.444896105,
       0.463434311e-2,-.178591076e-5,-.268938672e-5,
       0,-1.344541576008661e-3,-5.062194923138079e3,
       -6.833142732779555e-5,1.449382650173157e-8,
       -4.268463211410861,2.098334687947376e-1,
       -6.397251492537153e-08,3.824589508329281e2,
       -4.376060460948886e-09)
names(y) <- c("phi1","phi2","x3","q1","q2","q3","q4",
              "vphi1","vphi2","vx3","vq1","vq2","vq3","vq4",
              "aphi1","aphi2","ax3","aq1","aq2","aq3","aq4",
              "la1","la2","la3")
yprime <- c(y[8:21],rep(0,10))
crank2 <- mebdfi(y=y,dy=yprime,times=times,res="crankres", nind=nind,
          dllname="crank",initfunc=NULL, parms=NULL, ipar=ipar,
          hini=Tol*1e-2,atol=Tol,rtol=Tol,maxsteps=100000)
plot(crank2, type="l", lwd=2,
  which=c("phi2","x3","q1","q2","q3","q4","la1","la2","la3"))


ipar <- c(1,1)
crank3 <- mebdfi(y=y,dy=yprime,times=times,res="crankres", nind=nind,
          dllname="crank",initfunc=NULL, parms=NULL, ipar=ipar,
          hini=Tol*1e-2,atol=Tol,rtol=Tol,maxsteps=100000)
plot(crank3, type="l", lwd=2,
  which=c("phi2","x3","q1","q2","q3","q4","la1","la2","la3"))
