## THIS DOES NOT WORK ... yet

## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        NAND gate
##        index 0 IDE of dimension 14
##
##
##     The implementation is as a DLL; code in "nandgate.f"
##
#---------------------------------------------------------------------------
# before trying this code, the fortran programme has to be compiled
# this can be done in R:
# dyn.unload("nandgate.dll");system("R CMD SHLIB nandgate.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

require(daeSolve)
dyn.load("nandgate.dll")

#-----------------------------------------------------------------------
# initialising
VBB    <-  -2.5

parms <- c(RGS = 4, RGD = 4, RBS = 10, RBD = 10,
      CGS = .6e-4, CGD = .6e-4, CBD = 2.4e-5, CBS = 2.4e-5,
      C9  = .5e-4, DELTA = 0.2e-1, CURIS = 1.e-14, VTH = 25.85,
      VDD = 5., VBB = VBB)

yini   <- c(5,5,VBB,VBB,5,3.62385,5,VBB,VBB,3.62385,0,3.62385,VBB,VBB)
yprime <- rep(0,14)

# index of the system
DLLres(res="gateres",time=0.,y=yini,dy=yprime,parms= parms,
  initfunc="initgate", dllname="nandgate")

times <-seq(0,80,by=1)         # time: from 0 to 80 hours, steps of 1 hour

# integrate the model: low tolerances to restrict integration time
out<- mebdfi(y=yini, dy=yprime, times, res="gateres", initfunc="initgate",
             parms=parms, dllname="nandgate", hini=1e-6, rtol=1e-6, atol=1e-6)

plot(out,type="l")
