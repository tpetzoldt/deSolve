## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Andrews'' squeezing mechanism (in index 3 formulation)
##        index 3 DAE of dimension 27
##
##     The implementation is as a DLL; code in "andrews.f"
##
## =============================================================================

#---------------------------------------------------------------------------
# before trying this code, the fortran programme has to be compiled
# this can be done in R:
# dyn.unload("andrews.dll");system("R CMD SHLIB andrews.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

require(daeSolve)
dyn.load("andrews.dll")

yini <- c(-0.0617138900142764496358948458001, 0,
           0.455279819163070380255912382449, 0.222668390165885884674473185609,
           0.487364979543842550225598953530,-0.222668390165885884674473185609,
           1.23054744454982119249735015568 ,0,     0,0,     0,0,          0,0,
           14222.4439199541138705911625887,-10666.8329399655854029433719415,
           0,0,   0,0,  0,98.5668703962410896057654982170,
           -6.12268834425566265503114393122,0,         0,0,   0)

yprime <- rep(0,27)
yprime[1:14] <- yini[8:21]

parameter <- c(m1=.04325, m2=.00365, m3=.02373 ,m4=.00706 ,
                m5=.07050 ,m6=.00706 ,m7=.05498 ,
                xa=-.06934 ,ya=-.00227 ,
                xb=-0.03635 ,yb=.03273 ,
                xc=.014 ,yc=.072 ,c0=4530 ,
                i1=2.194e-6,i2=4.410e-7,i3=5.255e-6,i4=5.667e-7,
                i5=1.169e-5,i6=5.667e-7,i7=1.912e-5,
                d=28e-3,da=115e-4,e=2e-2,ea=1421e-5,
                rr=7e-3,ra=92e-5,l0=7785e-5,
                ss=35e-3,sa=1874e-5,sb=1043e-5,sc=18e-3,sd=2e-2,
                ta=2308e-5,tb=916e-5,u=4e-2,ua=1228e-5,ub=449e-5,
                zf=2e-2,zt=4e-2,fa=1421e-5,mom=33e-3)

DLLres(res="andres", time=0., y=yini, dy=yprime, parms=parameter,
         initfunc="andinit", dllname="andrews")

times <- seq(0,0.03,by=0.001)
ind <- c(7,7,13)

# including the jacobian
print(system.time(
AndOut <- mebdfi(y=yini, dy=yprime, times=times, res="andres", nind=ind,
          dllname="andrews", jacres="andjac", initfunc="andinit",
          parms=parameter, hini=0.01, atol=1e-10, rtol=1e-10,
          jactype="fullusr", maxsteps=100000)
))

# numerically differenced jacobian
print(system.time(
AndOut2 <- mebdfi(y=yini, dy=yprime, times=times, res="andres", nind=ind,
          dllname="andrews", parms=parameter, initfunc = "andinit",
          hini=0.01, atol=1e-10, rtol=1e-10,
          jactype="fullint", maxsteps=100000)
))

plot(AndOut,which=1:9,type="l",lwd=2)
