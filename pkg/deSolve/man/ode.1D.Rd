\name{ode.1D}
\alias{ode.1D}
\title{Steady-state solver for multicomponent 1-D ordinary differential equations}
\description{Solves the steady-state of a system of ordinary differential equations. 
    \cr Assumes a banded or sparse jacobian matrix. To create the banded matrix, 
    state variables are rearranged. 
    \cr Suitable for 1-D models that include transport only between adjacent layers and that model many species (e.g. multi-component reactive-transport models).}
\usage{ode.1D(y, times, func, nspec, method="lsode", ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{time sequence for which output is wanted; the first value of \code{times} must be the initial time}
  \item{func }{either a user-supplied function that computes the values of the derivatives in the ODE system (the model definition) at time t, or a character string giving the name of a compiled function in a dynamically loaded shared library. 
       		If \code{func} is a character string then integrator \code{lsodes} will be used. See details}
  \item{nspec }{the number of *species* (components) in the model. }
  \item{method }{the integrator to use, one of "vode", "lsode", "lsoda", "lsodar", "lsodes" }
  \item{... }{additional arguments passed to the integrator}
}
\value{
  A matrix with up to as many rows as elements in times and as many columns as elements in \code{y} plus the number of "global" values returned 
      in the second element of the return from \code{func}, plus an additional column (the first) for the time value. 
      There will be one row for each element in \code{times} unless the integrator returns with an unrecoverable error. 
      If \code{y} has a names attribute, it will be used to label the columns of the output value. 
      
      The output will have the attributes \code{istate}, and \code{rstate}, two vectors with several useful elements. 
      The first element of istate returns the conditions under which the last call to the integrator returned. Normal is istate = 2.
      If \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen. See the help for the selected integrator for details.
  
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# example 1

  #=======================================================
  # a predator and its prey diffusing on a flat surface
  # in concentric circles
  # 1-D model with using cylindrical coordinates
  # Lotka-Volterra type biology
  #=======================================================

  #==================#
  # Model equations  #
  #==================#

  lvmod <- function (time, state, parms,N,rr,ri,dr,dri)

  {
    with (as.list(parms),{
    PREY <- state[1:N]
    PRED <- state[(N+1):(2*N)]
    
    # Fluxes due to diffusion 
    # at internal and external boundaries: zero gradient
    FluxPrey <- -Da * diff(c(PREY[1],PREY,PREY[N]))/dri   
    FluxPred <- -Da * diff(c(PRED[1],PRED,PRED[N]))/dri   

    # Biology: Lotka-Volterra model
    Ingestion     <- rIng * PREY*PRED
    GrowthPrey    <- rGrow* PREY*(1-PREY/cap)
    MortPredator  <- rMort* PRED

    # Rate of change = Flux gradient + Biology   
    dPREY    <- -diff(ri * FluxPrey)/rr/dr   +
                GrowthPrey - Ingestion
    dPRED    <- -diff(ri * FluxPred)/rr/dr   +
                Ingestion*assEff -MortPredator

    return (list(c(dPREY,dPRED)))
   }) 
  }
  
  #==================#
  # Model application#
  #==================#
  # model parameters: 

  R  <- 20                    # total radius of surface, m
  N  <- 100                   # 100 concentric circles
  dr <- R/N                   # thickness of each layer
  r  <- seq(dr/2,by=dr,len=N) # distance of center to mid-layer
  ri <- seq(0,by=dr,len=N+1)  # distance to layer interface
  dri<- dr                    # dispersion distances

  parms <- c( Da     =0.05,   # m2/d, dispersion coefficient 
              rIng   =0.2,    # /day, rate of ingestion
              rGrow  =1.0,    # /day, growth rate of prey
              rMort  =0.2 ,   # /day, mortality rate of pred
              assEff =0.5,    # -, assimilation efficiency
              cap    =10  )   # density, carrying capacity

  # Initial conditions: both present in central circle (box 1) only
  state <- rep(0,2*N)
  state[1] <- state[N+1] <- 10
                  
  # RUNNING the model:   #
  times  <-seq(0,200,by=1)   # output wanted at these time intervals           

  out    <- ode.1D(y=state,times=times,func=lvmod,parms=parms,nspec=2,
                    N=N,rr=r,ri=ri,dr=dr,dri=dri)  

  #==================#
  # Plotting output  #
  #==================#
  # the data in 'out' consist of: 1st col times, 2-N+1: the prey
  # N+2:2*N+1: predators

  PREY   <- out[,2:(N  +1)]

  filled.contour(x=times,y=r,PREY,color= topo.colors,
                 xlab="time, days", ylab= "Distance, m",
                 main="Prey density")
 
  # Example 2.
  #=======================================================
  # Biochemical Oxygen Demand (BOD) and oxygen (O2) dynamics
  # in a river
  #=======================================================
  
  #==================#
  # Model equations  #
  #==================#
  O2BOD <- function(t,state,pars)
  
  {
    BOD <- state[1:N]
    O2  <- state[(N+1):(2*N)]
  
  # BOD dynamics
    FluxBOD <-  v*c(BOD_0,BOD)  # fluxes due to water transport
    FluxO2  <-  v*c(O2_0,O2)
    
    BODrate <- r*BOD            # 1-st order consumption
  
  #rate of change = flux gradient - consumption  + reaeration (O2)
    dBOD         <- -diff(FluxBOD)/dx  - BODrate
    dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)
  
    return(list(c(dBOD=dBOD,dO2=dO2)))
  
  }    # END O2BOD
   
   
  #==================#
  # Model application#
  #==================#
  # parameters
  dx      <- 25        # grid size of 25 meters
  v       <- 1e3       # velocity, m/day
  x       <- seq(dx/2,5000,by=dx)  # m, distance from river
  N       <- length(x)
  r       <- 0.05      # /day, first-order decay of BOD
  p       <- 0.5       # /day, air-sea exchange rate 
  O2sat   <- 300       # mmol/m3 saturated oxygen conc
  O2_0    <- 200       # mmol/m3 riverine oxygen conc
  BOD_0   <- 1000      # mmol/m3 riverine BOD concentration
  
  # initial conditions:
  state <- c(rep(200,N),rep(200,N))
  times     <- seq(0,20,by=1)
  
  # running the model
  #  step 1  : model spinup
  out       <- ode.1D (y=state,times,O2BOD,parms=NULL,nspec=2)
  
  #==================#
  # Plotting output  #
  #==================#
  # select oxygen (first column of out:time, then BOD, then O2
  O2   <- out[,(N+2):(2*N+1)]
  color= topo.colors
  
  filled.contour(x=times,y=x,O2,color= color,nlevels=50,
                 xlab="time, days", ylab= "Distance from river, m",
                 main="Oxygen")}
\details{This is the method of choice for multi-species 1-dimensional models, that are only subjected to transport between adjacent layers 
    \cr More specifically, this method is to be used if the state variables are arranged per species: 
    \cr A[1],A[2],A[3],....B[1],B[2],B[3],.... (for species A, B)) 
  
  
    For single-species 1-D models, use \code{\link{ode.band}}. 
  
  If state variables are arranged as (e.g. A[1],B[1],A[2],B[2],A[3],B[3],.... . 
  it is generally more efficient to use \code{\link{ode.band}} 
  
  If the integrator chosen is not \code{lsodes} then a banded Jacobian is created by rearranging state variables.
  \cr If the integrator chosen is \code{lsodes} then state variables are not rearranged, and the sparse Jacobian is selected.
  \cr If the model is specified in a DLL, then \code{lsodes} is used.

  See the selected integrator for the additional options}
\seealso{\item \code{\link{ode.band}}, for integrating, when the jacobian matrix is banded, and where the state variables need NOT be rearranged 
  	\item \code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}}.}
\keyword{math}

