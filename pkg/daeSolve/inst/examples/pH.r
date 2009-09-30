###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 8.                                ##
## Multiple times scales and equilibrium     ##
## Chapter 8.5.1, 8.5.2                      ##
## Solving pH in aquatic systems             ##
## Model of pH changes due to algal growth   ##
###############################################

# load package with the integration routine:
require(mebdfi)

#########################################################################
##                                                                     ##
## FUNCTIONS FOR pH CALCULATIONS                                       ##
##    Karline Soetaert                                                 ##
##                                                                     ##
#########################################################################

require(seacarb)

Salinity     <- 0
Temperature  <- 20
WDepth       <- 0

# the dissociation constants
k1           <- K1(Salinity,Temperature,WDepth)    # Carbonate k1
k2           <- K2(Salinity,Temperature,WDepth)    # Carbonate k2


#========================================================================
# solve for equilibrium pH for given alkalinity, DIC concentration.
# Uniroot is root solver routine (for one function) in R
#========================================================================

pHfunction <- function(pH, k1,k2, DIC, Alkalinity )
{
   H    <- 10^(-pH)
   HCO3 <- H*k1  /(H*(k1+H) + k1*k2)*DIC
   CO3  <- k1*k2 /(H*(k1+H) + k1*k2)*DIC

   EstimatedAlk  <- (- H) *1.e6  + HCO3 + 2*CO3

   return(EstimatedAlk  - Alkalinity)

}

Alkalinity <- 2200
DIC        <- 2100

sol <- uniroot(pHfunction,lower=0,upper=12, tol=1.e-20,
                k1=k1, k2=k2, DIC=DIC, Alkalinity=Alkalinity)
sol$root


#========================================================================
# pH changes due to primary production
#========================================================================

#----------------------#
# the model equations: #
#----------------------#
model<-function(time,state,parameters)
  {
with(as.list(c(state,parameters)),{

    PAR    <- 0.
    if(time%%24 < dayLength) PAR <- parDay

    Growth <- maxGrowth*DIN/(DIN+ksDIN)*PAR/(PAR+ksPAR)*ALGAE -
                respRate * ALGAE

    dDIN   <- -Growth                   # DIN is consumed
    dDIC   <- -Growth * CNratio         # DIC is consumed ~ CN ratio
    dALGAE <- Growth                    # algae increase by growth
    dALKALINITY <- Growth               #alkalinity production if nitrate


    # estimate the pH
    pH  <- uniroot(pHfunction,lower=0,upper=12,tol=1.e-20,
           k1=k1,k2=k2, DIC=DIC,Alkalinity=ALKALINITY)$root

   list(c(dDIN,dALGAE,dALKALINITY,dDIC),c(PAR=PAR,pH=pH)  )

    })
 }

#-----------------------#
# the model parameters: #
#-----------------------#

parameters<-c(maxGrowth   =0.125,      #molN/molN/hr  Maximal growth rate
              ksPAR       =100,        #µEinst/m2/s   Half-saturation ct for light-limited growth
              ksDIN       =1.0,        #mmolN/m3      Half-saturation ct of N uptake Phytoplankton
              respRate    =0.001,      #/h            Respiration rate
              CNratio     =6.5,        #molC/molN     carbon:Nitrogen ratio
              parDay      =250.,       #µEinst/m2/s   PAR during the light phase
              dayLength   =12.         #hours         Length of illuminated period (in one day)
              )
Salinity     <- 0
Temperature  <- 20
WDepth       <- 0

# the dissociation constants
k1           <- K1(Salinity,Temperature,WDepth)    # Carbonate k1
k2           <- K2(Salinity,Temperature,WDepth)    # Carbonate k2

#-------------------------#
# the initial conditions: #
#-------------------------#

state     <-c(DIN        =30,     #mmolN/m3
              ALGAE      =0.1,    #mmolN/m3
              ALKALINITY =2200,   #mmol/m3
              DIC        =2100)   #mmolC/m3

#----------------------#
# RUNNING the model:   #
#----------------------#

times <-seq(0,24*10,1)

require(deSolve)
out   <-as.data.frame(mebdfi(state,times,model,parameters))


#------------------------#
# PLOTTING model output: #
#------------------------#

par(mfrow=c(2,2), oma=c(0,0,3,0))         # set number of plots (mfrow) and margin size (oma)
plot (out$time,out$ALGAE,type="l",main="Algae",xlab="time, hours",ylab="mmol/m3")
polygon(out$time,out$PAR-10,col="lightgrey",border=NA)
box()
lines (out$time,out$ALGAE  ,lwd=2 )

plot (out$time,out$DIN ,type="l",main="DIN"  ,xlab="time, hours",ylab="mmolN/m3", lwd=2)
plot (out$time,out$DIC ,type="l",main="DIC"  ,xlab="time, hours",ylab="mmolC/m3", lwd=2)
plot (out$time,out$pH,type="l",main="pH"  ,xlab="time, hours",ylab="-", lwd=2)
mtext(outer=TRUE,side=3,"Algal growth and pH",cex=1.5)
