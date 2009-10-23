###############################################
## The AQUAPHY model                         ##
##                                           ##
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 2. Model formulation              ##
## Equations discussed in Chapter 2.9.2      ##
###############################################


# load package with the integration routine:

require(mebdfi)

#----------------------#
# the model equations: #
#----------------------#

model<-function(t,state,parameters)
  {
  with(as.list(c(state,parameters)),{  # unpack the state variables, parameters

    # PAR, on-off function depending on the hour within a day
    hourofday       <- t%%24
    PAR <- ifelse (hourofday  < dayLength, parMean , 0)

    ## the output variables
    PhytoC           <- PROTEIN + RESERVE + LMW       # all components contain carbon
    PhytoN           <- PROTEIN * rNCProtein          # only proteins contain nitrogen
    NCratio          <- PhytoN / PhytoC           
    Chlorophyll      <- PhytoN * rChlN
    TotalN           <- PhytoN + DIN
    ChlCratio        <- Chlorophyll / PhytoC

    ## the rates, in mmol/hr 
    PartLMW          <- LMW / PhytoC
    Limfac           <- max(0,min(1,(maxpLMW -PartLMW)/(maxpLMW-minpLMW)))
    PhotoSynthesis   <- maxPhotoSynt*Limfac*(1-exp(alpha*PAR/maxPhotoSynt)) * PROTEIN
    Exudation        <- pExudation * PhotoSynthesis 
    MonodQuotum      <- max(0,LMW / PROTEIN - minQuotum)
    ProteinSynthesis <- maxProteinSynt*MonodQuotum * DIN / (DIN+ksDIN)      * PROTEIN
    Storage          <- maxStorage    *MonodQuotum                          * PROTEIN
    Respiration      <- respirationRate * LMW + pResp*ProteinSynthesis 
    Catabolism       <- catabolismRate  * RESERVE

    ## the rates of change of state variables; includes dilution effects (last term)
    dLMW     <- ( PhotoSynthesis + Catabolism
                - Exudation - Storage  - Respiration - ProteinSynthesis 
                - dilutionRate * LMW)

    dRESERVE <-  Storage - Catabolism          - dilutionRate * RESERVE

    dPROTEIN <-  ProteinSynthesis              - dilutionRate * PROTEIN

    dDIN     <- -ProteinSynthesis * rNCProtein - dilutionRate * (DIN - inputDIN)


    ## the output, as a list
    list(c(dDIN,dPROTEIN,dRESERVE,dLMW),              ## the rate of change of state variables
           c(PAR               = PAR,                 ## the ordinary variables
             TotalN            = TotalN,
             PhotoSynthesis    = PhotoSynthesis,
             NCratio           = NCratio,
             ChlCratio         = ChlCratio,
             Chlorophyll       = Chlorophyll))
    })
 }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#

parameters<-c(maxPhotoSynt   =0.125,      #molC/molC/hr      Maximal protein C-specific rate of photsynthesis at 20 dg
              rMortPHY       =0.001,      #/hr               Mortality rate of Phytoplankton (lysis and zooplankton grazing)
              alpha          =-0.125/150, #µEinst/m2/s/hr    Light dependency factor
              pExudation     =0.0,        #-                 Part of photosynthesis that is exudated
              maxProteinSynt =0.136,      #molC/molC/hr      Maximal Biosynthetic C-specific N-uptake rate 
              ksDIN          =1.0,        #mmolN/m3          Half-saturation ct of N uptake Phytoplankton        
              minpLMW        =0.05,       #molC/molC         Minimum metabolite/totalC ratio in algae
              maxpLMW        =0.15,       #molC/molC         Maximum metabolite/totalC ratio in algae
              minQuotum      =0.075,      #molC/molC         Minimum metabolite/Protein ratio for synthesis
              maxStorage     =0.23,       #/h                Maximum storage rate for Phytoplankton                             
              respirationRate=0.0001,     #/h                Respiration rate of LMW
              pResp          =0.4,        #-                 Part of protein synthesis that is respired (cost of biosynthesis)
              catabolismRate =0.06,       #/h                Catabolism rate of Phytoplankton reserves 
              dilutionRate   =0.01,       #/h                dilution rate in chemostat
              rNCProtein     =0.2,        #molN/molC         Nitrogen/carbon ratio of proteins
              inputDIN       =10.0,       #mmolN/m3          DIN in inflowing water
              rChlN          =1,          #gChl/molN         Chl to nitrogen ratio
              parMean        =250.,       #µmolPhot/m2/s     PAR during the light phase
              dayLength      =15.         #hours             Length of illuminated period 
              )

#-------------------------#
# the initial conditions: #
#-------------------------#
 
# assume the amount of reserves = 50% amount of proteins
# 10% LMW

state     <-c(DIN     =6.,     #mmolN/m3
              PROTEIN =20.0,   #mmolC/m3
              RESERVE =5.0,    #mmolC/m3
              LMW     =1.0)    #mmolC/m3

#----------------------#
# RUNNING the model:   #
#----------------------#

times <-seq(0,24*10,1)
out   <-as.data.frame(mebdfi(state,times,model,parameters))


#------------------------#
# PLOTTING model output: #
#------------------------#

par(mfrow=c(2,2), oma=c(0,0,3,0))         # set number of plots (mfrow) and margin size (oma)
col <- grey(0.9)
ii <- 1:length(out$PAR)                   # output over entire period

plot (times[ii],out$Chlorophyll[ii],type="l",main="Chlorophyll",xlab="time, hours",ylab="µg/l")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$Chlorophyll[ii]  ,lwd=2 )


plot (times[ii],out$DIN[ii]        ,type="l",main="DIN"        ,xlab="time, hours",ylab="mmolN/m3")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$DIN[ii]  ,lwd=2 )


plot (times[ii],out$NCratio[ii]    ,type="n",main="NCratio"    ,xlab="time, hours",ylab="molN/molC")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$NCratio[ii]  ,lwd=2 )


plot (times[ii],out$PhotoSynthesis[ii],type="l",main="PhotoSynthesis",xlab="time, hours",ylab="mmolC/m3/hr")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$PhotoSynthesis[ii]  ,lwd=2 )

mtext(outer=TRUE,side=3,"AQUAPHY",cex=1.5)

#------------------------#
# SUMMARY  model output: #
#------------------------#
t(summary(out))
