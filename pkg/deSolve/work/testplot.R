require(deSolve)

LVmatrix <- function(t, n, parms) {
  with(parms, {
    dn <- r * n + n * (A %*% n)
    return(list(dn, sum = sum(n)))
  })
}
parms <- list(
  r = c(r1 = 0.1, r2 = 0.1, r3 = -0.1, r4 = -0.1),
  A = matrix(c(0.0, 0.0, -0.2, 0.01,      # prey 1
               0.0, 0.0, 0.02, -0.1,      # prey 2
               0.2, 0.02, 0.0, 0.0,       # predator 1; prefers prey 1
               0.01, 0.1, 0.0, 0.0),      # predator 2; prefers prey 2
               nrow = 4, ncol = 4, byrow=TRUE)
)
times <- seq(from = 0, to = 500, by = 1)
y     <- c(prey1 = 1, prey2 = 1, pred1 = 2, pred2 = 2)

out <- ode(y, times, LVmatrix, parms)

parms2 <- parms
parms2$r[1] <- 0.2
out2 <- ode(y, times, LVmatrix, parms2)

parms3 <- parms
parms3$r[2] <- 0.2
out3 <- ode(y, times, LVmatrix, parms3)




## Basic line plot
plot(out, type = "l")

## User-specified axis labels
plot(out, which = 1:4, type = "l", ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", main = "Time Series")

## add second output 
plot(out, out2, type = "l", ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2", "Sum"),
  xlab = "Time (d)", main = "Time Series", col = c("red", "blue"))

## add second output 
plot(out, out2, which = 1:4, type = "l", ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", main = "Time Series", col = c("red", "blue"))


## add third output
plot(out, out2, out3, type = "l", which = 1:4, 
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= 1:4, main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1)
  

 
## ylim is vector
plot(out, out2, out3, which = 1:5, type = "l",
  main = c("Prey 1", "Prey 2", "Pred 1", "Pred 2", "Sum"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= 1,  lwd=1,
  ylim = c(0,20))

  
## ylim is list
plot(out, out2, out3, which = 1:4,  type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim=list(c(0,20)))
  
## individual ylim
plot(out, out2, out3, which = 1:4,  type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim=list(c(0,15), c(0, 12), c(0,20), c(0, 25)))
  
## individual ylim, some not specified   - THIS  ....DOES NOT WORK...
plot(out, out2, out3, which = 1:4,  type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim=list(c(0,15), c(0, 12), NULL, NULL))


## incomplete ylim list
plot(out, out2, out3, which = 1:4,  type = "l",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), lwd=1,
  ylim=list(c(0,15), c(0, 12), c(0,20)))
  
## type"p"
plot(out, out2, out3, which = 1:4,  type = "p",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), cex=0.5)

plot(out, out2, out3, which = 1:4,  type = "p",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), 
  cex=0.8, cex.axis=.5, cex.lab=1.5, cex.main=1)

plot(out, out2, out3, which = 1:4,  type = "p",
  ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", col = c("red", "blue", "forestgreen"),
  lty= c(1,1,1), main=c("Prey 1", "Prey 2", "Prey 3", "Prey 4"), 
  cex=c(1, 0.5, 0.3), cex.axis=c(0.5, 0.2, 0.8), cex.lab=c(0.5, 0.2, 0.8), 
  cex.main=c(0.5, 0.2, 0.8))

  
## ThPe: todo
## - list version of cex  
## KS: WHY????
## ThPe: to encode additional information; done (experimentally), see example:

obs <- as.data.frame(out[out[,1] %in% seq(10, 500, by=30), ])

nobs <- nrow(obs)
rbd <- rainbow(nobs)
plot(out, type = "l", obs = obs,
  obspar = list(col=list(rbd, "red", "blue", rbd),
  pch=list(18, 19, 1:2, 1:nobs), cex=list(2, seq(0.5, 3, length=nobs), 1.2, 1.2)))


################################################################################
  
## add second output 
plot(out, out2, which = c("prey1", "prey2"), 
type = "l", ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", main = "Time Series", col = c("red", "blue"))

plot(out, out2, which = "prey2",
type = "l", ylab = c("Prey 2", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", main = "Time Series", col = c("red", "blue"))

timeobs <- seq (0,300,50)
prey2   <- (1 - 0.7*sin(2 * pi * timeobs / 60))*(1+0.1*runif(length(timeobs)))
obsdat <- cbind (time=timeobs, prey2 = prey2)

obsdat[3,2] <- NA
obsdat <- cbind (obsdat, sum = NA)
obsdat <- rbind (obsdat, c(100,NA, 2), c(300,NA, 3.5))

# plotting withe different parameter settings
plot(out, out2, type = "l", xlab = "Time (d)", col = c("red", "blue"),
  obs = obsdat, obspar = list(pch = 16, cex = 2, col = c("orange", "darkblue")))

plot(out, out2, which = 5:1, type = "l", xlab = "Time (d)", col = c("red", "blue"),
  obs = obsdat, obspar = list(pch = 16, cex = 2, col = c("orange", "darkblue")))


### ============================================================================
### The histogram function
### ============================================================================
hist(out, col = c("blue", "red", "green"))

hist(out, col = c("blue", "red", "green"), 
  xlim = c(0,6))

hist(out, col = c("blue", "red", "green"), 
  xlim = list(c(0.,2.0),c(0.4,2),c(0.4,2),c(0.4,2),c(0,6)))

# This does also work....
hist(out, col = c("blue", "red", "green"), 
  xlim = list(c(0.,2),c(0.4,2),NULL,NULL,NULL))

###
### images and filled.contours
###


## ================
## Model equations
## ================
require (ReacTran)

lvmod <- function (time, state, parms, N, rr, ri, dr, dri) {
  with (as.list(parms), {
    PREY <- state[1:N]
    PRED <- state[(N+1):(2*N)]

    ## transport
    tranPrey <- tran.1D(C = PREY, D = Da, dx = dri)
    tranPred <- tran.1D(C = PRED, D = Da, dx = dri)

    ## Biology: Lotka-Volterra model
    Ingestion     <- rIng  * PREY * PRED
    GrowthPrey    <- rGrow * PREY * (1-PREY/cap)
    MortPredator  <- rMort * PRED

    ## Rate of change = Flux gradient + Biology
    dPREY    <- tranPrey$dC + GrowthPrey - Ingestion
    dPRED    <- tranPred$dC + Ingestion * assEff - MortPredator

    return (list(c(dPREY, dPRED)))
  })
}

## ==================
## Model application
## ==================

## model parameters:

R  <- 20                        # total radius of surface, m
N  <- 100                       # 100 concentric circles
dr <- R/N                       # thickness of each layer
r  <- seq(dr/2,by = dr,len = N) # distance of center to mid-layer
ri <- seq(0,by = dr,len = N+1)  # distance to layer interface
dri <- dr                       # dispersion distances

parms <- c(Da     = 0.05,       # m2/d, dispersion coefficient
           rIng   = 0.2,        # /day, rate of ingestion
           rGrow  = 1.0,        # /day, growth rate of prey
           rMort  = 0.2 ,       # /day, mortality rate of pred
           assEff = 0.5,        # -, assimilation efficiency
           cap    = 10)         # density, carrying capacity

## Initial conditions: both present in central circle (box 1) only
state    <- rep(0, 2 * N)
state[1] <- state[N + 1] <- 10
                
## RUNNING the model:
times  <- seq(0, 200, by = 1)   # output wanted at these time intervals

## the model is solved by the two implemented methods:
## 1. Default: banded reformulation
print(system.time(
  out <- ode.1D(y = state, times = times, func = lvmod, parms = parms,
                nspec = 2, names = c("PREY", "PRED"),
                N = N, rr = r, ri = ri, dr = dr, dri = dri)
))

## 2. Using sparse method
print(system.time(
  out2 <- ode.1D(y = state, times = times, func = lvmod, parms = parms,
                 nspec = 2, names = c("PREY","PRED"), 
                 N = N, rr = r, ri = ri, dr = dr, dri = dri,
                 method = "lsodes")
))

## ================
## Plotting output
## ================
# the data in 'out' consist of: 1st col times, 2-N+1: the prey
# N+2:2*N+1: predators

image(out, which = "PREY", grid = r, xlab = "time, days", 
      ylab = "Distance, m", main = "Prey density")

# zoom in
image(out, which = c("PREY","PREY"), grid = r, xlab = "time, days", 
      ylab = "Distance, m", main = "Prey density",
      xlim = list(NULL,c(0,20)), ylim = list(NULL, c(0,5)))



