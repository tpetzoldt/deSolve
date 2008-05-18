oask   <- devAskNewPage(dev.interactive(orNone = TRUE))
oldpar <- par(mfrow=c(2,2))

## a helper function for plotting
plotIt <- function(col="red") {

  plot (out1$time, out1$s,  type="l",   ylim=c(0,3))
  lines(out2$time, out2$s, col=col,   lty="dotted", lwd=3)
  lines(out3$time, out3$s, col="green", lty="dotted")

  plot (out1$time, out1$p, type="l",    ylim=c(0,3))
  lines(out2$time, out2$p, col=col,   lty="dotted", lwd=3)
  lines(out3$time, out3$p, col="green", lty="dotted")

  plot (out1$time, out1$k, type="l",    ylim=c(0,3))
  lines(out2$time, out2$k, col=col,   lty="dotted", lwd=3)
  lines(out3$time, out3$k, col="green", lty="dotted")

  plot (out1$p, out1$k, type="l")
  lines(out2$p, out2$k, col=col,   lty="dotted", lwd=3)
  lines(out3$p, out3$k, col="green", lty="dotted")
}


## A simple resource limited Lotka-Volterra-Model
lvmodel <- function(t, x, parms) {
  s <- x[1] # substrate
  p <- x[2] # producer
  k <- x[3] # consumer
  with(as.list(parms),{
    import <- approx(signal$times, signal$import, t, rule=2)$y
    ds <- import - b*s*p + g*k
    dp <- c*s*p  - d*k*p
    dk <- e*p*k  - f*k
    res<-c(ds, dp, dk)
    list(res)
  })
}

## vector of timesteps
times  <- seq(0, 100, length=101)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = 0:200,
                           import = rep(0,length(0:200))))

signal$import[signal$times >= 10 & signal$times <=11] <- 0.2

## Parameters for steady state conditions
parms  <- c(a=0.0, b=0.0, c=0.1, d=0.1, e=0.1, f=0.1, g=0.0)

## Start values for steady state
y<-xstart <- c(s=1, p=1, k=1)



## Classical RK4 with fixed time step
system.time(out1  <- as.data.frame(rk4(xstart, times, lvmodel, parms)))

## LSODA
system.time(out3 <- as.data.frame(lsoda(xstart, times, lvmodel,
  hmax=1, parms)))

## same: rk4 fixed step, generalized implementation
system.time(out2 <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=1,
  method = rkMethod("rk4")))) 
  
plotIt("blue")
  
## rk2, note smaller time step !! 
system.time(out2x <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=.5,
  method = rkMethod("rk2"))))

plotIt()

## Euler, note smaller time step !! 
system.time(out2x <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=.1,
  method = rkMethod("euler"))))

plotIt("blue")

## rk23bs
system.time(out2 <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=1,
  method = rkMethod("rk23bs"))))

plotIt()
 
## Runge-Kutta-Fehlberg  45
system.time(out2 <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=1,
  hmax = 1, method = rkMethod("rk45f"))))

plotIt("blue")

## Prince-Dormand  5(4)7m
system.time(out2x <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=1,
  hmax = 1, hmin=0.1, method = rkMethod("rk45dp7"))))

plotIt()
  
##other syntax 
system.time(out2x <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=1,
  hmax = 1, hmin=0.1, method = "rk45dp7")))  

plotIt("blue")

## tolerance values can also be set for single state variables
## if all tolerance values are zero, hini is retained as fixed step
## Prince-Dormand  5(4)7m
system.time(out2x <- as.data.frame(rk(xstart, times, lvmodel, parms, hini=1,
  hmax = 10, method = rkMethod("rk45dp7"), atol=c(0, 0, 0), rtol=c(0, 0, 0))))

plotIt()

## clean up
par(oldpar)
devAskNewPage(oask)




