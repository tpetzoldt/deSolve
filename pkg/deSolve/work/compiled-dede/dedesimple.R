### Simple DDE, adapted version of ?dede example from package deSolve

library(deSolve)

derivs <- function(t, y, parms) {
  with(as.list(parms), {
    if (t < tau)
      ytau <- 1
    else
      ytau <- lagvalue(t - tau)

    dy <- k * ytau
    list(c(dy))
  })
}

yinit <- 1
times <- seq(0, 30, 0.1)
parms <- c(tau = 1, k = -1)

#yout <- dede(y = yinit, times = times, func = derivs, parms = parms)

#plot(yout, main = "dy/dt = -y(t-1)")


system("R CMD SHLIB dedesimple.c")
dyn.load("dedesimple.dll")

yout2 <- dede(yinit, times = times, func = "derivs", parms = parms,
  dllname = "dedesimple", initfunc = "initmod", nout = 1)

dyn.unload("dedesimple.dll")


plot(yout2, main=c("y", "ytau"))
