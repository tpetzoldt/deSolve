

parms  <- c(a=0.0, b=0.0, c=0.1, d=0.1, e=0.1, f=0.1, g=0.0)

## The model
#lvmodel <- function(t, x, parms, input)  {
lvmodel <- function(t, x, parms)  {
    with(as.list(c(parms,x)),  {
      import <- sigimp(t) ## global variable
      dS <- import - b*S*P + g*K
      dP <- c*S*P  - d*K*P
      dK <- e*P*K  - f*K
      res<-c(dS, dP, dK)
      list(c(res), c(a=a, dS=dS))
    })
  }

## vector of timesteps
times  <- seq(0, 100, length=1001)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = times,
                            import = rep(0,length(times))))

signal$import[signal$times >= 10 & signal$times <=11] <- 0.2
sigimp <- approxfun(signal$times, signal$import, rule=2)

## Start values for steady state
xstart <- c(S=1, P=1, K=1)


Func    <- function(time,state,parms) {
             attr(state,"names") <- Ynames
             lvmodel(time, state, parms)[1]
}

Ynames <- names(xstart)

rho <- environment(lvmodel)
# SEXP call_rkAuto(SEXP xstart, SEXP times, SEXP func, SEXP parms, SEXP rho,
#     SEXP rtol, SEXP atol, SEXP tcrit, //SEXP verbose,
#     SEXP hmin, SEXP hmax, SEXP hini, SEXP method, SEXP maxsteps)
