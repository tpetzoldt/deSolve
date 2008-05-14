
### ode.1D, ode.band -- solves banded ordinary differential equation systems 
### ode.1D is designed for solving multi-component 1-D reaction-transport models 
### ode.band is designed for solving single-component 1-D reaction-transport models
### they offer the choice between the integrators vode, lsode, lsoda, lsodar and lsodes.

ode    <- function (y,
                    times,
                    func,
                    parms,
                    method= c("lsoda","lsode","lsodes","lsodar","vode","daspk"),
                    ...)
{
  if (is.null(method)) method <- "lsoda"
  if ( is.list(method)) 
    out <- NULL # Thomas: 
  else if (is.function(method))
    out <- method(y,times,func,parms,...)
  else
   out <- switch(match.arg(method),
    lsoda = lsoda(y,times,func,parms,...),
    vode  = vode(y,times,func,parms,...),
    lsode = lsode(y,times,func,parms,...),
    lsodes=lsodes(y,times,func,parms,...),
    lsodar=lsodar(y,times,func,parms,...),
    daspk = daspk(y,times,func,parms,...))

  return(out)
}

ode.1D    <- function (y,
                       times,
                       func,
                       parms,
                       nspec,
                       method= "lsode",
                       ...)
{
  if (hasArg(jacfunc)) stop ("cannot run ode.1D with jacfunc specified - remove jacfunc from call list")
  if (is.null(nspec)  ) stop ("cannot run ode.1D: nspec is not specified")
  N     <- length(y)

  if (is.character(func) || method=="lsodes")
  {
    if ( method != "lsodes") warning("ode.1D: R-function specified in a DLL-> integrating with lsodes") 
    # use lsodes
    out <- lsodes(y=y,times=times,func=func,parms,...)                    
 
  } else {
  # internal function #
  bmodel <- function (time,state,pars,model,...)
  {
     Modconc <-  model(time,state[ij],pars,...)   # ij: reorder state variables
     list(c(Modconc[[1]][ii]),Modconc[-1])        # ii: reorder rate of change     
  }

  if (is.character(func))stop ("cannot run ode.1D with R-function specified in a DLL")

  ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec
  ij    <- as.vector(t(matrix(nrow=nspec,1:N)))   # from ordering per spec -> per slice
    
  bmod  <- function(time,state,pars,...) {bmodel(time,state,pars,func,...)}
  if (is.null(method)) method <- "lsode" 
  if (method == "vode")
   out <- vode(y[ii],times,func=bmod,parms=parms,bandup=nspec,banddown=nspec,jactype="bandint",...) 
  else if (method == "lsode")
   out <- lsode(y[ii],times,func=bmod,parms=parms,bandup=nspec,banddown=nspec,jactype="bandint",...) 
  else if (method == "lsoda")
   out <- lsoda(y[ii],times,func=bmod,parms=parms,bandup=nspec,banddown=nspec,jactype="bandint",...) 
  else if (method == "lsodar")
   out <- lsodar(y[ii],times,func=bmod,parms=parms,bandup=nspec,banddown=nspec,jactype="bandint",...) 
   
  else
   stop ("cannot run ode.1D: method should be one of vode, lsoda, lsodar, lsode")   
  out[,(ii+1)] <- out[,2:(N+1)]
  }
  return(out)
}

ode.band  <- function (y,
                       times,
                       func,
                       parms,
                       nspec=NULL,
                       bandup=nspec,
                       banddown=nspec,
                       method = "lsode",
                       ...)
{
if (is.null(bandup)  ) stop ("cannot run ode.band: bandup is not specified")
if (is.null(banddown)) stop ("cannot run ode.band: banddown is not specified")
if (is.null(method)) method <- "lsode"
  if (method == "vode")
   vode(y,times,func,parms=parms,bandup=bandup,banddown=banddown,jactype="bandint",...) 
  else if (method == "lsode")
   lsode(y,times,func,parms=parms,bandup=bandup,banddown=banddown,jactype="bandint",...) 
  else if (method == "lsoda")
   lsoda(y,times,func,parms=parms,bandup=bandup,banddown=banddown,jactype="bandint",...) 
  else if (method == "lsodar")
   lsodar(y,times,func,parms=parms,bandup=bandup,banddown=banddown,jactype="bandint",...) 
  else
   stop ("cannot run ode.band: method should be one of vode, lsoda, lsodar or lsode")   
  
}
