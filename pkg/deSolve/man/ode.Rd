\name{ode}
\alias{ode}
\title{General solver for ordinary differential equations}
\description{Solves a system of ordinary differential equations.}
\usage{ode(y, times, func, parms, 
method=c("lsoda","lsode","lsodes","lsodar","vode","daspk", "euler", "rk4", 
         "ode23", "ode45"), ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{time sequence for which output is wanted; the first value of \code{times} must be the initial time}
  \item{func }{either an R-function that computes the values of the derivatives in the ODE system (the model definition) at time t, or a character string giving the name of a compiled function in a dynamically loaded shared library. 
      If \code{func} is an R-function, it must be called as:
      \code{yprime = func(t, y, parms, ...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ODE system.  \code{parms} is a vector of parameters.
  
      The return value of \code{func} should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose next elements are global values that are required at
      each point in \code{times}. }
  \item{parms }{parameters passed to \code{func}}
  \item{method }{the integrator to use, either a string ("lsoda","lsode","lsodes",
    "lsodar","vode", "daspk", "euler", "rk4", "ode23" or "ode45") or a function
     that performs integration, or a list of class \code{rkMethod}.}
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
\details{This is simply a wrapper around the various ode solvers.
\cr See \code{\link{lsoda}} for information about specifying the model in compiled code.
\cr See the selected integrator for the additional options}
\seealso{\item \code{\link{ode.band}}, and \code{\link{ode.1D}} when the jacobian is banded
 	\item \code{\link{aquaphy}}, \code{\link{ccl4model}}, where \code{ode} is used

	\item \code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}},
  \item  \code{\link{rk}}, \code{\link{rkMethod}}}
\keyword{math}
\examples{

  ## A simple resource limited Lotka-Volterra-Model
  ## The parameters 
  parms  <- c(a=0.0, b=0.0, c=0.1, d=0.1, e=0.1, f=0.1, g=0.0)


  ## The model
  lvmodel <- function(t, x, parms, input)  {
      with(as.list(c(parms,x)),  {
      import <- input(t)
      dS <- import - b*S*P + g*K
      dP <- c*S*P  - d*K*P
      dK <- e*P*K  - f*K
      res<-c(dS, dP, dK)
      list(res)                   })
    }
  
  ## vector of timesteps
  times  <- seq(0, 100, length=101)
  
  ## external signal with rectangle impulse
  signal <- as.data.frame(list(times = times,
                              import = rep(0,length(times))))
  
  signal$import[signal$times >= 10 & signal$times <=11] <- 0.2
  
  sigimp <- approxfun(signal$times, signal$import, rule=2)
  
  
  ## Start values for steady state
   xstart <- c(S=1, P=1, K=1)
  
  ## Solve model
  out <- as.data.frame(ode(y=xstart,times= times, 
                       func=lvmodel, parms, input =sigimp))
  plot(out$P,out$K,type="l",lwd=2)
}

