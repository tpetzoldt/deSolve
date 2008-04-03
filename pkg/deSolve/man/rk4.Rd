\name{rk4}
\alias{rk4}
\title{Solve System of ODE (ordinary differential equation)s by
    classical Runge-Kutta 4th order integration.}
\description{Solving initial value problems for 
    systems of first-order ordinary differential equations
    (ODEs) using the classical Runge-Kutta 4th order integration.
    The system of ODE's may be written as an \R function (which
    may, of course, use \code{\link{.C}}, \code{.Fortran},
    \code{\link{.Call}}, etc., to call foreign code).
    A vector of parameters is passed to the ODEs, so the solver may
    be used as part of a modeling package for ODEs, or for parameter
    estimation using any appropriate modeling tool for non-linear models
    in \R such as \code{\link{optim}}, \code{\link[nls]{nls}},
    \code{\link{nlm}} or \code{\link[nlme]{nlme}}.}
\usage{rk4(y, times, func, parms, ...)}
\arguments{
  \item{y }{the initial values for the ode system.  If \code{y} has a
      name attribute, the names will be used to label the output matrix.}
  \item{times }{times at which explicit estimates for \code{y} are
      desired.  The first value in \code{times} must be the initial time.}
  \item{func }{a user-supplied function that computes the values of the
      derivatives in the ode system (the \emph{model defininition}) at time
      t.
      The user-supplied function  \code{func} must be called as:
      \code{yprime = func(t, y, parms)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ode system, and \code{parms} is a vector of parameters (which
      may have a names attribute, desirable in a large system).
  
      The return value of func should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose second element is a vector (possibly with a
      \code{\link{names}} attribute) of global values that are required at
      each point in \code{times}.}
  \item{parms }{vector or list holding the parameters used in \code{func} 
     that should be modifiable without rewriting the function.}
  \item{... }{additional arguments, allowing this to be a generic function}
}
\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}
\examples{## A simple resource limited Lotka-Volterra-Model
  lvmodel <- function(t, x, parms) {
    s <- x[1] # substrate
    p <- x[2] # producer
    k <- x[3] # consumer
    with(as.list(parms),{
      import <- approx(signal$times, signal$import, t)$y
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
  signal <- as.data.frame(list(times = times,
                              import = rep(0,length(times))))
  
  signal$import[signal$times >= 10 & signal$times <=11] <- 0.2
  
  ## Parameters for steady state conditions
  parms  <- c(a=0.0, b=0.0, c=0.1, d=0.1, e=0.1, f=0.1, g=0.0)
  
  ## Start values for steady state
  y<-xstart <- c(s=1, p=1, k=1)
  
  ## Classical RK4 with fixed time step
  out1  <- as.data.frame(rk4(xstart, times, lvmodel, parms))
  
  ## LSODA (default step size)
  out2 <- as.data.frame(lsoda(xstart, times, lvmodel, parms))
  
  ## LSODA: with fixed maximum time step
  out3 <- as.data.frame(lsoda(xstart, times, lvmodel, parms, hmax=1))
  
  par(mfrow=c(2,2))
  plot (out1$time, out1$s,  type="l",   ylim=c(0,3))
  lines(out2$time, out2$s, col="red",   lty="dotted")
  lines(out3$time, out3$s, col="green", lty="dotted")
  
  plot (out1$time, out1$p, type="l",    ylim=c(0,3))
  lines(out2$time, out2$p, col="red",   lty="dotted")
  lines(out3$time, out3$p, col="green", lty="dotted")
  
  plot (out1$time, out1$k, type="l",    ylim=c(0,3))
  lines(out2$time, out2$k, col="red",   lty="dotted")
  lines(out3$time, out3$k, col="green", lty="dotted")
  
  plot (out1$p, out1$k, type="l")
  lines(out2$p, out2$k, col="red",   lty="dotted")
  lines(out3$p, out3$k, col="green", lty="dotted")}
\references{Press, W. H.,  Teukolsky, S. A., Vetterling, W. T. and
             Flannery, B. P. (1992) Numerical Recipes in C. Cambridge
  	   University Press.}
\details{The method is implemented primarily for didactic purposes. Please use
    \link{lsoda}, \link{lsode} or another integrator for your real work!}
\seealso{\code{\link{lsoda}}}
\keyword{math}

