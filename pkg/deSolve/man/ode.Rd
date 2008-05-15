\name{ode}
\alias{ode}
\title{General solver for ordinary differential equations}
\description{Solves a system of ordinary differential equations.}
\usage{ode(y, times, func, parms, 
method=c("lsoda","lsode","lsodes","lsodar","vode","daspk"), ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{time sequence for which output is wanted; the first value of \code{times} must be the initial time}
  \item{func }{either a user-supplied function that computes the values of the derivatives in the ODE system (the model definition) at time t, or a character string giving the name of a compiled function in a dynamically loaded shared library. 
     		See details}
  \item{parms }{parameters passed to the integrator}
  \item{method }{the integrator to use, either a string ("lsoda","lsode","lsodes","lsodar","vode",or "daspk") or a function that performs integration, or
    a list of class \code{rkMethod}.}
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
See the selected integrator for the additional options}
\seealso{\item \code{\link{ode.band}}, and \code{\link{ode.1D}} when the jacobian is banded
 	\item \code{\link{aquaphy}}, \code{\link{ccl4model}}, where \code{ode} is used

	\item \code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}},
  \item  \code{\link{rk}}, \code{\link{rkMethod}}}
\keyword{math}

