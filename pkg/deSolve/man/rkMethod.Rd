\name{rkMethod}
\alias{rkMethod}

\title{Collection of parameter sets (Butcher arrays) for common methods
  of the Runge-Kutta family of ODE solvers.}
\description{
  This function returns a list specifying coefficients and properties of
  ODE solver methods from the Runge-Kutta family.
}
\usage{
rkMethod(method = NULL, ...)
}

\arguments{
  \item{method}{a string constant naming one of the
    pre-defined methods of the Runge-Kutta family of solvers.
    The most common methods are the fixed-step methods "euler", "rk2", "rk4" 
    or the variable step methods "rk34", "rk45f" or "rk45m7".}
  \item{\dots}{specification of user-defined solver, see example below.}
}
\details{
  This function is normally called by \code{\link{rk4}}
}
\value{
  A list with the following elements:
  \item{ID}{name of the method (character)}
  \item{varstep}{boolean value specifying if the method allows for variable time step (\code{TRUE}) or not (\code{FALSE}).}
  \item{A}{coefficient matrix of the method. 
    As \code{link{rk}} supports only explicit methods, this matrix must be lower triangular.
    \code{A} can also be a vector if only the subdiagonal values are different from zero}
  \item{b1}{weighting coefficients for averaging the function evaluations of method 1}
  \item{b2}{weighting coefficients for averaging the function evaluations of method 2 (for embedded methods that allow variable time step)}
  \item{c}{coefficients for calculating the intermediate time steps}  
  \item{stage}{number of function evaluations needed (corresponds to number of rows in A)}  
  \item{Qerr}{global error order of the method, important for automatic time-step adjustment.}    
  ...
}
\references{
  Press, W. H.,  Teukolsky, S. A., Vetterling, W. T. and
  Flannery, B. P. (1992) Numerical Recipes in C. Cambridge
  University Press.
  
  %% add more references here; textbooks and original references
  }
\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}

\seealso{\code{\link{rk}}}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }

