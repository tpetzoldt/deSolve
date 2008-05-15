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
  This function is normally called by \code{\link{rk}}
}
\value{
  A list with the following elements:
  \item{ID}{name of the method (character)}
  \item{varstep}{boolean value specifying if the method allows for variable time step (\code{TRUE}) or not (\code{FALSE}).}
    \item{FSAL}{(first step as last) boolean value specifying if the method allows re-use of the last function evaluation (\code{TRUE}) or not (\code{FALSE}).}
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
  Runge, C. (1895) Ueber die numerische Aufloesung von Differentialgleichungen,
  Math. Ann. 46, 167-178.

  Kutta, W. (1901) Beitrag zur naeherungsweisen Integration totaler Differentialgleichungen,
  Z. Math. Phys. 46, 435-453.

  Dormand, J. R. and Prince, P. J. (1980) A family of embedded Runge-Kutta formulae,
  J. Comput. Appl. Math. 6(1), 19-26.

  Dormand, J. R. and Prince, P. J. (1981) High order embedded Runge-Kutta formulae,
  J. Comut. Appl. Math. 7(1), 67-75.

  Fehlberg, E. (1967) Klassische Runge-Kutta-Formeln fuenfter and siebenter Ordnung mit
  Schrittweiten-Kontrolle, Computing (Arch. Elektron. Rechnen) 4, 93-106.

  Butcher, J. C. (1987) The numerical analysis of ordinary differential equations,
  Runge-Kutta and general linear methods, Wiley, Chichester and New York.

  Engeln-Muellges, G. and Reutter, F. (1996) Numerik Algorithmen:
  Entscheidungshilfe zur Auswahl und Nutzung. VDI Verlag, Duesseldorf.
}
\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}

\seealso{\code{\link{rk}}}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }

