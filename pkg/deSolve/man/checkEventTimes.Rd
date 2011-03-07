\name{checkTimes}
\alias{checkTimes}
\title{
Check if Time Steps Are Contained in Events
}
\description{
This function checks a numeric vector (times) if it containes the values of 
another vector (events), avoiding values that are 'too close' too each other.
}
\usage{
checkTimes(times, events, eps = 1e-14, reldist = FALSE, silent = FALSE)
}
\arguments{
  \item{times}{
     a vector of time steps, that are checked for their distance to another vector
     (\code{events}),
}
  \item{events}{
     the vector with the events,
}
  \item{eps}{
     minimal difference that is allowed between any time step and an event,
}
  \item{reldist}{
    if \code{TRUE}, distance to the closest event is calculated as the ratio 
    between the absolute distance and the actual time step,
}
  \item{silent}{\code{FALSE} if true, text messages (warnings) are suppressed.
  }
}
\value{
  Vector with time steps containing all events, but without such time steps that 
  are very close (specified by \code{eps} to any event).
}
\details{
  In floating point arithmetics, problems can occur when values have to be 
  compared for 'equality', but also if two values are too close to each other.
  
  Function \code{checkEvents} returns a vector with time steps containing all 
  events, but without such time steps that are very close to an event. This means
  that all values of 'events' are contained but only the subset of 'times' that
  has no close neighbors in 'events'.
}
\author{
Thomas Petzoldt
}

\seealso{
  \code{\link{ode}}
}
\examples{

events <- sort(c(2,3,4 + 1e-12,5,7-1e-12, 7 + 6e-15, 7.5, 9, 24.9999, 25, 80))
times  <- sort(c(1:7, 4.5, 6.75, 7.5, 9.2, 9.0001, 25, 1e3))

options(digits=15)

checkTimes(times, events)
checkTimes(times, events, silent=TRUE)
checkTimes(times, events, reldist = TRUE)
checkTimes(times, events, eps=1e-2, reldist = TRUE)


}
\keyword{ misc }
