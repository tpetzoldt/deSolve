\name{checkEventTimes}
\alias{checkEventTimes}
\title{
Check if All Event Times Are Contained in Vector Times
}
\description{
This function checks if all values of a vector (events) are contained in
another vector (times), avoiding values that are 'too close' too each other.
}
\usage{
checkEventTimes(events, times, eps = 1e-12, reldist = TRUE, silent = FALSE)
}
\arguments{
  \item{events}{
    the vector with the events,
  }
  \item{times}{
    a vector of time steps, which is checked whether the times of all
    (\code{events}) are contained,
  }
  \item{eps}{
    tolerance that is allowed for comparing \code{events} and \code{times},
  }
  \item{reldist}{
    if \code{TRUE}, distance to the closest event is calculated as the ratio 
    between the absolute distance and the actual time step, otherwise
    absolute tolerance is used.
  }
  \item{silent}{if \code{TRUE}, text messages (warnings) are suppressed.
  }
}
\value{
  Returns a list with two elements:
  \itemize{
    \item{value}{that is \code{TRUE} if all events where in times and
      \code{FALSE} otherwise,}
    \item{times}{the vector with the new time steps containing all
      events, but without such time steps that 
      are very close (specified by \code{eps} to any event).}
  }
}
\details{
  In floating point arithmetics, problems can occur when values have to be 
  compared for 'equality', but also if two values are too close to each other.
  
  This utility function can be used to check if all \code{events} are exactly
  contained in times and then tries to solve it by returning a vector 
  that contains all events, but without such time steps that are very close
  to an event. This means that all values of \code{events} are contained
  but only the subset of \code{times} that have no close neighbors in
  \code{events}.
}
\author{
Thomas Petzoldt
}

\seealso{
  \code{\link{ode}}
}
\examples{

events <- sort(c(2, 3, 4 + 1e-10, 5, 7 - 1e-10,
                 7 + 6e-15, 7.5, 9, 24.9999, 25, 80))
times  <- sort(c(1:7, 4.5, 6.75, 7.5, 9.2, 9.0001, 25, 1e3))

options(digits = 15)

checkEventTimes(events, times, reldist = FALSE)
checkEventTimes(events, times, silent = TRUE)
checkEventTimes(events, times) # default: reldist = TRUE

## use rather big value for tolerance and return only the new time steps
checkEventTimes(events, times, eps = 1e-2, silent = TRUE)$times

## ====================================================================
## Example how checkEvents is used in practise 
## ====================================================================

eventmod <- function(t, var, parms) {
  list(dvar = rep(0, 2))
}

yini <- c(v1 = 1, v2 = 2)

## assume a vector of times with a small numerical error
times <- seq(0, 10, by = 0.1) + 1e - 12

## events that contains the 'exact' event times
eventdat <- data.frame(var = c("v1", "v2", "v2", "v1"),
                       time = c(1, 1, 5, 9) ,
                       value = c(1, 2, 3, 4),
                       method = c("add", "mult", "rep", "add"))
eventdat

checked <- checkEventTimes(eventdat$time, times, eps = 1e-10, silent = TRUE)
cat("The check returned:", checked$value, "\n")

## simulation with events not exactly contained in times
## would stop with an error
\dontrun{
out <- ode(yini, times, eventmod, parms = NULL, 
  events = list(data = eventdat))
}

## simulation with fixed time steps
out <- ode(yini, checked$times, eventmod, parms = NULL,
  events = list(data = eventdat))

plot(out, type = "l")

}
\keyword{ misc }
