checkEventTimes <- function (events, times, eps = 1e-12, reldist = TRUE, silent = FALSE) {
  events  <- unique(events) # remove double events first
  nevents <- length(events)
  ntimes  <- length(times)
  value   <- TRUE # assume, all events are in times

  if (events[1] <= times[1])
      stop("first time step must be less or equal then first event")
  if (events[nevents] >= times[ntimes])
      stop("last time step must be equal or greater then last event")

  if (any(!(events %in% times))) {
    value <- FALSE
    if (! silent) {
      warning ("Times did not contain all events, so they were included.")
    }
    x <- numeric(length(times)) # reserve memory for output
    ## this .C function finds nearest event for each value in times
    ## and returns its relative resp. absolute distance.
    xout <- .C("checkeventtimes",
            as.integer(nevents),
            as.integer(ntimes),
            as.double(events),
            as.double(times),
            as.integer(reldist), # relative = TRUE, absolute = FALSE
            x=as.double(x))$x
    ## thpe: code for debugging, remove later
    df <- data.frame(times=times, xout=xout) # see how it works
    View(df)
    times <- sort(c(times[xout > eps], events))
  } else {
    ## else return times unchanged
    cat("o.k., all events are contained in times.\n")
  }
  return(list(value = value, times = times))
}







