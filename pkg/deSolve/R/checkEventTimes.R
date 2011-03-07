checkTimes <- function (times, events, eps = 1e-14, reldist = FALSE, silent=FALSE) {
  events  <- unique(events) # remove double events first
  nevents <- length(events)
  ntimes  <- length(times)
  
  if (events[1] <= times[1]) stop("first time step must occur before first event")
  if (events[nevents] >= times[ntimes]) stop("last time step must occur after last event")

  if (any(!(events %in% times))) {
    ## thpe: improve warning message, mention numerical issues and what we do
    if (! silent) {
      warning ("Not all event times are exactly in output 'times'; they will be included")
    }
    x <- numeric(length(times)) # reserve memory for output
    xout <- .C("checktimes", 
            as.integer(ntimes),
            as.integer(nevents),
            as.double(times),
            as.double(events),
            as.integer(reldist), # relative = TRUE, absolute = FALSE
            x=as.double(x))$x
    # df <- data.frame(times=times, xout=xout) # see how it works
    # View(df)
    times <- sort(c(times[xout > eps], events))
  }
  # else return times unchanged
  return(times)
}


      




