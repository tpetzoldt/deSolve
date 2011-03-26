## find nearest event for each time step
nearestEvent <- function(times, events) {
  events  <- unique(events) # remove double events first
  ## sorting does not cost much if already sorted
  times <- sort(times)
  events <- sort(events)
  ## find index of events where time is between
  inearest <- findInterval(times, events)
  ## special care for smallest and biggest element
  lower <- events[pmax(inearest, 1)]
  upper <- events[pmin(inearest + 1, length(events))]
  nearest <- ifelse(times - lower < upper - times, lower, upper)
  return(nearest)
}

## remove times that are too close to an event
removeTooClose <- function(times, events, eps = .Machine$double.eps  * 10) {
  ## sorting does not cost much if already sorted
  ## sort times to ensure match of returned "nearest" value
  times <- sort(times)
  nearest <- nearestEvent(times, events)
  ## use bigger of the two numbers
  div <- pmax(times, nearest)
  ## special handling of zero
  div <- ifelse(div == 0, 1, div)
  reldiff <- abs(times - nearest) / div
  tooClose <- reldiff < eps
  times[!tooClose]
}