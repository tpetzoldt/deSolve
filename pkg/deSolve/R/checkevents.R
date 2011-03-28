### ============================================================================
### Check events data set 
### ============================================================================

checkevents <- function (events, times, vars, dllname, root = FALSE) {

  if (is.null(events)) return(list())
  if (is.null(events$data) && is.null(events$func)) return(list())
  # only effective if lsodar, lsode,... "root" triggers an event, does not stop 
  if (root) {  # check if root should trigger an event...
    Root <- events$root
    if (is.null(Root)) Root <- 0
    Root <- as.integer(Root)
  } else Root <- as.integer(0)

  Rootsave <- events$maxroot
  if (is.null(Rootsave)) Rootsave <- 100  # number of roots to save.
  if (Rootsave < 0)
    stop("events$Rootsave should be > 0 in events")

  funevent <- events$func
  if (!is.null(funevent)) {
    if (is.character(funevent)){ 
     if (is.null(dllname))
       stop("'dllname' should be given if 'events$func' is a string")
     if (is.loaded(funevent, PACKAGE = dllname, type = "") ||
     is.loaded(funevent, PACKAGE = dllname, type = "Fortran")) {
       funevent <- getNativeSymbolInfo(funevent, PACKAGE = dllname)$address
     } else
       stop(paste("'events$func' should be loaded ",funevent))
       Type <- 3  
    } else {
      Type <- 2  # SHOULD ALSO CHECK THE FUNCTION if R-function....
      if (!is.null(dllname))
       stop("'events$func' should be a string, events specified in compiled code if 'dllname' is not NULL")
    }
    if (Root == 0) {
      if (is.null(events$time)) 
        stop("'events$time' should be given and contain the times of the events, if 'events$func' is specified and no root function")
      Time <- as.double(events$time) 
      # Karline: added this extra check ....
      #if (prod(Time %in% times) != 1)
      if (any(!(Time %in% times))) {      #thpe changed "prod" test to "any ..."
        #stop ("Not all event times 'events$times' are in output 'times'; include event$times in 'times'")
        warning("Not all event times 'events$times' where in output 'times' so they are automatically included.")
        uniqueTimes <- cleanEventTimes(times, Time)
        if (length(uniqueTimes) < length(times))
          warning("Some time steps were very close to events - only the event times are used in these cases.")
        times <- sort(c(uniqueTimes, Time))
      }
    } else Time <- min(times) - 1  # never reached....
      return (list (Time = Time, SVar = NULL, Value = NULL, 
        Method = NULL, Type = as.integer(Type), func = funevent,
        Rootsave = as.integer(Rootsave), Root = Root))

  }  ## Check the event data series
  event <- events$data
  if (is.matrix(event)) event <- as.data.frame(event) 

  if (ncol(event) < 3) 
    stop("'event' should have at least 3 columns: state variable, time, value")

  if (!is.data.frame(event)) 
    stop("'event' should be a data.frame with 3(4) columns: state variable, time, value, (method)")
    
  ## thpe: added the following check; makes check < 3 columns obsolete
  evtcols <-  c("var", "time", "value", "method")
  if (!all(evtcols %in% names(eventdat)))
    stop("structure of events does not match specification, see help('events')")
  
  ## thpe: make sure that event data frame has correct order
  eventdat <- eventdat[evtcols]

## variables, 1st column should be present
  if (is.factor(event[,1])) 
    event[,1] <- as.character(event[,1])

  if (is.character(event[,1]))  {
    vv <- match(event[,1], vars)
    if (any(is.na(vv)))
      stop("unknown state variable in 'event': ", paste(event[,1][which(is.na(vv))], ","))
    event[,1] <- vv  
  } else if (max(event[,1]) > length(vars))
      stop("too many state variables in 'event'; should be < ", paste(length(vars)))

## 2nd and 3rd columns should be numeric
  if (!is.numeric(event[,2])) 
      stop("times in 'event', 2nd column should be numeric")

  if (!is.numeric(event[,3]))  
      stop("values in 'event', 3rd column should be numeric")

## Times in 'event' should be embraced by 'times'
  rt <- range(times)
  ii <- c(which(event[,2] < rt[1]), which(event[,2] > rt[2]))
  if (length(ii) > 0) 
    event <- event [-ii,]

  if (any(!(event[,2] %in% times))) {
    warning("Not all event times 'events$times' where in output 'times' so they are automatically included.")
    uniqueTimes <- cleanEventTimes(times, event[,2])
    if (length(uniqueTimes) < length(times))
      warning("Some time steps were very close to events - only the event times are used in these cases.")
    times <- sort(c(uniqueTimes, event[,2]))
  }  


## 4th column: method; if not available: "replace" = method 1 - to date: 3 methods
  if (ncol(event) ==3) 
    event$method <- rep(1,nrow(event))
  else if (is.numeric(event[,4])) { 
      if (max(event[,4]) > 3 | min(event[,4]) < 1)
        stop("unknown method in 'event': should be >0 and < 4") 
  } else {
    vv <- charmatch(event[,4],c("replace", "add", "multiply"))
    if (any(is.na(vv)))
      stop("unknown method in 'event': ", paste(event[,3][which(is.na(vv))],","),
        " should be one of 'replace', 'add', 'multiply'")
    event$method <- vv
  }

## Check the other events elements (see optim code)
  con <- list(ties = "notordered", time = NULL, data = NULL, func = NULL, root = NULL)
  nmsC <- names(con)
  con[(namc <- names(events))] <- events
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
     warning("unknown names in events: ", paste(noNms, collapse = ", "))

## Check what needs to be done in case the time series is not "ordered"

  if (!identical(con$ties, "ordered")) { # see approx code

## first order with respect to time (2nd col), then to variable (1st col)
    if(length(x <- unique(event[,1:2])) < nrow(event)){
      ties <- mean
      if (missing(ties))
        warning("collapsing to unique 'x' values")
      event <- aggregate(event[,c(3, 4)], event[,c(1, 2)], ties)
    }
  }

  return (list (Time = as.double(event[,2]), SVar = as.integer(event[,1]), 
    Value = as.double(event[,3]), Method = as.integer(event[,4]), 
    Rootsave = as.integer(Rootsave), 
    Type = as.integer(1), Root = Root, newTimes = times))
}

