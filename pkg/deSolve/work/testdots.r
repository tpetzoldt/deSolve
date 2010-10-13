dots <- list(par = TRUE, lty =1:2, lwd = 2, las=c(NULL,1,2))

nm      <- names(dots)

setdots <- function (dots, np) {
  dd <- list()
  id <- 0 
  for (i in 1: length(dots)) {
   if (! is.null(dots[[i]])) {
    dd[[id<-id+1]] <- rep(dots[[i]], length = np)
    names(dd)[id] <- names(dots)[i]
   }
  }  
  return(dd)
} 
# works but not if there is somewhere a "NULL", e.g. as for las
ddots <- setdots(dots,3) 

extractdots <- function (dots, ii) {
  Dots <- list()
  for (i in 1:length(dots))
   if (! is.null(dots[[i]][ii])) {
    Dots[i] <- dots[[i]][ii]
    names(Dots)[i] <- names(dots)[i]   
   }
  Dots
}    
extractdots(ddots,1)
extractdots(ddots,2)