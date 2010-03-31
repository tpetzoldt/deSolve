require(deSolve)
set.seed(123)

func1 <- function(t,y,p) {
  cat("----------------\n")
  cat(t, "outer: ", timestep(TRUE), timestep(FALSE),"\n")
  ff <- function(t,y,p) {
    cat(t, "inner ", timestep(TRUE), timestep(FALSE), "\n")
    list(-0.1*y)
   }
   #AA <- ode(y=5,times=c(0.15, 0.4), p=NULL, func=ff, method="euler")
   AA <- ode(y=5,times=c(0.15, 0.4, 0.77), p=NULL,
           func=ff, method="lsoda") # stack imbalance !!!
   #AA <- euler(y=5,times=sort(sample(1:10, 3)/10), p=NULL, func=ff) # special "euler"

  list(-0.1*y)
}
out <- rk(y=8,times=sort(sample(1:20, 5)), p=NULL,func=func1, method="euler")
## can be called repeatedly without problems

out <- ode(y=8, times=sort(sample(1:20, 5)), p=NULL,func=func1)
## ---> crashes R after repeated calls

