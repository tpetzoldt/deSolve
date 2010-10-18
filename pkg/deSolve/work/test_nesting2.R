require(deSolve)
set.seed(123)

func1 <- function(t,y,p) {
  cat("----------------\n")
  cat(t, "outer: ", timestep(TRUE), timestep(FALSE),"\n")
  ff <- function(t,y,p) {
    cat(t, "inner ", timestep(TRUE), timestep(FALSE), "\n")
    list(-0.1 * y)
   }
   oo <- ode(y = y,times=1:3, p=NULL, func=ff, method="euler")
   oo <- as.data.frame(oo)
   y2 <- unlist(oo[3, 2:4])
   list(0.1 * y2)
}
out <- rk(y=c(-1,0,1), times=1:20, p=NULL,func=func1, method="rk4")
plot(out)


out <- ode(y=c(-1,0,1), times=1:20, p=NULL,func=func1, hmax=1, method="daspk")
plot(out)