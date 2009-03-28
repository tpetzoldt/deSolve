

### rk4, euler
### Shortcut functions for Euler and Runge-Kutta 4th order integration



rk4 <- function(y, times, func, parms, hini=min(diff(times)), verbose = FALSE, ...) {
  rk(y, times, func, parms, hmin = hini, hini = hini,
    method = rkMethod("rk4", ...), verbose = verbose, maxsteps = 1e6, ...)
}

euler <- function(y, times, func, parms, hini=min(diff(times)), verbose = FALSE, ...) {
  rk(y, times, func, parms, hmin = hini, hini = hini,
    method = rkMethod("euler", ...), verbose = verbose, maxsteps = 1e6, ...)
}
