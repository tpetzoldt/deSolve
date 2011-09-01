library(deSolve)

ball <- function(t, y, parms) {
  dy1 <- y[2]
  dy2 <- -9.8     
  list(c(dy1, dy2))
}

yini  <- c(height = 0, velocity = 10)
times <- seq(from = 0, to = 20, by = 0.01)

rootfunc <- function(t, y, parms) return (y[1])

eventfunc <- function(t, y, parms) {
  y[1] <- 0
  y[2] <- -0.9 * y[2]
  return(y)
}

out <- ode(times = times, y = yini, func = ball,
           parms = NULL, rootfun = rootfunc,
           events = list(func = eventfunc, root = TRUE))

for (i in seq(1, 2001, 5)) {
  plot(out, which = "height", type = "l", lwd = 1,
       main = "", xlab = "Time", ylab = "Height")
  points(t(out[i, 1:2]), pch = 21, lwd = 1, col = 1, cex = 2,
         bg = rainbow(30, v = 0.6)[20 - abs(out[i, 3]) + 1])
  Sys.sleep(0.1-abs(out[i,3]/100)) # depends on velocity...
}
