## =============================================================================
## A bouncing ball; ode with event location
## =============================================================================
require(deSolve)
#-----------------------------
# the model function
#-----------------------------
ballode<- function(t, y, parms) {
  dy1 <- y[2]
  dy2 <- -9.8
  list(c(dy1, dy2))
}

#-----------------------------
# the root and event function
#-----------------------------
# event triggered when the ball hits the ground (height =0)
root <- function(t, y, parms) y[1]

# bouncing
event <- function(t, y, parms) {
  y[1] <- 0
  y[2] <- -0.9 * y[2]
 return(y)
}

#-----------------------------
# initial values and times
#-----------------------------
yini  <- c(height = 0, v = 20)
times <- seq(0, 20, 0.01)

#-----------------------------
# solve the model
#-----------------------------
out   <- lsodar(times = times, y = yini, func = ballode, parms = NULL,
  events = list(func = event, root = TRUE), rootfun = root)

out2   <- lsode(times = times, y = yini, func = ballode, parms = NULL,
  events = list(func = event, root = TRUE), rootfun = root, verbose=TRUE)

attributes(out)$troot
attributes(out2)$troot
#-----------------------------
# display, plot results
#-----------------------------
for (i in seq(1, 2001, 10)) {
  #png(filename=paste("image", i+1000, ".png", sep=""), width=300, height=100)
  #par(mar=rep(.1,4))
  plot(out,which = "height", type="l", lwd=1,
       #axes=FALSE, main = "", xlab="", ylab = ""
       main = "Bouncing Ball", xlab="Time", ylab = "Height"
  )
  box()
  points(t(out[i,1:2]), pch=21, , lwd=1, col=1, cex=2,
  bg=rainbow(30, v=0.6)[20-abs(out[i,3])+1])
  Sys.sleep(.01)
  #dev.off()
}
