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
times <- seq(0, 35, length.out=2001)

#-----------------------------
# solve the model
#-----------------------------
out   <- lsodar(times = times, y = yini, func = ballode, parms = NULL,
  events = list(func = event, root = TRUE), rootfun = root)

attributes(out)$troot

#------------------------------------------
# make time steps proportional to 1/speed
#------------------------------------------
rel   <- cumsum(pmin(abs(1/out[,"v"]), 1)) # 1: limit time spent at top position
f     <- approxfun(rel, times, rule=2)
t.rel <- f(seq(0, max(rel), length.out=101))

# add roots to make animation nicer (may be made even better ...)
t.rel <- sort(c(t.rel, attributes(out)$troot))


#------------------------------------------
# run again and show results
#------------------------------------------
out2  <- lsodar(times = t.rel, y = yini, func = ballode, parms = NULL,
                events = list(func = event, root = TRUE), rootfun = root)
for (i in 1:length(t.rel)) {
  #png(filename=paste("bball/image", i+1000, ".png", sep=""),
  #  type="cairo", width=300, height=100, antialias="subpixel")
  par(mar=rep(.1, 4))
  plot(out,which = "height", type="l", lwd=2, col="grey",
       axes=FALSE, main = "", xlab="", ylab = ""
       #main = "Bouncing Ball", xlab="Time", ylab = "Height"
  )
  box()
  points(t(out2[i,1:2]), pch=21, lwd=1, col=1, cex=2,
         bg=rainbow(30, v=0.6)[20-abs(out2[i,3])+1])
  Sys.sleep(.1)
  #dev.off()
}
