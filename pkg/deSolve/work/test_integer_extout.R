library(deSolve)

m <- function(t, y, p) {
  dy <- c(0.1, 0.2)
  #a <- seq(0, 1, 0.2) # works
  a <- 1:5             # works for lsoda but was error for rk
  b <- rnorm(5)
  list(dy = dy, a = a, b = b)
}

out <- ode(c(x1 = 0, x2 = 0), seq(0, 1, 0.1), m, NULL, method = "rk4")

#out <- rk(c(x1 = 0, x2 = 0), seq(0, 1, 0.1), m, NULL, method = "rk4")

head(out)

subset(out, select="a")