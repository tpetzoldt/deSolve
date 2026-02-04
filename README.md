## R Package deSolve

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/deSolve)](https://CRAN.R-project.org/package=deSolve)
[![Dev version](https://img.shields.io/github/r-package/v/tpetzoldt/deSolve?label=dev%20version)](https://github.com/tpetzoldt/deSolve)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/deSolve)](https://cran.r-project.org/package=deSolve)
[![CRAN downloads per month](https://cranlogs.r-pkg.org/badges/deSolve)](https://cran.r-project.org/package=deSolve)
<!-- badges: end -->

Solvers for Initial Value Problems of Differential Equations (ODE, DAE, DDE)

### Overview

**deSolve** is a comprehensive R package for solving initial value problems of differential equations. It provides robust numerical solvers for:

- Ordinary Differential Equations (ODE)
- Differential Algebraic Equations (DAE)
- Partial Differential Equations (PDE) using the method of lines
- Delay Differential Equations (DDE)

The package can be used in scientific computing, ecological modeling, pharmacokinetics, and many other fields requiring numerical integration of dynamic systems.

### Features

* Multiple solver algorithms (lsoda, lsode, lsodes, vode, ode45, rk4, euler, ...)
* Support for compiled code (C, Fortran) for improved performance
* Forcing functions, event handling, and plotting

### Installation

**Stable version from CRAN (recommended)**

```
install.packages("deSolve")
```

**Development version from Github (unstable)**

```
# install.packages("remotes")
remotes::install_github("tpetzoldt/deSolve")
```


### Quick Example

```
library(deSolve)

## Chaos in the Atmosphere
Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <-  a * X + Y * Z
    dY <-  b * (Y - Z)
    dZ <- -X * Y + c * Y - Z
    list(c(dX, dY, dZ))
  })
}

parameters <- c(a = -8/3, b = -10, c = 28)
state      <- c(X = 1, Y = 1, Z = 1)
times      <- seq(0, 100, by = 0.01)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)

plot(out)
```

## License 

**deSolve** is Free and Open Source Software, released under the
[GPL 2.0](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html) or 
[GPL 3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) license.


### Documentation

For comprehensive documentation and tutorials, visit:

- Package website: https://cran.r-project.org/package=deSolve
- Example models collection: https://github.com/tpetzoldt/dynamic-R-models
- and more

### Citation

If you use deSolve in your research, please cite:

Soetaert, K., Petzoldt, T., and Setzer, R. W. (2010). Solving Differential Equations in R: Package deSolve. 
Journal of Statistical Software, 33(9), 1-25. https://doi.org/10.18637/jss.v033.i09 

## Contributing

This is the active development repository for deSolve. Previous development occured on R-Forge at http://desolve.r-forge.r-project.org/

----

Thomas & Karline  2026-02-04

