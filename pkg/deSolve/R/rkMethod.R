### ============================================================================
### Butcher tables for selected explicit ODE solvers of Runge-Kutta type
### Note that for fixed step methods A is a vector (the subdiagonal of matrix A)
###   For variable time step methods, A must be strictly lower triangular.
###   The underlying rk code is currently restricted to explicit methods.
### ============================================================================

rkMethod <- function(method = NULL, ...) {
  methods <- list(
    euler = list(ID = "euler",
        varstep = FALSE,
          A      = c(0),
          b1     = c(1),
          c      = c(0),
          stage  = 1,
          Qerr   = 1
    ),
    ## Heun's method
    rk2 = list(ID = "rk2",
        varstep = FALSE,
          A      = c(0, 1),
          b1     = c(0.5, 0.5),
          c      = c(0, 1),
          stage  = 2,
          Qerr   = 1
    ),
    ## classical Runge-Kutta 4th order method
    rk4 = list(ID = "rk4",
        varstep = FALSE,
          A      = c(0, .5, .5, 1),
          b1     = c(1/6, 1/3, 1/3, 1/6),
          c      = c(0, .5, .5, 1),
          stage  = 4,
          Qerr   = 4
    ),
    ## One of the numerous RK23 formulae
    rk23 = list(ID = "rk23",
      varstep = TRUE,
      FSAL    = FALSE,
      A  = matrix(c(0, 0, 0,
                  1/2, 0, 0,
                  -1, 2, 0), 3, 3, byrow = TRUE),
      b1 = c(0, 1, 0),
      b2 = c(1/6, 2/3, 1/6),
      c  = c(0, 1/2, 2),
      stage = 3,
      Qerr  = 2
    ),
    ## Bogacki & Shampine
    rk23bs = list(ID = "rk23bs",
      varstep = TRUE,
      FSAL    = TRUE,
      A  = matrix(c(0, 0, 0, 0,
                  1/2, 0, 0, 0,
                  0, 3/4, 0, 0,
                  2/9, 1/3, 4/9, 0), 4, 4, byrow = TRUE),
      b1 = c(7/24, 1/4, 1/3, 1/8),
      b2 = c(2/9, 1/3, 4/9, 0),
      c  = c(0, 1/2, 3/4, 1),
      stage = 4,
      Qerr  = 2
    ),
    ## RK-Fehlberg 34
    rk34f = list(ID = "rk34f",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0,
                      2/7, 0, 0, 0,
                      77/900, 343/900, 0, 0,
                      805/1444, -77175/54872, 97125/54872, 0,
                      79/490, 0, 2175/3626, 2166/9065),
                      5, 4, byrow = TRUE),
         b1 = c(79/490, 	0, 2175/3626, 2166/9065, 	0),
         b2 = c(229/1470, 0, 1125/1813, 13718/81585, 1/18),
         c  = c(0,	2/7, 	7/15, 35/38, 	1),
         stage = 5,
         Qerr  = 3
    ),
    ## RK-Fehlberg Method 45
    rk45f = list(ID = "rk45f",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/4, 0, 0, 0, 0,
                      3/32, 9/32, 0, 0, 0,
                      1932/2197, -7200/2197, 7296/2197, 0, 0,
                      439/216, -8, 3680/513, -845/4104, 0,
                      -8/27, 2, -3544/2565, 1859/4104, -11/40),
                      6, 5, byrow = TRUE),
         b1 = c(25/216, 	0, 	1408/2565, 	2197/4104, 	-1/5, 	0),
         b2 = c(16/135, 	0, 	6656/12825, 	28561/56430, 	-9/50, 	2/55),
         c  = c(0,	1/4, 	3/8, 	12/13, 	1, 	1/2),
         stage = 6,
         Qerr  = 4
    ),
    ## Cash-Karp method - td = 2: dense output type '2'
    rk45ck = list(ID = "rk45ck",
         varstep = TRUE,
         FSAL = TRUE,
         A = matrix(c(0,    0,       0,         0,            0,
                      1/5,  0,       0,         0,            0,
                      3/40, 9/40,    0,         0,            0,
                      3/10, -9/10,   6/5,       0,            0,
                    -11/54, 5/2,    -70/27,     35/27,        0,
                   1631/55296, 175/512, 575/13824, 44275/110592, 253/4096),
                      6, 5, byrow = TRUE),
         b1 = c(2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4),
         b2 = c(37/378, 0, 250/621, 125/594, 0, 512/1771),
         c = c(0, 1/5, 3/10, 3/5,  1, 7/8),
         densetype = 2,        # special dense output type 2
         stage = 6,
         Qerr = 4),
    ## England Method
    rk45e = list(ID = "rk45e",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/2, 0, 0, 0, 0,
                      1/4, 1/4, 0, 0, 0,
                      0, -1, 2, 0, 0,
                      7/27, 10/27, 0, 1/27, 0,
                      28/625, -125/625, 546/625, 54/625, -378/625),
                      6, 5, byrow = TRUE),
         b1 = c(1/6, 	0, 4/6, 1/6, 	0, 	0),
         b2 = c(14/336, 0, 0,	35/336, 162/336, 125/336),
         c  = c(0,	1/2, 	1/2, 	1, 	2/3, 	1/5),
         stage = 6,
         Qerr  = 4
    ),
    ## Prince-Dormand 5(4)6m
    rk45dp6 = list(ID = "rk45dp6",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/5, 0, 0, 0, 0,
                      3/40, 9/40, 0, 0, 0,
                      3/10, -9/10, 6/5, 0, 0,
                      226/729, -25/27, 880/729, 55/729, 0,
                      -181/270, 5/2, -266/297, -91/27, 189/55),
                      6, 5, byrow = TRUE),
         b1 = c(31/540, 	0, 190/297, -145/108, 351/220, 1/20),
         b2 = c(19/216, 0, 1000/2079,	-125/216, 81/88, 5/56),
         c  = c(0,	1/5, 	3/10, 3/5, 	2/3, 	1),
         stage = 6,
         Qerr  = 4
    ),
    ## Prince-Dormand 5(4)7m -- recommended by the Octave developers
    rk45dp7 = list(ID = "rk45dp7",
         varstep = TRUE,
         FSAL    = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0, 0,
                      1/5, 0, 0, 0, 0, 0,
                      3/40, 9/40, 0, 0, 0, 0,
                      44/45, -56/15, 32/9, 0, 0, 0,
                      19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,
                      9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,
                      35/384, 0, 500/1113, 125/192, -2187/6784, 11/84),
                      7, 6, byrow = TRUE),
         b1 = c(5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40),
         b2 = c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0),
         c  = c(0, 1/5, 3/10, 4/5, 8/9, 1, 1),
         d  = c(-12715105075.0/11282082432.0, 0, 87487479700.0/32700410799.0,
                -10690763975.0/1880347072.0, 701980252875.0/199316789632.0,
                -1453857185.0/822651844.0, 69997945.0/29380423.0),
         densetype = 1, # default type of dense output formula, if available
         stage = 7,
         Qerr  = 4
    ),
    ## Prince-Dormand 78 method
    rk78dp = list(ID = "rk78dp",
         varstep = TRUE,
         FSAL = FALSE,
         ## ToDo: use common fractions instead of decimals
         A = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,
           5.55555555555555555555555555556e-2,0,0,0,0,0,0,0,0,0,0,0,
           2.08333333333333333333333333333e-2,6.25e-2,0,0,0,0,0,0,0,0,0,0,
           3.125e-2,0,9.375e-2,0,0,0,0,0,0,0,0,0,3.125e-1,0,-1.171875,1.171875,
           0,0,0,0,0,0,0,0,3.75e-2,0,0,1.875e-1,1.5e-1,0,0,0,0,0,0,0,
           4.79101371111111111111111111111e-2,0,0,1.12248712777777777777777777778e-1,
           -2.55056737777777777777777777778e-2,1.28468238888888888888888888889e-2,0,0,0,0,0,0,
           1.6917989787292281181431107136e-2,0,0,3.87848278486043169526545744159e-1,
           3.59773698515003278967008896348e-2,1.96970214215666060156715256072e-1,
           -1.72713852340501838761392997002e-1,0,0,0,0,0,
           6.90957533591923006485645489846e-2,0,0,-6.34247976728854151882807874972e-1,
           -1.61197575224604080366876923982e-1,1.38650309458825255419866950133e-1,
           9.4092861403575626972423968413e-1,2.11636326481943981855372117132e-1,0,0,0,0,
           1.83556996839045385489806023537e-1,0,0,-2.46876808431559245274431575997e0,
           -2.91286887816300456388002572804e-1,-2.6473020233117375688439799466e-2,
           2.84783876419280044916451825422e0,2.81387331469849792539403641827e-1,
           1.23744899863314657627030212664e-1,0,0,0,
           -1.21542481739588805916051052503,0,0,
           1.66726086659457724322804132886e1,9.15741828416817960595718650451e-1,
           -6.05660580435747094755450554309,-1.60035735941561781118417064101e1,
           1.4849303086297662557545391898e1,-1.33715757352898493182930413962e1,
           5.13418264817963793317325361166,0,0,
           2.58860916438264283815730932232e-1,0,0,-4.77448578548920511231011750971,
           -4.3509301377703250944070041181e-1,-3.04948333207224150956051286631,
           5.57792003993609911742367663447,6.15583158986104009733868912669,
           -5.06210458673693837007740643391,2.19392617318067906127491429047,
           1.34627998659334941535726237887e-1,0,
           8.22427599626507477963168204773e-1,0,0,-1.16586732572776642839765530355e1,
           -7.57622116690936195881116154088e-1,7.13973588159581527978269282765e-1,
           1.20757749868900567395661704486e1,-2.12765911392040265639082085897,
           1.99016620704895541832807169835,-2.34286471544040292660294691857e-1,
           1.7589857770794226507310510589e-1,0), nrow = 13, ncol = 12,
           byrow = TRUE),

         b1 = c(2.9553213676353496981964883112e-2,0,0,0,0,
           -8.28606276487797039766805612689e-1,3.11240900051118327929913751627e-1,
           2.46734519059988698196468570407,-2.54694165184190873912738007542,
           1.44354858367677524030187495069,7.94155958811272872713019541622e-2,
           4.44444444444444444444444444445e-2,0.),

         b2 = c(4.17474911415302462220859284685e-2,0,0,0,0,
           -5.54523286112393089615218946547e-2,2.39312807201180097046747354249e-1,
           7.0351066940344302305804641089e-1,-7.59759613814460929884487677085e-1,
           6.60563030922286341461378594838e-1,1.58187482510123335529614838601e-1,
           -2.38109538752862804471863555306e-1, 2.5e-1),

         c = c(0,5.55555555555555555555555555556e-2,
           8.33333333333333333333333333334e-2,1.25e-1,3.125e-1,3.75e-1,
           1.475e-1,4.65e-1,5.64865451382259575398358501426e-1,6.5e-1,
           9.24656277640504446745013574318e-1,1,1), stage = 13, Qerr = 7
    ),

    ## Runge-Kutta-Fehlberg 78 method
    rk78f = list(ID = "rk78f",
        varstep = TRUE,
        FSAL    = FALSE,
        A  = matrix(
         c(rep(0,12),
         2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0,
         0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0,
         -25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0,
         31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0,
         2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0,
         -91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0,
         2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0,
         3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0,
         -1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1
        ), nrow=13, ncol=12, byrow = TRUE),
        b1 = c(41/840, 0,0,0,0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840, 0, 0),
        b2 = c(0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840),
        c  = c(0, 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1),
        stage = 13,
        Qerr  = 7
    ),

    ## KS -> ThPe: these are implicit methods 

    ## Radau order 3
    irk3r = list(ID = "irk3r",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c(5/12, -1/12,
             3/4, 1/4),
             nrow=2, ncol=2, byrow = TRUE),
      b1 = c(3/4, 1/4) ,
      c = c(1/3, 1/4),
      stage = 2,
      Qerr = 3
    ),
    
    ## Radau IIA order 5
    irk5r = list(ID = "irk5r",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c((88-7*sqrt(6))/360, (296-169*sqrt(6))/1800,(-2+3*sqrt(6))/225,
             (296+169*sqrt(6))/1800, (88+7*sqrt(6))/360,(-2-3*sqrt(6))/225,
             (16-sqrt(6))/36, (16+sqrt(6))/36,1/9),
             nrow=3, ncol=3, byrow = TRUE),
      b1 = c((16-sqrt(6))/36,(16+sqrt(6))/36, 1/9), 
      c = c(0.4-sqrt(6)/10, 0.4+sqrt(6)/10, 1),
      stage = 3,
      Qerr = 6
    ),

    ## Hammer - Hollingsworth coefficients , order 4
    irk4hh = list(ID = "irk4hh",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c(1/4, 1/4-sqrt(3)/6,
             1/4+sqrt(3)/6, 1/4),
             nrow=2, ncol=2, byrow = TRUE),
      b1 = c(1/2, 1/2), 
      c = c(0.5-sqrt(3)/6, 0.5+sqrt(3)/6),
      stage = 2,
      Qerr = 4
    ),

    ## Kuntzmann and Butcher order 6
    irk6kb = list(ID = "irk6kb",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c(5/36, 2/9-sqrt(15)/15, 5/36 - sqrt(15)/30,
             5/36+sqrt(15)/24,2/9,5/36-sqrt(15)/24,
             5/36+sqrt(15)/30,2/9+sqrt(15)/15,5/36),
             nrow=3, ncol=3, byrow = TRUE),
      b1 = c(5/18, 4/9, 5/18), 
      c = c(1/2-sqrt(15)/10, 1/2, 1/2+sqrt(15)/10),
      stage = 3,
      Qerr = 6
    ),

    ## Lobatto order 4
    irk4l = list(ID = "irk4l",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c(0,  0,  0,
                   1/4,1/4,0,
                   0,  1,  0),
                   nrow=3, ncol=3, byrow = TRUE),
      b1 = c(1/6, 2/3, 1/6) ,
      c  = c(0,   1/2, 1),
      stage = 3,
      Qerr = 4
    ),

    ## Lobatto order 6
    irk6l = list(ID = "irk6l",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c(0,0,0,0,
             (5+sqrt(5))/60,  1/6,   (15-7*sqrt(5))/60,0,
              (5-sqrt(5))/60, (15+7*sqrt(5))/60,1/6,   0,
              1/6,(5-sqrt(5))/12, (5+sqrt(5))/12   ,   0),
              nrow=4, ncol=4, byrow = TRUE),
      b1 = c(1/12, 5/12, 5/12, 1/12) ,
      c = c(0,(5-sqrt(5))/10,(5+sqrt(5))/10,1),
      stage = 4,
      Qerr = 6
    )
  )
  ## look if the method is known; ode23 and ode45 are used as synonyms
  knownMethods <- c(lapply(methods,"[[", "ID"), "ode23", "ode45")

  if (!is.null(method)) {
    method <- unlist(match.arg(method, knownMethods))
    if (method == "ode23")
      method <- "rk23bs"
    else if (method == "ode45")
      method <- "rk45dp7"

    out <- methods[[method]]
  } else {
    out <- vector("list", 0)
  }

  ## modify a known or add a completely new method)
  ldots <- list(...)
  out[names(ldots)] <- ldots

  ## return the IDs of the methods if called with an empty argument list
  if (is.null(method) & length(ldots) == 0) {
    out <- as.vector(unlist(knownMethods))
  } else {
    ## check size consistency of parameter sets
    sl    <- lapply(out, length)
    stage <- out$stage
    if (is.matrix(out$A)) {
      if (nrow(out$A) != stage | ncol(out$A)  < stage -1 | ncol(out$A) > stage)
        stop("Size of matrix A does not match stage")
    } else {
      if (length(out$A) != stage) stop("Size of A does not match stage")
    }
    if (stage != sl$b1 | stage != sl$c)
      stop("Wrong rkMethod, length of parameters do not match")
    if (out$varstep & is.null(out$b2))
      stop("Variable stepsize method needs non-empty b2")
    if (!is.null(out$b2))
      if (sl$b2 != stage)
        stop("Wrong rkMethod, length of b2 must be empty or equal to stage")
    if (!is.null(out[["d"]])) # exact argument matching!
      if (sl[["d"]] != stage)
        stop("Wrong rkMethod, length of d must be empty or equal to stage")
    class(out) <- c("list", "rkMethod")
  }

  out
}
