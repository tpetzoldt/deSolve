c----------------------------------------------------------------------
c initialising the parameter common block
c----------------------------------------------------------------------
      SUBROUTINE carpar(daeparms)

      EXTERNAL daeparms
      double precision parms(7)
      common /myparms/parms

      call daeparms(7, parms)
      return
      end

c----------------------------------------------------------------------
c the residual function
c----------------------------------------------------------------------
      SUBROUTINE carres (T,Y,YPRIME,CJ,DELTA,ier,FOUT,IPAR)
      IMPLICIT NONE
      INTEGER N, IPAR(3), IER, I
      PARAMETER (N=10)
      DOUBLE PRECISION T,Y(N), DELTA(N),YPRIME(N), FOUT(10)

      double precision eps, m,  L, L0, r, w, g
      common / myparms/  eps, m, L, L0, R, W, G

      double precision k, cj, lam1, lam2
      double precision yl, yr, yb, xl, xr, xb, Lr, Ll
C
      if (ipar(1) < 10) call rexit("nout should be at least 10")
      k = m*eps*eps/2d0

      yb  = r*sin(w*t)
      xb  = sqrt(L*L-yb*yb)

      do 10 i=1,4
         FOUT(i) = y(i+4)
   10 continue

      xl   = y(1)
      yl   = y(2)
      xr   = y(3)
      yr   = y(4)
      lam1 = y(9)
      lam2 = y(10)

      Ll = sqrt(xl**2+yl**2)
      Lr = sqrt((xr-xb)**2+(yr-yb)**2)

      FOUT(5)  =(L0-Ll)*xl/Ll +lam1*xb+2d0*lam2*(xl-xr)
      FOUT(6)  =(L0-Ll)*yl/Ll +lam1*yb+2d0*lam2*(yl-yr)-k*g
      FOUT(7)  =(L0-Lr)*(xr-xb)/Lr -2d0*lam2*(xl-xr)
      FOUT(8)  =(L0-Lr)*(yr-yb)/Lr -2d0*lam2*(yl-yr)-k*g

      FOUT(9)  = xb*xl+yb*yl
      FOUT(10) = (xl-xr)**2+(yl-yr)**2-L*L
      do i=1,4
         delta(i) =yprime(i)-FOUT(i)
      enddo
      do i=5,8
         delta(i) = k*yprime(i)- FOUT(i)
      enddo
      do i=9,10
         delta(i) = -FOUT(i)
      enddo

      return
      end
C-----------------------------------------------------------------------


