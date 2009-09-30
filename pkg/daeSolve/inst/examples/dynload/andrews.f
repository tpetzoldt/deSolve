c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is derived from the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Andrews'' squeezing mechanism (in index 3 formulation)
c        index 3 DAE of dimension 27
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/andrews.f
c
c     This is revision
c     $Id: andrews.F,v 1.2 2006/10/02 10:29:13 testset Exp $
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c the parameter initialisation function
c-----------------------------------------------------------------------

      SUBROUTINE andinit(daeparms)
      EXTERNAL daeparms
      
      double precision parms(42)
      common / myparms/ parms

        CALL daeparms(42,parms)
        
      END SUBROUTINE andinit

c-----------------------------------------------------------------------
c the model residual function
c-----------------------------------------------------------------------

      SUBROUTINE andres(T,Y,YPRIME,CJ,DELTA,IERR,RPAR,IPAR)
      
      implicit none
      DOUBLE PRECISION T, Y(27), DELTA(27), YPRIME(27),RPAR(*), CJ
      INTEGER I, J, IERR, N, IPAR(*)
      character (len=100) msg
C
      N = 27
      CALL Fand(N,T,Y,DELTA,IERR)
C
      DO J = 1,14
         DELTA(J) = YPRIME(J)-DELTA(J)
      ENDDO
      DO I=15,N
         DELTA(I) = -DELTA(I)
      ENDDO
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE Fand(NEQN,T,Y,DY,IIERR)
      implicit none

      DOUBLE PRECISION T, Y(27), DY(27)
      integer i,j,neqn,IIERR
      DOUBLE PRECISION m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom
      common/myparms/ m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom

      DOUBLE PRECISION sibe,sith,siga,siph,side,siom,siep,
     +     cobe,coth,coga,coph,code,coom,coep,
     +     sibeth,siphde,siomep,cobeth,cophde,coomep,
     +     bep,thp,php,dep,omp,epp,
     +     m(7,7),ff(7),gp(6,7),g(6),xd,yd,lang,force,fx,fy

      sibe = dsin(y(1))
      sith = dsin(y(2))
      siga = dsin(y(3))
      siph = dsin(y(4))
      side = dsin(y(5))
      siom = dsin(y(6))
      siep = dsin(y(7))
c
      cobe = dcos(y(1))
      coth = dcos(y(2))
      coga = dcos(y(3))
      coph = dcos(y(4))
      code = dcos(y(5))
      coom = dcos(y(6))
      coep = dcos(y(7))
c
      sibeth = dsin(y(1)+y(2))
      siphde = dsin(y(4)+y(5))
      siomep = dsin(y(6)+y(7))
c
      cobeth = dcos(y(1)+y(2))
      cophde = dcos(y(4)+y(5))
      coomep = dcos(y(6)+y(7))
c
      bep = y(8)
      thp = y(9)
      php = y(11)
      dep = y(12)
      omp = y(13)
      epp = y(14)
c
      do 20 j = 1,7
         do 10 i = 1,7
            m(i,j) = 0d0
 10      continue
 20   continue
c
      m(1,1) = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 + i2
      m(2,1) = m2*(da**2-da*rr*coth) + i2
      m(2,2) = m2*da**2 + i2
      m(3,3) = m3*(sa**2+sb**2) + i3
      m(4,4) = m4*(e-ea)**2 + i4
      m(5,4) = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
      m(5,5) = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)
     &     + i4 + i5
      m(6,6) = m6*(zf-fa)**2 + i6
      m(7,6) = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
      m(7,7) = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)
     &     + i6 + i7

      do 40 j=2,7
         do 30 i=1,j-1
            m(i,j) = m(j,i)
 30      continue
 40   continue
c
      xd = sd*coga + sc*siga + xb
      yd = sd*siga - sc*coga + yb
      lang  = dsqrt ((xd-xc)**2 + (yd-yc)**2)
      force = - c0 * (lang - l0)/lang
      fx = force * (xd-xc)
      fy = force * (yd-yc)
      ff(1) = mom - m2*da*rr*thp*(thp+2*bep)*sith
      ff(2) = m2*da*rr*bep**2*sith
      ff(3) = fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga)
      ff(4) = m4*zt*(e-ea)*dep**2*coph
      ff(5) = - m4*zt*(e-ea)*php*(php+2*dep)*coph
      ff(6) = - m6*u*(zf-fa)*epp**2*coom
      ff(7) = m6*u*(zf-fa)*omp*(omp+2*epp)*coom
c
      do 60 j=1,7
         do 50 i=1,6
            gp(i,j) = 0d0
 50      continue
 60   continue
c
      gp(1,1) = - rr*sibe + d*sibeth
      gp(1,2) = d*sibeth
      gp(1,3) = - ss*coga
      gp(2,1) = rr*cobe - d*cobeth
      gp(2,2) = - d*cobeth
      gp(2,3) = - ss*siga
      gp(3,1) = - rr*sibe + d*sibeth
      gp(3,2) = d*sibeth
      gp(3,4) = - e*cophde
      gp(3,5) = - e*cophde + zt*side
      gp(4,1) = rr*cobe - d*cobeth
      gp(4,2) = - d*cobeth
      gp(4,4) = - e*siphde
      gp(4,5) = - e*siphde - zt*code
      gp(5,1) = - rr*sibe + d*sibeth
      gp(5,2) = d*sibeth
      gp(5,6) = zf*siomep
      gp(5,7) = zf*siomep - u*coep
      gp(6,1) = rr*cobe - d*cobeth
      gp(6,2) = - d*cobeth
      gp(6,6) = - zf*coomep
      gp(6,7) = - zf*coomep - u*siep
c
      g(1) = rr*cobe - d*cobeth - ss*siga - xb
      g(2) = rr*sibe - d*sibeth + ss*coga - yb
      g(3) = rr*cobe - d*cobeth - e*siphde - zt*code - xa
      g(4) = rr*sibe - d*sibeth + e*cophde - zt*side - ya
      g(5) = rr*cobe - d*cobeth - zf*coomep - u*siep - xa
      g(6) = rr*sibe - d*sibeth - zf*siomep + u*coep - ya
c
      do 70 i=1,14
         dy(i) = y(i+7)
   70 continue

      do 100 i=15,21
         dy(i) = -ff(i-14)
         do 80 j=1,7
            dy(i) = dy(i)+m(i-14,j)*y(j+14)
   80    continue
         do 90 j=1,6
            dy(i) = dy(i)+gp(j,i-14)*y(j+21)
   90    continue
  100 continue
      do 110 i=22,27
         dy(i) = g(i-21)
  110 continue

      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE andjac(T,Y,YPRIME,DFDY,CON,RPAR,IPAR)
      INTEGER NEQN,MN
      PARAMETER (NEQN=27, MN = 27)
      DOUBLE PRECISION T, Y(NEQN), DFDY(MN,NEQN),CON


c-----------------------------------------------------------------------
c     the Jacobian computed here is an approximation, see p. 540 of
c     Hairer & Wanner `solving ordinary differential equations II'
c-----------------------------------------------------------------------
      DOUBLE PRECISION m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom
      common/myparms/ m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom

      DOUBLE PRECISION sibe,siga,siph,side,siom,siep,
     +     cobe,coth,coga,code,coep,
     +     sibeth,siphde,siomep,cobeth,cophde,coomep,
     +     m(7,7),gp(6,7)


      sibe = dsin(y(1))
      siga = dsin(y(3))
      siph = dsin(y(4))
      side = dsin(y(5))
      siom = dsin(y(6))
      siep = dsin(y(7))
c
      cobe = dcos(y(1))
      coth = dcos(y(2))
      coga = dcos(y(3))
      code = dcos(y(5))
      coep = dcos(y(7))
c
      sibeth = dsin(y(1)+y(2))
      siphde = dsin(y(4)+y(5))
      siomep = dsin(y(6)+y(7))
c
      cobeth = dcos(y(1)+y(2))
      cophde = dcos(y(4)+y(5))
      coomep = dcos(y(6)+y(7))
c
      do 51 j = 1,7
         do 52 i = 1,7
            m(i,j) = 0d0
 52      continue
 51      continue
c
      m(1,1) = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 +i2
      m(2,1) = m2*(da**2-da*rr*coth) + i2
      m(2,2) = m2*da**2 + i2
      m(3,3) = m3*(sa**2+sb**2) + i3
      m(4,4) = m4*(e-ea)**2 + i4
      m(5,4) = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
      m(5,5) = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)
     +         + i4 + i5
      m(6,6) = m6*(zf-fa)**2 + i6
      m(7,6) = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
      m(7,7) = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)
     +         + i6 + i7

      do 40 j=2,7
         do 30 i=1,j-1
            m(i,j) = m(j,i)
   30    continue
   40 continue
c
      do 60 j=1,7
         do 50 i=1,6
            gp(i,j) = 0d0
   50    continue
   60 continue
c
      gp(1,1) = - rr*sibe + d*sibeth
      gp(1,2) = d*sibeth
      gp(1,3) = - ss*coga
      gp(2,1) = rr*cobe - d*cobeth
      gp(2,2) = - d*cobeth
      gp(2,3) = - ss*siga
      gp(3,1) = - rr*sibe + d*sibeth
      gp(3,2) = d*sibeth
      gp(3,4) = - e*cophde
      gp(3,5) = - e*cophde + zt*side
      gp(4,1) = rr*cobe - d*cobeth
      gp(4,2) = - d*cobeth
      gp(4,4) = - e*siphde
      gp(4,5) = - e*siphde - zt*code
      gp(5,1) = - rr*sibe + d*sibeth
      gp(5,2) = d*sibeth
      gp(5,6) = zf*siomep
      gp(5,7) = zf*siomep - u*coep
      gp(6,1) = rr*cobe - d*cobeth
      gp(6,2) = - d*cobeth
      gp(6,6) = - zf*coomep
      gp(6,7) = - zf*coomep - u*siep
c
      do 80 j=1,neqn
         do 70 i=1,neqn
            dfdy(i,j) = 0d0
 70      continue
 80   continue
      do 90 i=1,14
         dfdy(i,i+7) = 1d0
 90   continue
      do 110 i=1,7
         do 100 j=1,7
            dfdy(14+j,14+i) = m(j,i)
 100     continue
 110  continue
      do 130 i=1,6
         do 120 j=1,7
            dfdy(14+j,21+i) = gp(i,j)
 120     continue
 130  continue
      do 150 i=1,7
         do 140 j=1,6
            dfdy(21+j,i) = gp(j,i)
 140     continue
 150  continue
C
      do i=1,neqn
         do  j=1,neqn
            dfdy(i,j) = -dfdy(i,j)
         enddo
      enddo
c compute pd = -df/dy + con*M
      do j=1,14
         dfdy(j,j) = 1.0d0/con+dfdy(j,j)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
