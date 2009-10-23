c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is derived from the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        NAND gate
c        index 0 IDE of dimension 14
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/nand.f
c
c     This is revision
c     $Id: nand.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine initgate (daeparms)
      EXTERNAL daeparms

      Double precision parms (14)
      common/myparms/parms
      CALL daeparms(14,parms)
      
      END SUBROUTINE
      
c-----------------------------------------------------------------------
      subroutine gateres(t,y,yprime,CJ,f,IERR,RPAR,IPAR)
      integer neqn,ierr,ipar(*)
      parameter(neqn = 14)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      integer i,j
      double precision am(14,14),fy(14),dum

      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI

      call   CAP(14,y,AM)
      call   FCN(14,t,y,fy,ierr)

      if(ierr.eq.-1)return

      do 20 i=1,14
         dum = -fy(i)
         do 10 j=1,14
            dum = dum+AM(i,j)*yprime(j)
   10    continue
         f(i) = dum
   20 continue

      return
      end

c-----------------------------------------------------------------------
      subroutine nandsoln(neqn,t,y)
      integer neqn
      double precision t,y(14)
c
c computed at Cray C90, using Cray double precision:
C Solving NAND gate using PSIDE
C
C User input:
C
C give relative error tolerance: 1d-16
C give absolute error tolerance: 1d-16
C
C
C Integration characteristics:
C
C    number of integration steps       22083
C    number of accepted steps          21506
C    number of f evaluations          308562
C    number of Jacobian evaluations      337
C    number of LU decompositions       10532
C
C CPU-time used:                         451.71 sec

      y(  1) =  0.4971088699385777d+1
      y(  2) =  0.4999752103929311d+1
      y(  3) = -0.2499998781491227d+1
      y(  4) = -0.2499999999999975d+1
      y(  5) =  0.4970837023296724d+1
      y(  6) = -0.2091214032073855d+0
      y(  7) =  0.4970593243278363d+1
      y(  8) = -0.2500077409198803d+1
      y(  9) = -0.2499998781491227d+1
      y( 10) = -0.2090289583878100d+0
      y( 11) = -0.2399999999966269d-3
      y( 12) = -0.2091214032073855d+0
      y( 13) = -0.2499999999999991d+1
      y( 14) = -0.2500077409198803d+1

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE FCN(N,T,Y,F,ierr)
C ---------------------------------------------------------------------
C
C Right-hand side f(X,t) for the network equation
C             C(Y) * Y' - f(Y,t) = 0
C describing the nand gate
C
C ---------------------------------------------------------------------
C
C Input parameters:
C          N......number of node potentials (14)
C          T......time point t
C          Y......node potentials at time point t
C Output parameter:
C          F......right-hand side f(Y,t)
C
C External reference:
C          IDS: Drain-source current
C          IBS: Nonlinear current characteristic for diode between
C               bulk and source
C          IBD: Nonlinear current characteristic for diode between
C               bulk and drain
C          PULSE: Input signal in pulse form
C
C ---------------------------------------------------------------------

      INTEGER N,ierr
      double precision T,Y(N),F(N) ,IDS,IBS,IBD, V1,V2,V1D,V2D
      EXTERNAL IDS, IBS, IBD, PULSE

      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI

      CALL PULSE(T,V1,V1D,0.D0,5.D0,5.D0,5.D0,5.D0,5.D0,20.D0)
      CALL PULSE(T,V2,V2D,0.D0,5.D0,15.D0,5.D0,15.D0,5.D0,40.D0)


      F(1)=-(Y(1)-Y(5))/RGS-IDS(1,Y(2)-Y(1),Y(5)-Y(1),Y(3)-Y(5),
     *       Y(5)-Y(2),Y(4)-VDD,ierr)
      F(2)=-(Y(2)-VDD)/RGD+IDS(1,Y(2)-Y(1),Y(5)-Y(1),Y(3)-Y(5),
     *       Y(5)-Y(2),Y(4)-VDD,ierr)
      F(3)=-(Y(3)-VBB)/RBS + IBS(Y(3)-Y(5))
      F(4)=-(Y(4)-VBB)/RBD + IBD(Y(4)-VDD)
      F(5)=-(Y(5)-Y(1))/RGS-IBS(Y(3)-Y(5))-(Y(5)-Y(7))/RGD-
     *       IBD(Y(9)-Y(5))
      F(6)=CGS*V1D-(Y(6)-Y(10))/RGS-
     *       IDS(2,Y(7)-Y(6),V1-Y(6),Y(8)-Y(10),V1-Y(7),Y(9)-Y(5),ierr)
      F(7)=CGD*V1D-(Y(7)-Y(5))/RGD+
     *       IDS(2,Y(7)-Y(6),V1-Y(6),Y(8)-Y(10),V1-Y(7),Y(9)-Y(5),ierr)
      F(8)=-(Y(8)-VBB)/RBS + IBS(Y(8)-Y(10))
      F(9)=-(Y(9)-VBB)/RBD + IBD(Y(9)-Y(5))
      F(10)=-(Y(10)-Y(6))/RGS-IBS(Y(8)-Y(10))-
     *         (Y(10)-Y(12))/RGD-IBD(Y(14)-Y(10))
      F(11)=CGS*V2D-Y(11)/RGS-IDS(2,Y(12)-Y(11),V2-Y(11),Y(13),
     *       V2-Y(12),Y(14)-Y(10),ierr)
      F(12)=CGD*V2D-(Y(12)-Y(10))/RGD+
     *       IDS(2,Y(12)-Y(11),V2-Y(11),Y(13),V2-Y(12),Y(14)-Y(10),ierr)
      F(13)=-(Y(13)-VBB)/RBS + IBS(Y(13))
      F(14)=-(Y(14)-VBB)/RBD + IBD(Y(14)-Y(10))

      if(ierr.eq.-1)return
      RETURN
      END

      double precision FUNCTION IDS (NED,VDS, VGS, VBS, VGD, VBD, ierr)
C ---------------------------------------------------------------------------
C
C Function evaluating the drain-current due to the model of
C Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   NED  Integer parameter for MOSFET-type
C   VDS  Voltage between drain and source
C   VGS  Voltage between gate and source
C   VBS  Voltage between bulk and source
C   VGD  Voltage between gate and drain
C   VBD  Voltage between bulk and drain
C
C External reference:
C   GDSP, GDSM Drain function for VDS > 0 gevalp. VDS < 0
C
C ---------------------------------------------------------------------------

      INTEGER NED,ierr
      double precision VDS, VGS, VBS, VGD, VBD,GDSP, GDSM
      EXTERNAL GDSP, GDSM

      IF ( VDS .GT. 0.D0 ) THEN
       IDS = GDSP (NED,VDS, VGS, VBS,ierr)
      ELSE IF ( VDS .EQ. 0.D0) THEN
       IDS = 0.D0
      ELSE IF ( VDS .LT. 0.D0) THE N
       IDS = GDSM (NED,VDS, VGD, VBD,ierr)
      END IF

      if(ierr.eq.-1)return

      RETURN
      END

      double precision FUNCTION GDSP (NED,VDS, VGS, VBS, ierr)
      integer NED,ierr
      double precision  VDS, VGS, VBS,VTE

      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI
C

      IF(NED.EQ.1) THEN
C --- Depletion-type
      VT0=-2.43D0
      CGAMMA=.2D0
      PHI=1.28D0
      BETA=5.35D-4
      ELSE
C --- Enhancement-type
      VT0=.2D0
      CGAMMA=0.035D0
      PHI=1.01D0
      BETA=1.748D-3
      END IF

      if(phi-vbs.lt.0d0.or.phi.lt.0d0)then
         ierr=-1
         return
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBS) - DSQRT(PHI) )

      IF ( VGS-VTE .LE. 0.D0) THEN
       GDSP = 0.D0
      ELSE IF ( 0.D0 .LT. VGS-VTE .AND. VGS-VTE .LE. VDS ) THEN
       GDSP = - BETA * (VGS - VTE)**2.D0 * (1.D0 + DELTA*VDS)
      ELSE IF ( 0.D0 .LT. VDS .AND. VDS .LT. VGS-VTE ) THEN
       GDSP = - BETA * VDS * (2.D0*(VGS - VTE) - VDS) *
     *          (1.D0 + DELTA*VDS)
      END IF

      RETURN
      END

      double precision FUNCTION GDSM (NED,VDS, VGD, VBD, ierr)
      integer NED,ierr
      double precision VDS, VGD, VBD,VTE


      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI

      IF(NED.EQ.1) THEN
C --- Depletion-type
      VT0=-2.43D0
      CGAMMA=.2D0
      PHI=1.28D0
      BETA=5.35D-4
      ELSE
C --- Enhancement-type
      VT0=.2D0
      CGAMMA=0.035D0
      PHI=1.01D0
      BETA=1.748D-4
      END IF

      if(phi-vbd.lt.0d0.or.phi.lt.0d0)then
         ierr=-1
         return
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBD) - DSQRT(PHI) )

      IF ( VGD-VTE .LE. 0.D0) THEN
       GDSM = 0.D0
      ELSE IF ( 0.D0 .LT. VGD-VTE .AND. VGD-VTE .LE. -VDS ) THEN
       GDSM = BETA * (VGD - VTE)**2d0 * (1.D0 - DELTA*VDS)
      ELSE IF ( 0.D0 .LT. -VDS .AND. -VDS .LT. VGD-VTE ) THEN
       GDSM = - BETA * VDS * (2d0 *(VGD - VTE) + VDS) *
     *          (1.D0 - DELTA*VDS)
      END IF

      RETURN
      END


      double precision FUNCTION IBS (VBS)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C source due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and source
C
C ---------------------------------------------------------------------------

      double precision VBS
      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI

C

C
C     IBS = GBS (VBS)
C

      IF ( VBS .LE. 0.D0 ) THEN
       IBS = - CURIS * ( DEXP( VBS/VTH ) - 1.D0 )
      ELSE
       IBS = 0.D0
      END IF

      RETURN
      END

      double precision FUNCTION IBD (VBD)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C drain  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and drain
C
C ---------------------------------------------------------------------------

      double precision VBD
      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI


C

C
C     IBD = GBD (VBD)
C
      IF ( VBD .LE. 0.D0 ) THEN
       IBD = - CURIS * ( DEXP( VBD/VTH ) - 1.D0 )
      ELSE
       IBD = 0.D0
      END IF
      RETURN
      END

      SUBROUTINE PULSE(X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD)
C ---------------------------------------------------------------------------
C
C Evaluating input signal at time point X
C
C Structure of input signal:
C
C                -----------------------                       HIGH
C               /                       \
C              /                         \
C             /                           \
C            /                             \
C           /                               \
C          /                                 \
C         /                                   \
C        /                                     \
C  ------                                       ---------      LOW
C
C |DELAY|   T1  |         T2           |   T3  |
C |          P     E     R     I     O     D            |
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   X                      Time-point at which input signal is evaluated
C   LOW                    Low-level of input signal
C   HIGH                   High-level of input signal
C   DELAY,T1,T2,T3, PERIOD Parameters to specify signal structure
C
C Output parameter:
C   VIN    Voltage of input signal at time point X
C   VIND   Derivative of VIN at time point X
C
C ---------------------------------------------------------------------------

      double precision X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD,TIME

      TIME = DMOD(X,PERIOD)

      IF (TIME.GT.(DELAY+T1+T2+T3)) THEN
      VIN = LOW
      VIND= 0.D0
      ELSE IF (TIME.GT.(DELAY+T1+T2)) THEN
      VIN = ((HIGH-LOW)/T3)*(DELAY+T1+T2+T3-TIME) + LOW
      VIND= -((HIGH-LOW)/T3)
      ELSE IF (TIME.GT.(DELAY+T1)) THEN
      VIN = HIGH
      VIND= 0.D0
      ELSE IF (TIME.GT.DELAY) THEN
      VIN = ((HIGH-LOW)/T1)*(TIME-DELAY) + LOW
      VIND= ((HIGH-LOW)/T1)
      ELSE
      VIN = LOW
      VIND=0.D0
      END IF

      RETURN
      END

      SUBROUTINE CAP(N,Y,AM)
C ---------------------------------------------------------------------
C
C Voltage-dependent capacitance matrix C(Y) for the network equation
C             C(Y) * Y' - f(Y,t) = 0
C describing the nand gate
C
C ---------------------------------------------------------------------
C
C Input parameters:
C          N......number of node potentials (14)
C          Y......value of node potentials
C Output parameter:
C          AM.....voltage-dependent capacitance matrix
C
C External reference:
C          CBDBS: Voltage-dependent capacitance CBS(V) and CBD(V)
C
C ---------------------------------------------------------------------
      double precision CBDBS
      INTEGER N
      double precision Y(N), AM(N,N)
      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI

      EXTERNAL CBDBS
      integer I,J

      DO 10 I=1,N
        DO 20 J=1,N
          AM(I,J)=0d0
 20     CONTINUE
 10   CONTINUE

      AM(1,1)=CGS
      AM(1,5)=-CGS
      AM(2,2)=CGD
      AM(2,5)=-CGD
      AM(3,3)=CBDBS(Y(3)-Y(5))
      AM(3,5)=-CBDBS(Y(3)-Y(5))
      AM(4,4)=CBDBS(Y(4)-VDD)
      AM(5,1)=-CGS
      AM(5,2)=-CGD
      AM(5,3)=-CBDBS(Y(3)-Y(5))
      AM(5,5)=CGS+CGD-AM(5,3)+
     *          CBDBS(Y(9)-Y(5))+C9
      AM(5,9)=-CBDBS(Y(9)-Y(5))
      AM(6,6)=CGS
      AM(7,7)=CGD
      AM(8,8)=CBDBS(Y(8)-Y(10))
      AM(8,10)=-CBDBS(Y(8)-Y(10))
      AM(9,5)=-CBDBS(Y(9)-Y(5))
      AM(9,9)=CBDBS(Y(9)-Y(5))
      AM(10,8)=-CBDBS(Y(8)-Y(10))
      AM(10,10)=-AM(8,10)+CBDBS(Y(14)-Y(10))+C9
      AM(10,14)=-CBDBS(Y(14)-Y(10))
      AM(11,11)=CGS
      AM(12,12)=CGD
      AM(13,13)=CBDBS(Y(13))
      AM(14,10)=-CBDBS(Y(14)-Y(10))
      AM(14,14)=CBDBS(Y(14)-Y(10))

      RETURN
      END
      
      double precision FUNCTION CBDBS (V)
C ---------------------------------------------------------------------------
C
C Function evaluating the voltage-dependent capacitance between bulk and
C drain gevalp. source  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   V    Voltage between bulk and drain gevalp. source
C
C ---------------------------------------------------------------------------
      double precision V,PHIB

      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      double precision VT0, BETA, CGAMMA,  PHI
     *
      COMMON /PARMS/   RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 DELTA, CURIS, VTH,  VDD, VBB
      COMMON /CONST/   VT0, BETA, CGAMMA,  PHI

      PHIB=0.87D0

      IF ( V .LE. 0.D0 ) THEN
       CBDBS = CBD/DSQRT(1.D0-V/PHIB)
      ELSE
       CBDBS = CBD*(1.D0+V/(2.D0*PHIB))
      END IF

      RETURN
      END
