C--------------------------------------------------------------------
C IDE of dimension 17 - describes motion of a simple
C wheelset on a rail track 
C--------------------------------------------------------------------
      SUBROUTINE wheelres(X,Y,YPRIME,CJ,DELTA,IERR,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    	Integer M
	    parameter (N=17)

      DIMENSION Y(N), DELTA(N), YPRIME(N),RPAR(2),IPAR(2)
      DOUBLE PRECISION MR, LI1, LI2,PHI,SIP,COP
C
      MR = 16.08D0
      LI1 = 0.0605D0
      LI2 = 0.366D0      
      PHI = Y(5)
      SIP = DSIN(PHI)
      COP = DCOS(PHI)
C      
      CALL wheelF(N,T,Y,DELTA,IPAR,RPAR,IERR)
C
      DO I=1,5
         DELTA(I) = YPRIME(I) - DELTA(I)
      ENDDO
      DO I=6,8
         DELTA(I) = MR*YPRIME(I) -DELTA(I)
      ENDDO
      DELTA(9)  = LI2*COP*YPRIME(9) - DELTA(9)
      DELTA(10) = LI2*YPRIME(10) - DELTA(10)

      DELTA(11) = LI1*SIP*YPRIME(9) +LI1*YPRIME(11) - DELTA(11)
      DO I=12,N
         DELTA(I) = -DELTA(I)
      ENDDO
C
      RETURN
      END
C--------------------------------------------------------------------
      subroutine wheelF(NEQN,T,Y,DF,IPAR,RPAR,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer neqn,ierr,ipar(2)
      double precision t,y(neqn),df(neqn),rpar(2)
c
c     uses the routines
c
c          reswhs
c          wheelp
c          railp
c          creep
c          constm
c
      integer ires
      double precision y12,y13
c
c we interchange y(12)-y(16) and y(13)-y(17)
c
      y12=y(12)
      y(12)=y(16)
      y(16)=y12
      y13=y(13)
      y(13)=y(17)
      y(17)=y13
c
c     index-2 formulation: ipar(1)=0
c
      ipar(1)=0
      call reswhs(t,y,y,df,ires,rpar,ipar)
c
c we interchange y(12)-y(16) and y(13)-y(17)
c
      y12=y(12)
      y(12)=y(16)
      y(16)=y12
      y13=y(13)
      y(13)=y(17)
      y(17)=y13
c
      end
C-----------------------------------------------------------------------------
	
      SUBROUTINE RESWHS ( T, Y, DY, DELTA, IRES, RPAR, IPAR )
C====================================================================
C
C     RESWHS
C     ======
C
C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
C
C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
C
C     VERSION :        NOV. 1995
C     AUTHORS :        SIMEON, FUEHRER, RENTROP
C     PURPOSE :        RESIDUAL OF EQS. OF MOTION FOR
C                      SIMULATION OF WHEELSET WITH DASSL
C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
C                      SURV. MATH. IND. I: 1-37 (1991)
C     SUBROUTINES:     WHEELP   WHEEL PROFILE (CONE)
C                      RAILP     RAIL PROFILE  (CIRCLE)
C                      CREEP   CREEP FORCES (ANALYTICAL)
C                      CONSTM  CONSTRAINT JACOBIAN ("G"-MATRIX)
C
C====================================================================
C
C     PARAMETERS  (I INPUT , O OUTPUT)
C     ==========
C
C       T            SIMULATION TIME                    (I)
C       Y(1:17)      DEPENDENT VARIABLES IN THE ORDER   (I)
C                    ( P, V, BETA, LAMBDA, Q), CF. DESCRIPTION
C       DY(1:17)     DERIVATIVES OF Y                   (I)
C       RPAR         SYSTEM PARAMETERS (DUMMY)
C       IPAR(1)      SWITCH FOR THE CONSTRAINT          (I)
C                    (=0: INDEX-2 CONSTR., =1: INDEX-3 CONSTR.)
C       DELTA(1:17)  RESIDUALS OF EQS. OF MOTION        (O)
C       IRES         ERROR FLAG: 0 = THINGS WENT FINE   (O)
C                               -1 = THE INTEGRATION MUST BE
C                                    TERMINATED DUE TO AN ERROR
C                                    (TYPICALLY A DERAILMENT OF
C                                     THE WHEELSET)
C
C       THE NOTATION OF INTERNAL VARIABLES IS CONFORMAL TO THE
C       DESCRIPTION OF THE PROBLEM IN THE TESTSET, EXCEPT FOR
C       THE DYNAMIC VARIABLE BETA WHICH IS HERE CALLED BETAP.
C
C==============================================================
C
      DOUBLE PRECISION
     *           T, Y(1:17), DY(1:17), DELTA(1:17), RPAR
      INTEGER    IRES, IPAR(*)
      DOUBLE PRECISION
     *           XX, YY, ZZ, XXP, YYP, ZZP, XXPP, YYPP, ZZPP,
     *           TETA, PHI, TETAP, PHIP, TETAPP, PHIPP, BETAP, BETAPP,
     *           E1, E2, PSIL, PSIR, XRL, XRR, XGL, XGR, FNL, FNR,
     *           RXL, RXR, DRXL, DRXR, D2RXL, D2RXR, D3RXL, D3RXR,
     *           GXL, GXR, DGXL, DGXR, D2GXL, D2GXR, D3GXL, D3GXR,
     *           TL(1:3), TR(1:3), QL(1:5), QR(1:5),
     *           DELTAL, DELTAR, DETER,
     *           SIT, COT, SIP, COP, SIA, COA, SISL, COSL, SISR, COSR,
     *           SIDL,CODL,SIDR,CODR,W1,W2,SE,S0,
     *           ALPHA, S, MR, G, V, FA, OMEGA, RN0, LA, LI1, LI2,
     *           MA, HA, FQ1, FQ2, FQ3, LM1, LM2, LM3, TOL, MU ,TG,
     *           CX,CZ, XL
      INTEGER    IERR, I
      PARAMETER( TOL   = 0.00000001d0,
     *           MR    = 16.08d0,
     *           G     = 9.81000d0,
     *           V     = 30.000d0,
     *           RN0   = 0.10000d0,
     *           LI1   = 0.0605d0,
     *           LI2   = 0.366d0,
     *           MA    = 0.d0,
     *           HA    = 0.2d0,
C                FRICTION COEFFICIENT
     *           MU    = 0.120000d0 ,
     *           XL    = 0.19d0,
     *           CX   = 6400.d0,
     *           CZ   = 6400.d0 )
C
C =====================================================================
C
C     COORDINATES
C
      XX      = Y(1)
      YY      = Y(2)
      ZZ      = Y(3)
      TETA    = Y(4)
      PHI     = Y(5)
      XXP     = Y(6)
      YYP     = Y(7)
      ZZP     = Y(8)
      TETAP   = Y(9)
      PHIP    = Y(10)
      BETAP   = Y(11)
C     SCALE LAGRANGE MULTIPLIERS
      E1      = Y(12)*1.d+4
      E2      = Y(13)*1.d+4
      PSIL    = Y(14)
      XRL     = Y(15)
      PSIR    = Y(16)
      XRR     = Y(17)
C
      XXPP    = DY(6)
      YYPP    = DY(7)
      ZZPP    = DY(8)
      TETAPP  = DY(9)
      PHIPP   = DY(10)
      BETAPP  = DY(11)
C
C    STRAIGHT TRACK, NO ADDITIONAL PROPULSION FORCES
C
      TG      = 0.D0
      S0      = 0.D0
      SE      = 0.D0
      FA      = 0.D0
      LA      = 0.D0
      S       =  0.0D0
      ALPHA   =  0.0D0
C
      OMEGA = V/RN0
      SIT   = SIN(TETA)
      COT   = COS(TETA)
      SIP   = SIN(PHI)
      COP   = COS(PHI)
      SIA   = SIN(ALPHA)
      COA   = COS(ALPHA)
      SISL  = SIN(PSIL)
      COSL  = COS(PSIL)
      SISR  = SIN(PSIR)
      COSR  = COS(PSIR)
C
      IERR  = 0
      IRES  = 0
C
C    EVALUATION OF PROFILE FUNCTIONS
C    WHEEL, LEFT SIDE
      CALL WHEELP ( XRL, RXL, DRXL, D2RXL, D3RXL, IERR )
      IF (IERR .NE. 0) GO TO 999
      XGL   = XX + XRL*COT*COP + RXL*(COT*SIP*COSL - SIT*SISL)
C    RAIL,  LEFT SIDE
      CALL RAILP   ( XGL, GXL, DGXL, D2GXL, D3GXL, IERR )
      IF (IERR .NE. 0) GO TO 999
C    WHEEL, RIGHT SIDE
      CALL WHEELP ( XRR, RXR, DRXR, D2RXR, D3RXR, IERR )
      IF (IERR .NE. 0) GO TO 999
      XGR   = XX + XRR*COT*COP + RXR*(COT*SIP*COSR - SIT*SISR)
C    RAIL,  RIGHT SIDE
      CALL RAILP   ( XGR, GXR, DGXR, D2GXR, D3GXR, IERR )
      IF (IERR .NE. 0) GO TO 999
C
C     BUILD UP CONSTRAINT JACOBIAN MATRIX "G"
C     LEFT CONSTRAINT
      CALL CONSTM (XRL, RXL, DGXL,
     *             SIT, COT, SIP, COP, SISL, COSL, QL)
C     RIGHT CONSTRAINT
      CALL CONSTM (XRR, RXR,   DGXR,
     *             SIT, COT, SIP, COP, SISR, COSR, QR)
C
C     ANGLE OF CONTACT PLANE
C
      W1     = (DRXL*COP - SIP*COSL)*COT + SISL*SIT
      W2     = -DRXL*SIP - COSL*COP
      DELTAL = ATAN( W1/W2 )
      W1     = (DRXR*COP - SIP*COSR)*COT + SISR*SIT
      W2     =  DRXR*SIP + COSR*COP
      DELTAR = ATAN( W1/W2 )
      SIDL   = SIN(DELTAL)
      CODL   = COS(DELTAL)
      SIDR   = SIN(DELTAR)
      CODR   = COS(DELTAR)
C
C     NORMAL FORCES N(P,Q,LAMBDA)
C
      DETER    = -SIDL*CODR - SIDR*CODL
      IF (ABS(DETER) .LT. TOL)  THEN
          IERR = - 20
          GO TO 999
      END IF
      W1     = QL(1)*E1 + QR(1)*E2
      W2     = QL(2)*E1 + QR(2)*E2
      FNL    = ( CODR*W1 - SIDR*W2) / DETER
      FNR    = (-CODL*W1 - SIDL*W2) / DETER
C
C     CREEPAGE FORCES
C
      CALL CREEP ( Y, FNL, FNR, V, S, ALPHA, OMEGA,
     *              RXL, RXR, DRXL, DRXR, D2RXL, D2RXR,
     *              DGXL, DGXR, D2GXL, D2GXR,
     *              DELTAL, DELTAR, MU, TL, TR, IERR )
C
C     FORCES OF CHASSIS
C
      FQ1    = MA*G*( V*V*S/G - TAN(ALPHA) ) / COA
      FQ2    = -MA*G*COA*(V*V*S*TAN(ALPHA)/G + 1.0d0)
      FQ3    = -2.0d0*CZ*ZZ
      LM1    = 0.0d0
      LM2    = -2.0D0*XL**2*CZ*TETA
      LM3    = -HA*FQ1
C
C --------------------------------------------------------------------
C
C     KINEMATICS
C
      DELTA(1) = XXP 
      DELTA(2) = YYP 
      DELTA(3) = ZZP 
      DELTA(4) = TETAP 
      DELTA(5) = PHIP
C
C     DYNAMICS: NEWTON'S LAW
C
      DELTA(6) = MR*(V*V*S*COA*(1.0d0 + (XX*COA-YY*SIA)*S)
     *            + 2.0D0*V*S*COA*ZZP)
     *            + TL(1) + TR(1) + FQ1 - MR*G*SIA
     *            + QL(1)*E1 + QR(1)*E2 - 2.0d0*CX*XX
      DELTA(7) = MR*( -V*V*S*SIA*(1.0d0 + (XX*COA-YY*SIA)*S)
     *            -  2.0D0*V*S*SIA*ZZP)
     *            + TL(2) + TR(2) + FQ2 - MR*G*COA
     *            + QL(2)*E1 + QR(2)*E2
      DELTA(8)  = MR*( - 2.0d0*V*S*(XXP*COA-YYP*SIA)
     *                + V*V*S*S*ZZ)
     *            + TL(3) + TR(3) + FA + FQ3
     *            + QL(3)*E1 + QR(3)*E2
C
C     DYNAMICS: EULER'S LAW
C
      W1     = -(XRL*SIT+RXL*SISL*COT*COP)*TL(1) - RXL*SISL*SIP*TL(2)
     *            + (-XRL*COT+RXL*SISL*SIT*COP)*TL(3)
      W2     = -(XRR*SIT+RXR*SISR*COT*COP)*TR(1) - RXR*SISR*SIP*TR(2)
     *            + (-XRR*COT+RXR*SISR*SIT*COP)*TR(3)
      DELTA(9) =  - LI2*(  - TETAP*PHIP*SIP +
     *                    V*S*(PHIP*(SIA*COT*COP + COA*SIP)
     *                         - TETAP*SIA*SIT*SIP)        )
     *            - LI1*(OMEGA+BETAP)*(PHIP - V*S*SIT*SIA)
     *         - (LI1 - LI2)*(TETAP*SIP - V*S*(COT*COP*SIA + SIP*COA))
     *                         *(PHIP - V*S*SIT*SIA)
     *            + W1 + W2 + COP*LM2 - COT*SIP*LM1 + SIT*SIP*LM3
     *            + QL(4)*E1 + QR(4)*E2
      W1        = -(XRL*COT*SIP-RXL*COSL*COT*COP)*TL(1)
     *            + (XRL*COP+RXL*COSL*SIP)*TL(2)
     *            + (XRL*SIT*SIP-RXL*COSL*SIT*COP)*TL(3)
      W2        = -(XRR*COT*SIP-RXR*COSR*COT*COP)*TR(1)
     *            + (XRR*COP+RXR*COSR*SIP)*TR(2)
     *            + (XRR*SIT*SIP-RXR*COSR*SIT*COP)*TR(3)
      DELTA(10) = - LI2*( - TETAP*V*S*SIA*COT)
     *    + LI1*(OMEGA+BETAP)*(TETAP*COP + V*S*(COT*SIP*SIA-COP*COA))
     *         + (LI1-LI2)*(TETAP*SIP - V*S*(COT*COP*SIA + SIP*COA))
     *                    *(TETAP*COP + V*S*(COT*SIP*SIA - COP*COA))
     *         + W1 + W2 + LM3 + QL(5)*E1 + QR(5)*E2
      W1     = -RXL*(COSL*SIT+SISL*COT*SIP)*TL(1) + RXL*SISL*COP*TL(2)
     *            - RXL*(COSL*COT-SISL*SIT*SIP)*TL(3)
      W2     = -RXR*(COSR*SIT+SISR*COT*SIP)*TR(1) + RXR*SISR*COP*TR(2)
     *            - RXR*(COSR*COT-SISR*SIT*SIP)*TR(3)
      DELTA(11) =  - LI1*( TETAP*PHIP*COP
     *                   - V*S*(PHIP*(COA*COP-SIA*COT*SIP)
     *                          - TETAP*SIA*SIT*COP       ) )
     *            + W1 + W2 + SIP*LM2 + LA
C
C     CONSTRAINT EQUATIONS
C
C     CONTACT CONDITION "G_1(P,Q)"
C
      DO 50 I = 12,13
         DELTA(I) = 0.0D0
  50  CONTINUE
      IF (IPAR(1) .EQ.1)THEN
C
C     INDEX-3 FORMULATION
         DELTA(12) = GXL - YY - XRL*SIP + RXL*COP*COSL
         DELTA(13) = GXR - YY - XRR*SIP + RXR*COP*COSR
      ELSE
C
C     INDEX-2 FORMULATION
         DO 60 I = 1,5
            DELTA(12) = DELTA(12) + QL(I)*Y(5+I)
            DELTA(13) = DELTA(13) + QR(I)*Y(5+I)
  60     CONTINUE
      END IF
C
C     ADDITIONAL INDEX- 1 EQUATIONS "G_2(P,Q)"
C     (NORMAL VECTORS OF CONTACT PLANE ARE PARALLEL, NONINTERSECTION CONDITION)
C
      DELTA(14) = DGXL*(DRXL*SIP + COP*COSL) + DRXL*COT*COP
     *             - COT*SIP*COSL + SIT*SISL
      DELTA(15) = DRXL*SIT*COP - SIT*SIP*COSL - COT*SISL
      DELTA(16) = DGXR*(DRXR*SIP + COP*COSR) + DRXR*COT*COP
     *             - COT*SIP*COSR + SIT*SISR
      DELTA(17) = DRXR*SIT*COP - SIT*SIP*COSR - COT*SISR
C
C     ERROR HANDLING
C
 999  IF (IERR .LT. 0) THEN
         IRES = -1
      ELSE
         IRES = 0
      END IF
C
C --------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE CONSTM (XR, RX, DGX,
     *                   SIT, COT, SIP, COP, SIPS, COPS, Q)
      DOUBLE PRECISION XR, RX, DGX, SIT, COT, SIP, COP, SIPS, COPS, Q(5)
C=======================================================================
C
C     constm
C     ======
C
C     Part of the test set for Initial Value Problems, see
C
C        http://www.cwi.nl/cwi/projects/IVPtestset.html
C
C     Version :        Nov. 1995
C     Authors :        SIMEON, FUEHRER, RENTROP
C     Purpose :        Computation of constraint Jacobian for
C                      simulation of wheelset with DASSL
C     Description:     See the above mentioned testset
C     Reference:       B. SIMEON, C. FUEHRER, P. RENTROP:
C                      Introduction to differential-algebraic
C                      equations in vehicle system dynamics,
C                      Surv. Math. Ind. I: 1-37 (1991)
C
C     PARAMETER  (I=input, O=output)
C
C     XR                      displacement xi                       (I)
c     rx                      wheel profile                         (I)
C     DGX                     derivative of rail profile            (I)
C     SIT                     SIN(TETA),TETA=Y(4)                   (I)
C     COT                     COS(TETA)                             (I)
C     SIP                     SIN(PHI),PHI=Y(5)                     (I)
C     COP                     COS(PHI)                              (I)
C     SIPS                    SIN(PSIL/R),PSIL/R=Y(14/16)           (I)
C     COPS                    COS(PSIL/R)                           (I)
C     Q(1:5)                  constraint matrix (left/right)        (O)
C
C
C     constraint matrix
C=======================================================================
C
      Q(1) = DGX
      Q(2) = -1.0D0
      Q(3) = 0.0D0
      Q(4) = DGX*(RX*(-COT*SIPS - COPS*SIP*SIT) - COP*SIT*XR)
      Q(5) = -(COPS*RX*SIP) - COP*XR + DGX*(COP*COPS*COT*RX -
     *       COT*SIP*XR)
C
      RETURN
      END



      SUBROUTINE CREEP ( Y, FNL, FNR, V, S, ALPHA, OMEGA,
     *                    RXL, RXR, DRXL, DRXR, D2RXL, D2RXR,
     *                    DGXL, DGXR, D2GXL, D2GXR,
     *                    DELTAL, DELTAR, MU, TL, TR, IERR )
C====================================================================
C
C     CREEP
C     =====
C
C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
C
C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
C
C     VERSION :        NOV. 1995
C     AUTHORS :        SIMEON, FUEHRER, RENTROP
C     PURPOSE :        COMPUTATION OF CREEPAGE FORCES FOR
C                      SIMULATION OF WHEELSET WITH DASSL
C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
C                      SURV. MATH. IND. I: 1-37 (1991)
C
C     APPROXIMATION OF CREEPAGE FORCES DUE TO A. JASCHINSKI (DLR),
C     IN ORDER TO GUARANTEE SMOOTHNESS OF THE FUNCTIONS.
C     FOR THIS END SEE:
C     A. JASCHINSKI, DLR-REPORT DLR-FB 90-06, KÖLN, 1990
C
C     FOR THE DESCRIPTION OF THE PARAMETERS SEE RESRAD.
C
C ==================================================================
C
      DOUBLE PRECISION
     *           Y(1:17), FNL, FNR, V, S, ALPHA, OMEGA,
     *           RXL, RXR, DRXL, DRXR, D2RXL, D2RXR,
     *           DGXL, DGXR, D2GXL, D2GXR,
     *           DELTAL, DELTAR, MU, TL(1:3), TR(1:3)
      INTEGER    IERR
      DOUBLE PRECISION
     *           CL, CR, RHOG, RHOR, RR, A, B, E, G, SIGMA, GM, PI,
     *           C11, C22, C23, WVR, WVK1, WVK2, WVK3,WVK4,WABS,
     *           WNUX, WNUY, WPHIS, TX, TY, XX, YY, BETAP, XRL, XRR,
     *           ZZ,ZZP, WVK5, WVK6,
     *           SITL,SITR, COTL,COTR, SIP, COP, SIA, COA, SISL, COSL,
     *           SISR, COSR,
     *           SIDL, CODL, SIDR, CODR, WT, SQT,
     *           XXP, YYP, TETAP, PHIP
      PARAMETER( E     = 1.3537956D0,
     *           G     = 0.7115218D0,
     *           SIGMA = 0.28D0,
     *           GM    = 7.92D+10,
     *           PI    = 3.1415926D0,
     *           C11   = 4.72772197D0,
     *           C22   = 4.27526987D0,
     *           C23   = 1.97203505D0 )
C
      IERR  = 0
      XX    = Y(1)
      YY    = Y(2)
      ZZ    = Y(3)
C
      SIP   = SIN(Y(5))
      COP   = COS(Y(5))
C
      XXP   = Y(6)
      YYP   = Y(7)
      ZZP   = Y(8)
      TETAP = Y(9)
      PHIP  = Y(10)
C
      BETAP = Y(11)
C
      SISL  = SIN(Y(14))
      COSL  = COS(Y(14))
C
      XRL   = Y(15)
      SISR  = SIN(Y(16))
      COSR  = COS(Y(16))
      XRR   = Y(17)
C
      SIA   = SIN(ALPHA)
      COA   = COS(ALPHA)
      SIDL  = SIN(DELTAL)
      CODL  = COS(DELTAL)
      SIDR  = SIN(DELTAR)
      CODR  = COS(DELTAR)
      SITL   = SIN(Y(4)/CODL)
      COTL   = COS(Y(4)/CODL)
      SITR   = SIN(Y(4)/CODR)
      COTR   = COS(Y(4)/CODR)

C
C     CONTACT ELLIPSES
C
      RR    = RXL*SQRT(1.D0 + DRXL*DRXL)
      RHOG  = -D2GXL/(1.D0 + DGXL*DGXL)**1.5d0
      RHOR  =  D2RXL/(1.D0 + DRXL*DRXL)**1.5d0
      A     = 0.5D0/RR
      B     = 0.5D0*(RHOG + RHOR)
      WABS  = ABS(FNL)*3.D0
      CL    = ((WABS*(1.D0-SIGMA)*E) /
     *        (2.0D0*PI*(A+B)*GM*SQRT(G)))**(1.0d0/3.0d0)
      RR    = RXR*SQRT(1.D0 + DRXR*DRXR)
      RHOG  = -D2GXR/(1.D0 + DGXR*DGXR)**1.5d0
      RHOR  =  D2RXR/(1.D0 + DRXR*DRXR)**1.5d0
      A     = 0.5D0/RR
      B     = 0.5D0*(RHOG + RHOR)
      WABS  = ABS(FNR)*3.D0
      CR    = ((WABS*(1.D0-SIGMA)*E) /
     *        (2.0D0*PI*(A+B)*GM*SQRT(G)))**(1.0d0/3.0d0)
C
C     CREEPAGE LEFT CONTACT POINT
C
C     RELATIVE VELOCITY
      WVK1 = -(OMEGA+BETAP)*RXL*(SITL*COSL+COTL*SIP*SISL)
     *       + V*S*COA*( RXL*(SITL*SIP*COSL+COTL*SISL)
     *                   + XRL*SITL*COP   - ZZ             )
     *       + XXP - TETAP*(RXL*(SITL*SIP*COSL+COTL*SISL)
     *       + XRL*SITL*COP)
     *            - PHIP*COTL*(XRL*SIP - RXL*COP*COSL)
      WVK2 = (OMEGA+BETAP)*RXL*COP*SISL
     *       + V*S*SIA*(ZZ -XRL*SITL*COP
     *       - RXL*(SITL*SIP*COSL+COTL*SISL) )
     *       + YYP + PHIP*(XRL*COP + RXL*SIP*COSL)
      WVK3 = -(OMEGA+BETAP)*RXL*(COTL*COSL-SITL*SIP*SISL) + V + ZZP
     *       + V*S * ( XX*COA - YY*SIA
     *                + COA*(RXL*(COTL*SIP*COSL-SITL*SISL)+XRL*COTL*COP)
     *                + SIA*(RXL*COP*COSL - XRL*SIP)                 )
     *       - TETAP*(XRL*COTL*COP + RXL*(COTL*SIP*COSL-SITL*SISL))
     *       + PHIP*SITL*(XRL*SIP - RXL*COP*COSL)
C
C     ROLLING VELOCITY
C
      WVK4 = WVK1 - 2.D0*XXP + 2.D0*V*S*ZZ*COA
      WVK5 = WVK2 - 2.D0*YYP - 2.D0*V*S*ZZ*SIA
      WVK6 = WVK3 - 2.D0*ZZP - 2.D0*V*S*(XX*COA-YY*SIA) - 2.D0*V
      WVR  = 0.5d0*SQRT( WVK4*WVK4 + WVK5*WVK5 + WVK6*WVK6 )
C
C     CREEPAGE
C
      WNUX = ( SITL*WVK1 + COTL*WVK3 ) / WVR
      WNUY = ( COTL*CODL*WVK1 + SIDL*WVK2 - SITL*CODL*WVK3 ) / WVR
      WPHIS= (-SIDL*( OMEGA+BETAP - V*S*SIA )
     *        +CODL*( TETAP   - V*S*COA ) ) / WVR
C
C     CREEPAGE FORCES
C
      WT   =  MU*FNL
      TX   = -WT*TANH( GM*CL*CL*C11*WNUX/WT )
      TY   = -WT*TANH( GM*CL*CL*C22*WNUY/WT + GM*CL*CL*CL*C23*WPHIS/WT)
      SQT  =  SQRT ( TX*TX + TY*TY )
      IF ( SQT*SQT .GT. WT*WT ) THEN
C         NORMALIZE
          TX   = TX * ABS(WT) / SQT
          TY   = TY * ABS(WT) / SQT
          IERR = 3
      END IF
C
C     TRANSFORMATION TO NOMINAL REFERENCE FRAME
C
      TL(1) = SITL*TX + COTL*CODL*TY
      TL(2) = +SIDL*TY
      TL(3) = COTL*TX - SITL*CODL*TY
C
C     CREEPAGE RIGHT CONTACT POINT
C
C     RELATIVE VELOCITY
      WVK1 = -(OMEGA+BETAP)*RXR*(SITR*COSR+COTR*SIP*SISR)
     *       + V*S*COA*( RXR*(SITR*SIP*COSR+COTR*SISR)
     *                   + XRR*SITR*COP -ZZ               )
     *       + XXP - TETAP*(RXR*(SITR*SIP*COSR+COTR*SISR)
     *       + XRR*SITR*COP)
     *            - PHIP*COTR*(XRR*SIP - RXR*COP*COSR)
      WVK2 = (OMEGA+BETAP)*RXR*COP*SISR
     *       + V*S*SIA*(ZZ -XRR*SITR*COP
     *       - RXR*(SITR*SIP*COSR+COTR*SISR) )
     *       + YYP + PHIP*(XRR*COP + RXR*SIP*COSR)
      WVK3 = -(OMEGA+BETAP)*RXR*(COTR*COSR-SITR*SIP*SISR) + V +ZZP
     *       + V*S * ( XX*COA - YY*SIA
     *                + COA*(RXR*(COTR*SIP*COSR-SITR*SISR)+XRR*COTR*COP)
     *                + SIA*(RXR*COP*COSR - XRR*SIP)                 )
     *       - TETAP*(XRR*COTR*COP + RXR*(COTR*SIP*COSR-SITR*SISR))
     *       + PHIP*SITR*(XRR*SIP - RXR*COP*COSR)
C
C     ROLLING VELOCITY
C
      WVK4 = WVK1 - 2.D0*XXP + 2.D0*V*S*ZZ*COA
      WVK5 = WVK2 - 2.D0*YYP - 2.D0*V*S*ZZ*SIA
      WVK6 = WVK3 - 2.D0*ZZP - 2.D0*V*S*(XX*COA-YY*SIA) - 2.D0*V
      WVR  = 0.5d0*SQRT( WVK4*WVK4 + WVK5*WVK5 + WVK6*WVK6 )
C
C     CREEPAGE
C
      WNUX = ( SITR*WVK1 + COTR*WVK3 ) / WVR
      WNUY = ( COTR*CODR*WVK1 - SIDR*WVK2 - SITR*CODR*WVK3 ) / WVR
      WPHIS= (+SIDR*( OMEGA+BETAP - V*S*SIA )
     *        +CODR*( TETAP   - V*S*COA ) ) / WVR
C
C     CREEPAGE FORCES
C
      WT   =  MU*FNR
      TX   = -WT*TANH( GM*CR*CR*C11*WNUX/WT )
      TY   = -WT*TANH( GM*CR*CR*C22*WNUY/WT + GM*CR*CR*CR*C23*WPHIS/WT)
      SQT  =  SQRT ( TX*TX + TY*TY )
      IF ( SQT*SQT .GT. WT*WT ) THEN
C         NORMALIZE
          TX   = TX * ABS(WT) / SQT
          TY   = TY * ABS(WT) / SQT
          IERR = 4
      END IF
C
C     TRANSFORMATION TO NOMINAL REFERENCE FRAME
C
      TR(1) = SITR*TX + COTR*CODR*TY
      TR(2) = -SIDR*TY
      TR(3) = COTR*TX - SITR*CODR*TY

      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE RAILP ( X, SX, DSX, D2SX, D3SX, IERR )
C====================================================================
C
C     RAILP
C     =====
C
C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
C
C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
C
C     VERSION :        NOV. 1995
C     AUTHORS :        SIMEON, FUEHRER, RENTROP
C     PURPOSE :        EVALUATION OF PROFILE RAIL FOR
C                      SIMULATION OF WHEELSET WITH DASSL
C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
C                      SURV. MATH. IND. I: 1-37 (1991)
C
C     PARAMETER  (I=INPUT, O=OUTPUT)
C
C     X          DISPLACEMENT   (XI (LEFT OR RIGHT))              (I)
C     SX         VALUE OF RAIL PROFILE FUNKTION   R(XI)                (O)
C     DSX        ITS DERIVATIVE                                   (O)
C     D2SX       ITS SECOND DERIVATIVE                            (O)
C     D3SX       ITS THIRD DERIVATIVE (ONLY USED IN THE INDEX-1 CASE)
C     IERR       ERROR CODE                                       (O)
C                IERR =  0 : NO ERROR
C                     = -2 : OUT OF RANGE
C
C
C     CONSTANTS:
C
C     DELTA0     CONE ANGLE
C     RN0        NOMINAL RADIUS
C     AR         1/2 WHEEL DISTANCE
C     RS         RAIL RADIUS
C     EPS        TOLERANCE
C
C====================================================================
      DOUBLE PRECISION
     *           X, SX, DSX, D2SX, D3SX
      INTEGER    IERR
      DOUBLE PRECISION
     *           DELTA0, RN0, AR, RS, EPS, T, XABS, SIR
      PARAMETER( DELTA0 = 0.0262d0,
     *           RN0    = 0.1000d0,
     *           AR     = 0.1506d0,
     *           RS     = 0.06d0,
     *           EPS    = 0.00001d0 )
C
      SIR  = SIN(DELTA0)*RS
      XABS = ABS(X)
      IF ( ( (AR + SIR + RS - EPS) .LE. XABS) .OR.
     *     ( (AR + SIR - RS + EPS) .GE. XABS )    )  THEN
           IERR = -2
           WRITE(*,*) 'OUT OF RAIL PROFILE (DERAILMENT)'
      ELSE
           T    = SQRT ( RS*RS - (XABS-AR-SIR)**2 )
           SX   = T - RN0 - COS(DELTA0)*RS
           DSX  = SIGN(1.D0,X)*(-XABS+AR+SIR)/T
           D2SX = -RS*RS/(T*T*T)
           D3SX = 3*RS*RS/(T*T*T*T)*DSX
           IERR = 0
      END IF

      RETURN
      END
C
C -------------------  END  OF  SUBROUTINE  ------------------------
C
      SUBROUTINE WHEELP ( X, RX, DRX, D2RX,D3RX, IERR )
C====================================================================
C
C     WHEELP
C     ======
C
C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
C
C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
C
C     VERSION :        NOV. 1995
C     AUTHORS :        SIMEON, FUEHRER, RENTROP
C     PURPOSE :        EVALUATION OF PROFILE WHEEL FOR
C                      SIMULATION OF WHEELSET WITH DASSL
C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
C                      SURV. MATH. IND. I: 1-37 (1991)
C
C ===============================================================
C
C     PARAMETER  (I=INPUT, O=OUTPUT)
C
C     X          DISPLACEMENT   (XI (LEFT OR RIGHT))              (I)
C     RX         VALUE OF PROFILE FUNKTION   R(XI)                (O)
C     DRX        ITS DERIVATIVE                                   (O)
C     D2RX       ITS SECOND DERIVATIVE                            (O)
C     D3RX       ITS THIRD DERIVATIVE (ONLY USED IN THE INDEX-1 CASE)
C     IERR       ERROR CODE                                       (O)
C                IERR =  0 : NO ERROR
C                     = -1 : OUT OF RANGE
C
C
C     CONSTANTS:
C
C     DELTA0     CONE ANGLE
C     RN0        NOMINAL RADIUS
C     AR         1/2 WHEEL DISTANCE
C     B1         INNER WHEEL LIMIT
C     B2         OUTER WHEEL LIMIT
C
C====================================================================
C
      DOUBLE PRECISION
     *           X, RX, DRX, D2RX, D3RX
      INTEGER    IERR
      DOUBLE PRECISION
     *           DELTA0, RN0, AR, B1, B2, TA, XABS
      PARAMETER( DELTA0 = 0.0262D0,
     *           RN0    = 0.1000D0,
     *           AR     = 0.1506D0,
     *           B1     = 0.0000D0,
     *           B2     = 4.000D0 )
C
      TA   = TAN(DELTA0)
      XABS = ABS(X)
      IF ( (B1 .GE. XABS) .OR. (XABS .GE. B2) ) THEN
          IERR = -1
          WRITE(*,*) 'OUT OF WHEEL PROFILE (DERAILMENT)'
      ELSE
          RX   = RN0 + TA*(AR-XABS)
          DRX  = SIGN(1.D0,X)*(-TA)
          D2RX = 0.0d0
          D3RX = 0.0d0
          IERR = 0
      END IF

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE wheelsoln(NEQN,T,Y)
      DOUBLE PRECISION  T, Y(NEQN)
      
c
c     DASSL applied to Wheelset problem, tend = 10
c
c     ATOL = RTOL = 1d-9 for p, v, q and
c     ATOL = RTOL = 1d10 for lambda to exclude the Lagrange multipliers
c                            from error control
c
c     # steps =             58510
c     # steps accepted =    54229
c     # f-eval =            95533
c     # Jac-eval =          10218
c
      Y( 1) =  0.86355386965811D-02
      Y( 2) =  0.13038281022727D-04
      Y( 3) = -0.93635784016818D-04
      Y( 4) = -0.13642299804033D-01
      Y( 5) =  0.15292895005422D-02
      Y( 6) = -0.76985374142666D-01
      Y( 7) = -0.25151106429207D-03
      Y( 8) =  0.20541188079539D-02
      Y( 9) = -0.23904837703692D+00
      Y(10) = -0.13633468454173D-01
      Y(11) = -0.24421377661131D+00
      Y(16) = -0.10124044903201D-01
      Y(17) = -0.56285630573753D-02
      Y(14) =  0.37839614386969D-03
      Y(15) =  0.14173214964613D+00
      Y(12) = -0.33666751972196D-03
      Y(13) = -0.15949425684022D+00

      RETURN
      END










