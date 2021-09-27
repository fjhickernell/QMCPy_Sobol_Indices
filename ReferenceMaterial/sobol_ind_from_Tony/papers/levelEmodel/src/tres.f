      SUBROUTINE TRES(MXNELM,I,NLAYER,NPXL2,NPX,CA,CB,CBEF,DECAF,
     &                AF1,AF2,AF3,AF4,AF5,AF6,AF7,AF8,CATEMP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  SUBROUTINE TRES       ( TRIDIAGONAL EQUATION SYSTEM )               C
C  ================                                                    C
C                                                                      C
C  P. PRADO , JOINT RESEARCH CENTER OF ISPRA                           C
C  DIVISION OF RADIOCHEMISTRY                                          C
C  MARCH / 23 / 1989                                                   C
C  CONTRACT JRC/ISPRA - ENRESA (CIEMAT) TO ADAPT THE LISA CODE         C
C  TO FRACTURATED MEDIA.                                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  DESCRIPTION :                                                       C
C  ===========                                                         C
C  SUBROUTINE TO SET UP THE DIFFERENTIAL EQUATION SYSTEM AT ONE TIME   C
C  POINT FOR ALL INTERNAL  SPACE  POINTS  AND  FOR ALL NUCLIDES INTO   C
C  THE   DECAY   CHAIN  CONSIDERED. THE  IMPLICIT  METHOD IS THE       C
C  CRANK-NICOLSON SCHEME AND THE TRIDIAGONAL EQUATIONS                 C
C  SYSTEM IS SOLVED USING GAUSS'S ELIMINATION METHOD (WITHOUT PIVOTING)C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  DICTIONARY .                                                        C
C  ==========                                                          C
C   NLY     = NUMBER OF LAYERS CONSIDERED                              C
C   NPX(NLY)= NUMBER OF GRID SPACE IN EACH LAYER                       C
C   NPXL2   = TOTAL NUMBER OF INTERNAL GRID POINTS                     C
C   DT      = TIME STEP (Y)                                           C
C   CA               = NUCLIDE CONCENTRATION FOR X=0 AT TIME POINT    C
C                       UNDER COMPUTE.                                 C
C   CB               = NUCLIDE CONCENTRATION FOR X=0 AT BEFORE        C
C                       TIME POINT                                     C
C   CBEF  (IX)       = NUCLIDE CONCENTRATION AT BEFORE ITERATION.     C
C   DECAF (IX)       = DAUGTHER INCREASE FROM FATHER DECAY.           C
C             ------------------------------------                     C
C             ------------------------------------                     C
C   AF1   (NLY)      = 1ST. EQUATION COEFFICIENT (SP. POINT NOT BOUND)C
C   AF2   (NLY)      = 2ND. EQUATION COEFFICIENT (SP. POINT NOT BOUND)C
C   AF3   (NLY)      = 3RD. EQUATION COEFFICIENT (SP. POINT NOT BOUND)C
C   AF4   (NLY)      = 4RD. EQUATION COEFFICIENT (SP. POINT NOT BOUND)C
C                                                                      C
C   AF5   (NLY)      = 1ST. EQUATION COEFFICIENT (SP. POINT  BOUND)   C
C   AF6   (NLY)      = 2ND. EQUATION COEFFICIENT (SP. POINT  BOUND)   C
C   AF7   (NLY)      = 3RD. EQUATION COEFFICIENT (SP. POINT  BOUND)   C
C   AF8   (NLY)      = 4RD. EQUATION COEFFICIENT (SP. POINT  BOUND)   C
C   ALFA  (IX)       = INTERNAL VARIABLE USED TO SOLVE THE           C
C                       TRIDIAGONAL EQUATION SYSTEM.                  C
C   BETA  (IX)       = INTERNAL VARIABLE USED TO SOLVE THE           C
C                       TRIDIAGONAL EQUATION SYSTEM.                  C
C   GAMA  (IX)       = INTERNAL VARIABLE USED TO SOLVE THE           C
C                       TRIDIAGONAL EQUATION SYSTEM.                  C
C   CATEMP(IX)       = OUTPUT OF THIS ITERLATION                       C
C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PARAMETER ( MXNSPT =500 )
C     IMPLICIT  REAL*8(A-H,O-Z)
C
      DIMENSION NPX(*),CBEF(*),DECAF(*)
      DIMENSION AF1(MXNELM,*),AF2(MXNELM,*)
      DIMENSION AF3(MXNELM,*),AF4(MXNELM,*)
      DIMENSION AF5(MXNELM,*),AF6(MXNELM,*)
      DIMENSION AF7(MXNELM,*),AF8(MXNELM,*)
      DIMENSION CATEMP(*)
C
      DIMENSION ALFA(MXNSPT),BETA(MXNSPT),GAMA(MXNSPT)
C
C
      NPXL3 = NPXL2 - 1
C
C                   GAUSS ELIMINATION METHOD
C                   =========================
C
C                             FIRST INTERNAL POINT
C                             ====================
      LY = 1
      MPL = NPX(LY)
      ALFA(1) = AF3(I,LY)
      BETA(1) = AF1(I,LY)*(CA + CB) +
     &          AF4(I,LY)*CBEF(1) + AF2(I,LY)*CBEF(2) +
     &          DECAF(1)
      GAMA(1) = BETA(1)
C
C
      DO  10  IX = 2, NPXL3
          IF( IX.LT.MPL )  THEN
C
C                             INTERNAL POINTS
C                             ===============
C
              BETA(IX) = AF1(I,LY)*CBEF(IX-1) +
     &                   AF4(I,LY)*CBEF(IX) + AF2(I,LY)*CBEF(IX+1) +
     &                   DECAF(IX)
              TATE     = AF1(I,LY) / ALFA(IX-1)
              ALFA(IX) = AF3(I,LY) - (TATE*AF2(I,LY))
              GAMA(IX) = BETA(IX) + (TATE*GAMA(IX-1))
C
          ELSE IF( IX.EQ.MPL )  THEN
C
C                             INTERFACE POINT
C                             ===============
C
              BETA(IX) = AF5(I,LY)*CBEF(IX-1) +
     &                   AF8(I,LY)*CBEF(IX) + AF6(I,LY)*CBEF(IX+1) +
     &                   DECAF(IX)
              TATE     = AF5(I,LY) / ALFA(IX-1)
              ALFA(IX) = AF7(I,LY) - (TATE*AF2(I,LY))
              GAMA(IX) = BETA(IX) + (TATE*GAMA(IX-1))
C
              LY = LY + 1
          ELSE IF( IX.EQ.MPL+1 )  THEN
C
C                             FIRST INTERNAL POINT OF NEXT LAYER
C                             ==================================
C
              BETA(IX) = AF1(I,LY)*CBEF(IX-1) +
     &                   AF4(I,LY)*CBEF(IX) + AF2(I,LY)*CBEF(IX+1) +
     &                   DECAF(IX)
              TATE     = AF1(I,LY) / ALFA(IX-1)
              ALFA(IX) = AF3(I,LY) - (TATE*AF6(I,LY-1))
              GAMA(IX) = BETA(IX) + (TATE*GAMA(IX-1))
C
              MPL = MPL + NPX(LY)
          ENDIF
   10 CONTINUE
C
C                             LAST INTERNAL POINT
C                             BOUND: C(N) = C(N-1)
C                             ====================
C
      BETA(NPXL2) = AF1(I,LY)*CBEF(NPXL3) +
     &              (AF4(I,LY) + AF2(I,LY))*CBEF(NPXL2) +
     &              DECAF(NPXL2)
      TATE  = AF1(I,LY)/ALFA(NPXL3)
      BOUND = AF3(I,LY) - AF2(I,LY)
      ALFA(NPXL2) = BOUND - (TATE*AF2(I,LY))
      GAMA(NPXL2) = BETA(NPXL2) + (TATE*GAMA(NPXL3))
C
C                             SOLVE SYSTEM EQUATION
C                             =====================
C
      CATEMP(NPXL2) = GAMA(NPXL2)/ALFA(NPXL2)
      LY    = NLAYER
      NBACK = NPXL2 - NPX(LY) + 1
      DO  20  IX = NPXL3, 1, -1
          IF( IX.EQ.NBACK )  THEN
              LY    = LY - 1
              ANTER = GAMA(IX) + (CATEMP(IX+1)*AF6(I,LY))
              NBACK = NBACK - NPX(LY)
          ELSE
              ANTER = GAMA(IX) + (CATEMP(IX+1)*AF2(I,LY))
          ENDIF
          CATEMP(IX) = ANTER/ALFA(IX)
   20 CONTINUE
C
      RETURN
      END
