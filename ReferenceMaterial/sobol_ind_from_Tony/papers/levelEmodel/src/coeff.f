      SUBROUTINE COEFF(MXNCHN,MXNELM,NCHAIN,NELCHN,NLAYER,
     &                 DT,DXL,ALAMB,DISPH,VREAL,RET,
     &                 AF1,AF2,AF3,AF4,AF5,AF6,AF7,AF8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  SUBROUTINE COEFF :                                                  C
C  ================                                                    C
C  P. PRADO , JOINT RESEARCH CENTRE OF ISPRA                           C
C  DIVISION OF RADIOCHEMISTRY                                          C
C  MARCH / 23 / 1989                                                   C
C  CONTRACT JRC/ISPRA - ENRESA (CIEMAT) TO ADAPT THE LISA CODE         C
C  TO FRACTURATED MEDIA.                                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  DESCRIPTION :                                                       C
C  ===========                                                         C
C                COMPUTE THE COEFFICIENT VALUES FOR THE EQUATIONS      C
C  SYSTEM SOLVED INTO 'TRES SUBROUTINE'. CRANK-NICOLSON SCHEME IS      C
C  CONSIDERED TO SOLVE THE DIFFERENTIAL TRANSPORT EQUATION BY FINITE   C
C  DIFFERENCES TECNIQUE.                                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  DICTIONARY.                                                         C
C  ===========                                                         C
C   ICHAIN  = NUMBER OF CHAIN CONSIDERED                               C
C   IELM    = ELEMENT NUMBER INTO THE CHAIN                            C
C   NLY     = LAYER NUMBER INTO THE GEOSPHERE                          C
C   DT      = TIME STEP (Y)
C   DXL(LY) = SPACE STEP (M)
C             ------------------------------------                     C
C             ------------------------------------                     C
C  AF1   (NLY)      = 1ST. EQUATION COEFFICIENT (SP. POINT NOT BOUND)  C
C  AF2   (NLY)      = 2ND. EQUATION COEFFICIENT (SP. POINT NOT BOUND)  C
C  AF3   (NLY)      = 3RD. EQUATION COEFFICIENT (SP. POINT NOT BOUND)  C
C  AF4   (NLY)      = 4RD. EQUATION COEFFICIENT (SP. POINT NOT BOUND)  C
C                                                                      C
C  AF5   (NLY)      = 1ST. EQUATION COEFFICIENT (SP. POINT  BOUND)     C
C  AF6   (NLY)      = 2ND. EQUATION COEFFICIENT (SP. POINT  BOUND)     C
C  AF7   (NLY)      = 3RD. EQUATION COEFFICIENT (SP. POINT  BOUND)     C
C  AF8   (NLY)      = 4RD. EQUATION COEFFICIENT (SP. POINT  BOUND)     C
C                                                                      C
C  ALAMB (NCHAIN,I)   = NUCLIDE DECAY CONSTANT (1/Y)                   C
C  DISPH (LY)         = HYDRODINAMIC DISPERSION (M**2/Y)               C
C  RET   (NCHAIN,I,LY)=  RETENTION FACTOR                              C
C  VREAL (LY)         =  GROUNDWATER REAL VELOCITY INTO THE LAYER(M/Y) C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION DXL(*)
      DIMENSION ALAMB(MXNCHN,*),DISPH(*),VREAL(*),RET(MXNCHN,MXNELM,*)
      DIMENSION AF1(MXNELM,*),AF2(MXNELM,*),AF3(MXNELM,*),AF4(MXNELM,*)
      DIMENSION AF5(MXNELM,*),AF6(MXNELM,*),AF7(MXNELM,*),AF8(MXNELM,*)
C
C
      DO  20  I = 1, NELCHN
          DO  10  J = 1, NLAYER
              WF=DT/(RET(NCHAIN,I,J)*DXL(J)*2.)
              GF=WF*DISPH(J)/DXL(J)
              ZF=WF*VREAL(J)/2.
              AF1(I,J)=GF + ZF
              AF2(I,J)=GF - ZF
              AF3(I,J)=1. + (2.* GF + 0.5*ALAMB(NCHAIN,I)*DT)
              AF4(I,J)=1. - (2.* GF + 0.5*ALAMB(NCHAIN,I)*DT)
C
              IF( J.LT.NLAYER )  THEN
                  RT = DT/(RET(NCHAIN,I,J)+RET(NCHAIN,I,J+1))
                  A1 = RT*DISPH(J)/DXL(J)/(DXL(J)+DXL(J+1))*2.
                  A2 = RT*DISPH(J+1)/DXL(J+1)/(DXL(J)+DXL(J+1))*2.
                  B1 = RT*VREAL(J)/(2.*DXL(J))
                  B2 = RT*VREAL(J+1)/(2.*DXL(J+1))
                  AF5(I,J) = A1 + B1
                  AF6(I,J) = A2 - B2
                  AF7(I,J) = 1.+AF5(I,J)+AF6(I,J)+0.5*ALAMB(NCHAIN,I)*DT
                  AF8(I,J) = 1.-AF5(I,J)-AF6(I,J)-0.5*ALAMB(NCHAIN,I)*DT
              ENDIF
C
   10     CONTINUE
   20 CONTINUE
      RETURN
      END
