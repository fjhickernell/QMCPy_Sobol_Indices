      SUBROUTINE  GTM1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  SUBROUTINE GTM1 (FAR-FIELD MODULE)                                  C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INCLUDE 'par.h'
      INCLUDE 'cparam.h'
      INCLUDE 'cnucld.h'
      INCLUDE 'cfarf.h'
      INCLUDE 'cbios.h'
      INCLUDE 'ctime.h'
      INCLUDE 'coutp.h'
C
      DIMENSION AF1(MXNELM,MXNLYS),AF2(MXNELM,MXNLYS)
      DIMENSION AF3(MXNELM,MXNLYS),AF4(MXNELM,MXNLYS)
      DIMENSION AF5(MXNELM,MXNLYS),AF6(MXNELM,MXNLYS)
      DIMENSION AF7(MXNELM,MXNLYS),AF8(MXNELM,MXNLYS)
      DIMENSION CBEF(MXNELM,MXNSPT),CATEMP(MXNSPT)
      DIMENSION CBB(MXNSPT),DECAF(MXNSPT),CPAR(MXNSPT)
      DIMENSION ISIGN(MXNELM)
C
C
C                   TSV = TRUNCATION DOSE (SV)
C                   ==========================
C
      TSV = 1.E-15
      NELCHN = NEL(NCHAIN)
      NPXL1 = NXTOT - 1
      NPXL2 = NXTOT - 2

C
C
C                   INITIALIZE
C                   ==========
C
      DO  20  J = 1, NPXL2
          DO  10  I = 1, NELCHN
              CBEF(I,J) = 0.
   10     CONTINUE
   20 CONTINUE
      NXF=0
      DXF=0.
      DO  22  N = 1, NLAYER-1
          NXF=NXF + NPX(N)
          DXF=DXF + NPX(N)*DXL(N)
   22 CONTINUE
      NNN=XPATH(NLAYER)/DXL(NLAYER)
      NXF=NXF + NNN
C
C
C                   TIME  LOOP
C                   ===========
C
      DTB=0.
      NTFF = NTIME
      DO  110  J = 2, NTIME
          TNEW=TSTEP(J)
          DT=TNEW-TSTEP(J-1)
C
C                   COMPUTE SYSTEM COEFFICIENTS
C                   ===========================
C
          IF( DT.NE.DTB )  THEN
C
              C A L L  C O E F F (MXNCHN,MXNELM,NCHAIN,NELCHN,NLAYER,
     &                            DT,DXL,ALAMB,DISPH,VREAL,RET,
     &                            AF1,AF2,AF3,AF4,AF5,AF6,AF7,AF8)
C
              DTB=DT
          ENDIF
C
C
C                   ELEMENT LOOP
C                   ============
C
          DO  70  I = 1, NELCHN
C
C                   BOUNDARY VALUES AT
C                   BEFORE AND CURRENT TIME
C                   =======================
C
              CB = ADOS(I,J-1)
              CA = ADOS(I,J)
C
C                   INTERNAL VALUES
C                   AT BEFORE TIME
C                   ==============
C
              DO  30  IX = 1, NPXL2
                  CBB(IX) = CBEF(I,IX)
   30         CONTINUE
C
C                   DAUGHTER INCREASE
C                   DURING A TIME STEP
C                   ==================
C

              IF( I.EQ.1 )  THEN
                  DO  40  IX = 1,NPXL2
                      DECAF(IX) = 0.
   40             CONTINUE
              ELSE
                  LY  = 1
                  MPL = NPX(LY)
                  DO  50  IX = 1, NPXL2
                      DECAF(IX) = RET(NCHAIN,I-1,LY)/RET(NCHAIN,I,LY)
     &                            *ALAMB(NCHAIN,I-1)*CPAR(IX)*DT
                      IF (IX .EQ. MPL) THEN
                          LY  = LY + 1
                          MPL = MPL + NPX(LY)
                      ENDIF
   50             CONTINUE
              ENDIF
C
C
C                   CALL TRES(TRID.EQ.SYS.)
C                   =======================
C
              C A L L  T R E S (MXNELM,I,NLAYER,NPXL2,NPX,
     &                          CA,CB,CBB,DECAF,
     &                          AF1,AF2,AF3,AF4,AF5,AF6,AF7,AF8,CATEMP)
C
C     WRITE(6,9100)  NCHAIN,I,J,TSTEP(J)
C     WRITE(6,9101) (CATEMP(IK),IK=1,NPXL2)
 9100 FORMAT(' DEBUG **',3I5,1PE15.3)
 9101 FORMAT(1H ,1P,10E10.3)
C                   SAVE OUTPUTS FOR NEXT ITERATION
C                   ===============================
C

              DO  60  IK = 1, NPXL2
                 CPAR(IK)= CATEMP(IK)
                 CBEF(I,IK)= CATEMP(IK)
      IF( CATEMP(IK).LT.0. )  THEN
 9998     FORMAT(' ***** WARNING *****',3I10,1P,E15.3)
 9999     FORMAT(1H ,1P,10E10.3)
          ENDIF
   60         CONTINUE
C
C
C                   END OF ELEMENT LOOP
C                   ===================
C
   70     CONTINUE
C
C                   SAVE OUTPUT FOR EACH LAYER
C                   ==========================
C
          ITSIGN = 1
          DO  90  I = 1, NELCHN
              DO  80  K = 1, NLAYER
                  NN = NPX(K)
                  IF(K.EQ.NLAYER)  THEN
                      NN = NXF
                  ENDIF
                  ADOS1(I,J,K) = CBEF(I,NN)

   80         CONTINUE
              IF( ADOS1(I,J,NLAYER).GE.ADOS1(I,J-1,NLAYER) )  THEN
                  ISIGN(I) = 0
              ELSE
                  ISIGN(I) = 1
              ENDIF
              ITSIGN = ITSIGN*ISIGN(I)
   90     CONTINUE
          IF( ITSIGN.EQ.1 )  THEN
              TD = 0.
              DO  100  I = 1, NELCHN
                  TD = TD + ADOS1(I,J,NLAYER)*TTC(NCHAIN,I)
  100         CONTINUE
              IF( TD.LT.TSV )  THEN
                  NTFF = J
                  GO TO 120
              ENDIF
         ENDIF
  110 CONTINUE
C
  120 CONTINUE
C
      RETURN
      END
