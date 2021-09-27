      SUBROUTINE BATMAN(T,NELCHN,ALAMB,C0,RLEACH,TSLEAC,CIN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C SUBROUTINE BATMAN                                                    C
C A. SALTELLI, JOINT RESEARCH CENTRE OF ISPRA                          C
C DIVISION OF RADIOCHEMISTRY                                           C
C OCTOBER 30, 1988                                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C DESCRIPTION:                                                         C
C                                                                      C
C SUBROUTINE TO COMPUTE ACTIVITIES (MOL/KG OF CONCRETE) AFTER          C
C CHAIN DECAY USING BATEMAN EQUATIONS . SUBROUTINE FOR BATMAN          C
C EQUATIONS INCLUDING A EXPONENTIAL LEACHING RATE(RLEACH) STARTING     C
C AT 'CONTIM' TIME.                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C DICTIONARY:                                                          C
C                                                                      C
C T               = TIME PINT CONSIDERD                                C
C NELCHN          = NUMBER OF ELEMENT                                  C
C ALAMB(I)        = DECAY CONSTANT (1/A)                               C
C RLEACH(I)       = LEACHING RATE                                      C
C TSLEACH         = TIME WHEN LEACHING START                           C
C C0(I)           = NUCLIDE CONCENTRATION IN THE MATRIX AT TIME        C
C                   OF VAULT CLOSURE (MOL/KG)                          C
C CIN(I)          = NUCLIDE CONCENTRATION IN THE MATRIX AFTER DECAY    C
C                   (MOLES/KG) - IN THE VAULT  SUBROUTINE CIN ARE      C
C                   CONVERTED TO FLUX FROM THE BUFFER (MOLES/A)        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION ALAMB(*),C0(*),RLEACH(*)
      DIMENSION CIN(*)
      DIMENSION EX(10),RL(8)
C
      DOUBLE PRECISION XS,XSS,XP,XPP,EX
C
C
C
C                   INITIALISE
C                   ==========
C
      DO  5  I = 1, NELCHN
      IF( T.LT.TSLEAC )  THEN
          RL(I)=0.
      ELSE
          RL(I)=RLEACH(I)
      ENDIF
    5 CONTINUE
C
C                   BATEMAN EQ.S
C                   ============
C
      IF( T.LE.0. )  THEN
          DO  10  I = 1, NELCHN
              CIN(I)=C0(I)
   10     CONTINUE
      ELSE
          IF( NELCHN.EQ.1 )  THEN
              CIN(1)=C0(1)*EXP(-ALAMB(1)*T)*EXP(-RL(1)*(T-TSLEAC))
          ELSE
              DO  20  I = 1, NELCHN
                  EX(I)=EXP(-ALAMB(I)*T)*EXP(-RL(I)*(T-TSLEAC))
   20         CONTINUE
              DO  70  I = 1, NELCHN
                  XS=0.
                  DO  60  K = 1, I
                      XP=1.
                      IF( K.LT.I )  THEN
                          IMN1=I-1
                          DO  30  MM = K, IMN1
                              XP=XP*ALAMB(MM)
   30                     CONTINUE
                      ENDIF
                      XSS=0.
                      DO  50  NN = K, I
                          XPP=1.
                          IF( K.LT.I )  THEN
                              DO  40  JJ = K, I
                                  IF( JJ.NE.NN )  THEN
                                      XPP=XPP*(ALAMB(JJ)+RL(JJ)
     &                                       -(ALAMB(NN)+RL(NN)))
                                  ENDIF
   40                         CONTINUE
                          ENDIF
                          XSS=XSS+EX(NN)/XPP
   50                 CONTINUE
                      XS=XS+C0(K)*XP*XSS
   60             CONTINUE
                  CIN(I)=XS
   70         CONTINUE
          ENDIF
      ENDIF
      DO  80  I = 1, NELCHN
          IF( CIN(I).LT.0. )  THEN
              WRITE(6,9000) NELCHN,I,CIN(I) ,T
          ENDIF
   80 CONTINUE
      RETURN
 9000 FORMAT('***WARNING FROM BATMAN: CIN NEGAT. FOR ELEMENT
     & NUMBER ',I3,' CIN=',E12.4,' TIME=',E12.4,'****')
      END
