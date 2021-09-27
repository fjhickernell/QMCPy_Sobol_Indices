      SUBROUTINE NEARF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      NEARF                                                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INCLUDE 'par.h'
      INCLUDE 'cparam.h'
      INCLUDE 'cnucld.h'
      INCLUDE 'cnearf.h'
      INCLUDE 'ctime.h'
      INCLUDE 'coutp.h'
C
      DIMENSION  AL(MXNELM), CC(MXNELM), RL(MXNELM)
      DIMENSION  CIN(MXNELM)
C
C
      NELCHN = NEL(NCHAIN)
      DO  10  I = 1, NELCHN
          AL(I)=ALAMB(NCHAIN,I)
          CC(I)=C0(NCHAIN,I)
          RL(I)=RLEACH(NCHAIN,I)
   10 CONTINUE
      TSLEAC = TOCC
      IF( CONTIM.GT.TOCC )  TSLEAC = CONTIM
C
C
      NTNF=NTIME
      KDOS=0
      DO  30  J = 1, NTIME
          T = TSTEP(J)
C
C                   BATEMAN EQUATION
C                   ================
C
          C A L L  B A T M A N (T,NELCHN,AL,CC,RL,TSLEAC,CIN)
C
          IFLAG=0
          DO  20  I = 1, NELCHN
              IF( CIN(I).NE.0. )  THEN
                  KDOS=KDOS+1
                  ADOS(I,J) = CIN(I)*RL(I)
                  IFLAG=1
              ENDIF
   20     CONTINUE
          IF( IFLAG.EQ. 0 )  THEN
              NTNF=J
              GO TO 40
          ENDIF
   30 CONTINUE
C
   40 CONTINUE
      IF(KDOS.EQ.0)  THEN
          ISAVE(NCHAIN)=0
      ENDIF
C
      RETURN
      END
