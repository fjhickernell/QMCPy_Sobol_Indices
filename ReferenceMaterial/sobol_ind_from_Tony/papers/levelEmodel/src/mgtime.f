      SUBROUTINE  MGTIME
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     MGSPAC                                                           C
C                                                                      C
C     DESCRIPTION                                                      C
C         MANAGE TIME STEPS FOR DIFFERENTIAL EQUATION OF CONTAMINATE   C
C         TRANSPORT IN THE GEOSPHERE                                   C
C                                                                      C
C     DESCRIPTION OF VARIABLES                                         C
C           NXTOT  - TOTAL NUMBER OF SPACE POINTS                      C
C                                                                      C
C     CALLED BY MAIN                                                   C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INCLUDE  'par.h'
      INCLUDE  'cparam.h'
      INCLUDE  'cnucld.h'
      INCLUDE  'cnearf.h'
      INCLUDE  'cfarf.h'
      INCLUDE  'ctime.h'
C
C
      XUPPER=ALOG(2.)*2.**13.
      NELCHN = NEL(NCHAIN)
C
C
C                   MAXIMUN TRAVEL TIME
C                   ===================
C
      TRTMP = 0.
      DO  20  I =1, NELCHN
          TRVL = 0.
          DO  10  N = 1,NLAYER
              RV = VREAL(N)/RET(NCHAIN,I,N)
              TRVL = TRVL + (XPATH(N)/RV)
   10     CONTINUE
          IF( TRVL.GT.TRTMP )  TRTMP = TRVL
   20 CONTINUE
      TTMAX = (TRTMP + CONTIM)*1.5
C
C
C                   DTS < TRAVEL TIME OF ONE SPACE STEP
C                   ===================================
C
      DO  40  N = 1,NLAYER
          DO  30  I = 1,NELCHN
              DTTEMP = 3./4.*DXL(N)*RET(NCHAIN,I,N)/VREAL(N)
              IF( I.EQ.1 .AND. N.EQ.1 )  THEN
                  DTS = DTTEMP
              ELSE IF( DTTEMP.LT.DTS )  THEN
                  DTS = DTTEMP
              ENDIF
   30     CONTINUE
   40 CONTINUE
C
C
C                   DTS > ALAMB
C                   ===========
C
      DO  50  I = 1,NELCHN
          BANDA = 1./ALAMB(NCHAIN,I)
          IF(DTS.GE.BANDA)  THEN
              DTS = 3.*BANDA/4.
              if(sverb) then
                  WRITE(LOUT,6002) BNAME(NCHAIN,I),DTS
              endif
          ENDIF
   50 CONTINUE
C
C
C                   TIME SERIES
C                   ===========
C
      NTIME = MXNTPT - IOUT
      TINIT = TOCC
      IF( CONTIM.GT.TOCC ) TINIT = CONTIM
      TSTEP(1) = TINIT
      DTFAC = 1./100./RLEACH(NCHAIN,1)
C
      IF( DTFAC.GE.DTS )  THEN
          DT = DTS
          DO  55  J = 2, MXNTPT - IOUT
              TSTEP(J) = TSTEP(J-1) + DT
              IF(TSTEP(J).GT.TMAX)  THEN
                  NTIME=J
                  GO TO 70
              ENDIF
   55     CONTINUE
      ELSE
          JTEMP = 0
          TSTEP(2) = TSTEP(1) + DTFAC
          DO  60  J = 3, MXNTPT - IOUT
              IF( JTEMP.EQ.0 )  THEN
                  EQU = RLEACH(NCHAIN,1)*(TSTEP(J-1) - TINIT)
                  IF( EQU.GT.XUPPER )  EQU = XUPPER
                  DTP = DTFAC*EXP(EQU)
                  IF( DTP.LT.DTS )  THEN
                      DT = DTP
                  ELSE
                      DT = DTS
                      JTEMP=J
                  ENDIF
              ELSE
                  DT = DTS
              ENDIF
              TSTEP(J) = TSTEP(J-1) + DT
              IF(TSTEP(J).GT.TMAX)  THEN
                  NTIME=J
                  GO TO 70
              ENDIF
   60     CONTINUE
      ENDIF
C
   70 CONTINUE
C
C
C                   INCLUDE TOUT(IOUT) IN THE
C                   TIME SERIESE
C                   =========================
C
      IF( IOUT.NE.0 )  THEN
          K = 1
          DO  74  I = 1, IOUT
              TP = TOUT(I)
              IF( TP.GE.TSTEP(1) )  THEN
                  DO  72  J = 1, NTIME
                      IF( TP.EQ.TSTEP(J) )  THEN
                          GO TO 74
                      ENDIF
   72             CONTINUE
                  NN = NTIME + K
                  TSTEP(NN) = TP
                  K = K + 1
              ENDIF
   74     CONTINUE
          NTIME = NN
C
          CALL  S I F T (NTIME,TSTEP)
C
          DO  78  J = 1, NTIME
              TS = TSTEP(J)
              DO  76  I = 1, IOUT
                  IF( TOUT(I).EQ.TS )  THEN
                      ITPNT(I) = J
                  ENDIF
   76         CONTINUE
   78     CONTINUE
      ENDIF
c      if(sverb) then 
c          WRITE(33,6001) TMAX,TTMAX,DTTEMP,DTS
c          WRITE(33,6003) TINIT,DTFAC,JTEMP,NTIME
c      endif
      RETURN
 6001 FORMAT(' -- TIME -- TMAX,TTMAX,DTTEMP,DTS: ',1P4E11.3)
 6002 FORMAT(' -- TIME --  DTS = 1 / ALAMB : ',A8,'DTS =',1PE11.3)
 6003 FORMAT(' -- TIME -- TINIT,DTFAC,JTEMP,NTIME: ',1P2E11.3,2I6)
      END
