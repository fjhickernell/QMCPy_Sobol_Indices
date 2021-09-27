      SUBROUTINE SETZER
C
      INCLUDE 'par.h'
      INCLUDE 'cparam.h'
      INCLUDE 'cnucld.h'
      INCLUDE 'ctime.h'
      INCLUDE 'coutp.h'
      INCLUDE 'cspop.h'
C
      NELCHN = NEL(NCHAIN)
      DO  40  K = 1, MXNTPT
          TSTEP(K) = 0.
          DO  30  I = 1, NELCHN
              ADOS(I,K) = 0.
              DO  10  L = 1, MXNLYS
                  ADOS1(I,K,L) = 0.
   10         CONTINUE
              DO  20  JPTH = 1, NPTH
                  DOSRT(JPTH,NCHAIN,I,K)=0.
   20         CONTINUE
   30     CONTINUE
   40 CONTINUE
C
C
      DO  50  J = 1, IOUT
          CFLUX(J) = 0.
          AFLUX(J) = 0.
          CDOSE(NCHAIN,J) = 0.
   50 CONTINUE
      DO  70  I = 1, NELCHN
          DO  60  J = 1, IOUT
              FDOSE(NCHAIN,I,J) = 0.
   60     CONTINUE
   70 CONTINUE
      RETURN
      END
