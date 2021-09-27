      SUBROUTINE BIOS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     BIOS                                                             C
C                                                                      C
C     DESCRIPTION                                                      C
C         THIS BIOSPHERE MODEL IS VERY SIMPLE. THE FULL GEOSPHERE FLUX C
C         IS ASSUMED TO ENTER A STREAM WHICH IS USED FOR DRINKING      C
C         WATER.                                                       C
C                                                                      C
C     CALLED BY MAIN                                                   C
C                                                                      C
C     USES WRDOSE                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DICTIONARY                                                       C
C         DOSRT(1,N,I,K) --- INGESTION DOSE FROM N-TH CHAIN, I-TH      C
C                            ELEMENT AT K-TH TIME                      C
C         ADOS1(I,K,L)   --- GEOSPHERE FLUX                            C
C         TTC(N,I)       --- DOSE CONVERSION FACTOR OF N-TH CHAIN,     C
C                            I-TH ELEMENT                              C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INCLUDE  'par.h'
      INCLUDE  'cparam.h'
      INCLUDE  'cnucld.h'
      INCLUDE  'cbios.h'
      INCLUDE  'ctime.h'
      INCLUDE  'coutp.h'
      INCLUDE  'cspop.h'
C
C
      NELCHN = NEL(NCHAIN)
      DO  20  K = 1, NTIME
          DO  10  I = 1, NELCHN
C
C             SV/A  =  (SV/MOL)  *  MOL/A
C             DOSRT      ADOS        TTC
C
              DOSRT(1,NCHAIN,I,K) = ADOS1(I,K,NLAYER)*TTC(NCHAIN,I)

   10     CONTINUE
   20 CONTINUE
C
C
      TI=TSTEP(1)
      TF=TSTEP(NTIME)
      IF( IOUT.NE.0 )  THEN
          DO  40  IO = 1, IOUT
              T = TOUT(IO)
              IF( T.GE.TI .AND. T.LE.TF )  THEN
                  IT = ITPNT(IO)
                  DO  30  I = 1, NELCHN
                      FDOSE(NCHAIN,I,IO)=DOSRT(1,NCHAIN,I,IT)
   30             CONTINUE
              ENDIF
   40     CONTINUE
      ENDIF
C
C
      RETURN
      END
