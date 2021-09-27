      SUBROUTINE  MGSPAC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     MGSPAC                                                           C
C                                                                      C
C     DESCRIPTION                                                      C
C         MANAGE SPACE STEPS FOR DIFFERENTIAL EQUATION OF CONTAMINATE  C
C         TRANSPORT IN THE GEOSPHERE                                   C
C                                                                      C
C     DESCRIPTION OF ARGUMENTS                                         C
C         INPUT                                                        C
C           MXNSPT - MAXIMUM NUMBER OF SPACE POINTS                    C
C           NLAYER - NUMBER OF GEOSPHERE LAYERS                        C
C           XPATH(N) - GEOSPHERE PATH LENGTH OF EACH LAYER (M)         C
C           DISPC(N) - DISPERSION LENGTH OF EACH LAYER (M)             C
C                                                                      C
C         OUTPUT                                                       C
C           DXL(N) - SPACE LENGTH OF EACH LAYER (M)                    C
C           NPX(N) - NUMBER OF SPACE POINTS OF EACH LAYER              C
C           NXTOT  - TOTAL NUMBER OF SPACE POINTS                      C
C                                                                      C
C     CALLED BY MAIN                                                   C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INCLUDE  'par.h'
      INCLUDE  'cparam.h'
      INCLUDE  'cfarf.h'
      INCLUDE 'cvaria.h'
      INCLUDE 'var.h'
C
C
      NXTIN = 0
      DO  10  N = 1,NLAYER
          DXL(N) = 1.5*DISPC(N)
          XPP = XPATH(N)
          IF(N.EQ.NLAYER)  THEN
              XPP = 1.5*XPP
          ENDIF
          NPX(N) = XPP/DXL(N)  
          DXL(N) = XPP/NPX(N)
          NXTIN = NXTIN + NPX(N)
   10 CONTINUE
C

      IF( NXTIN+1.GT.MXNSPT )  THEN
          SPREL = FLOAT(MXNSPT-1)/FLOAT(NXTIN)
          NXTIN = 0
          DO  20  N = 1,NLAYER
              NPX(N) = NPX(N)*SPREL
              DXL(N) = XPATH(N)/NPX(N)
              IF(N.EQ.NLAYER)  THEN
                  DXL(N) = 1.5*XPATH(N)/NPX(N)
              ENDIF
              NXTIN = NXTIN + NPX(N)
   20     CONTINUE
      ENDIF
C
C
      NN      = XPATH(NLAYER)/DXL(NLAYER)
      AZUR    = XPATH(NLAYER) - DXL(NLAYER)*NN
      AZUR    = AZUR/NN
      DXL(NLAYER) = DXL(NLAYER) + AZUR
      NXTOT = NXTIN + 1
      RETURN
 6001 FORMAT(' -- SPACE -- N,XPATH,NPX,DXL: ',I2,1PE11.3,I3,E11.3)
      END
