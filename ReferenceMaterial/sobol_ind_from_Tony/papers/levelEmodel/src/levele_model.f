      SUBROUTINE levele_model
            
      REAL*4  ELAPS,USRTIM,SYSTIM
      
      INTEGER L,I
      INCLUDE 'par.h'
      INCLUDE 'cparam.h'
      INCLUDE 'cvaria.h'
      INCLUDE 'cname.h'
      INCLUDE 'cnucld.h'
      INCLUDE 'cnearf.h'
      INCLUDE 'cfarf.h'
      INCLUDE 'ctime.h'
      INCLUDE 'cbios.h'
      INCLUDE 'coutp.h'
      INCLUDE 'cspop.h'
      INCLUDE 'const.h'
      INCLUDE 'var.h'
      INCLUDE 'rel.h'

	DO   I = 1, NDC
	  NELCHN = NEL(I)
	  DO   J = 1, NELCHN
	    AM(I,J)=ALAMB(I,J)*6.02252E+23/3.1536E+07
	  enddo	  
	enddo
		
C                   CHECK IF ANY CHAIN IS EMPTY
C                    => KSAVE(NCHAIN)=0
C                   ===========================
C
      DO  20  N = 1, NDC
          NELCHN = NEL(N)
          KSAVE(N) = 0
          DO  10  I = 1, NELCHN
              IF( C0(N,I).NE.0. )  THEN
                  KSAVE(N) = 1
              ENDIF
   10     CONTINUE
   20 CONTINUE
	
C
C                   SET TOTAL DOSE TO ZERO
C                   ================
C
          DO  30  I = 1, IOUT
              TDOSE(I) = 0.
   30     CONTINUE
	
C                   START CHAINS LOOP
C                   =================
C
          DO  80  NCHAIN = 1, NDC
	
	    NELCHN = NEL(NCHAIN)
            ISAVE(NCHAIN) = 1
C
C                   SET DOSRT, ADOS AND ADOS1 TO ZERO
C                   =================================
C
              C A L L   S E T Z E R
C
C
C                   JUMP CHAINS WITH NO SOURCE TERM
C                   ===============================
C
              IF( KSAVE(NCHAIN).EQ.0 )  GO TO 80
C

              C A L L   M G S P A C
C
              C A L L   M G T I M E
C
              C A L L   N E A R F
C
C                   NO OUTPUT FROM NEAR FIELD
C                   JUMP TO END OF CHAIN LOOP
C                   =========================
C
              IF( ISAVE(NCHAIN).EQ.0 )  GO TO 80
C
              C A L L   G T M 1
C
              C A L L   B I O S
C
C
C                   SECTION TO SUM DOSES
C                   OVER THE ISOTOPES AND PATHWAYS
C                   ==============================
C
              DO  50  J = 1, IOUT
                  DO  40  I = 1, NELCHN 
                      CDOSE(NCHAIN,J) = CDOSE(NCHAIN,J)
     &                                  + FDOSE(NCHAIN,I,J)
   40             CONTINUE
   		  UU(J,NCHAIN,NRUN)=CDOSE(NCHAIN,J)
  50         CONTINUE
C
C                   END CHAINS LOOP
C                   ===============
C
   80     CONTINUE	
C
C                   TOTAL DOSE COMPUTATION
C                   ======================
C
          DO  92  J = 1, IOUT
              DO  90  I = 1, NDC
                  TDOSE(J) = TDOSE(J) + CDOSE(I,J)
   90         CONTINUE
	      UU(J,MXNCHN+1,NRUN)=TDOSE(J)

   92     CONTINUE
C
C                    STORE DOSE TIME SERIES
C                    ======================
C
c        do i=1,iout
c          WRITE(33,*) tout(i),TDOSE(I)
c        end do
C
		
      RETURN
      END
