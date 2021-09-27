       program levele
          
C
C       THE LEVELE MODEL.
C       
C       for the SIMLAB software
C-----------------------------------------------------------------------
C     INPUT:   
C		 ..\LevelE.sam 
C
C	OUTPUT:
C		 ..\LevelE.out 
C                (total dose at a number of time points)
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL (A-H,P-Z)

      INCLUDE 'par.h'
      INCLUDE 'cparam.h'                                                        
      INCLUDE 'cvaria.h' 
      INCLUDE 'cfast.h'
      include 'ctime.h'                                                      
      INCLUDE 'cspop.h'                                                            

      DIMENSION VAR(MXNVAR),U(MXNRUN)

      IOUT=MXNTIM

       DATA (TOUT(j), j=1,MXNTIM)  /
     *  2.000E04,  3.000E04,  4.000E04,  5.000E04,  6.000E04,  7.000E04,  
     *  8.000E04,  9.000E04,  1.000E05,  2.000E05,  3.000E05,  4.000E05,   
     *  5.000E05,  6.000E05,  7.000E05,  8.000E05,  9.000E05,  1.000E06,
     *  2.000E06,  3.000E06,  4.000E06,  5.000E06,  6.000E06,  7.000E06,   
     *  8.000E06,  9.000E06 /

c      Open Input sample file

	 open(55,file='.\LevelE.sam')

c      Open Model output file

c       open(66,file='.\LevelEsimlab.out')
       open(77,file='.\LevelE.out')

c     read dummy variable
c	read(55,*) ndum0

C     read number of runs of the model
c	read(55,*) nrunmax

C     read number of parameters in the model
c	read(55,*) nvars

c     read dummy variable
c	read(55,*) dum0

      nrunmax=57344
	nvars=12

C     write header on the output file

c      write(66,*)	'3'
c      write(66,'(A11)')	'IODINE_DOSE'
c      write(66,'(A14)')	'NEPTUNIUM_DOSE'
c      write(66,'(A10)')	'TOTAL_DOSE'
c      write(66,'(A10)')	'time = yes'      
c	 write(66,*) nrunmax

	DO nrun=1,nrunmax                                                	 

         WRITE(*,*)     '# RUN = ',NRUN
              
c        read input variables
	   READ(55,*)  (XVAR(NVAR),NVAR=1,nvars)
	
         CALL levele_model  !input ---> xvar(nvars)
C			              !output  ---> UU matrix

c        write results to file
c         write(66,'(A4,1x,i10)') 'RUN',NRUN-1
c         write(66,*) MXNTIM
c         DO ITIME=1,IOUT
c	     write(66,*) int(TOUT(ITIME)),
c     +(UU(ITIME,II,NRUN),II=1,3)  
c         END DO

c        itime are the time points
c        3 means the total dose (Np + chain)
c        nrun is the i-th run

         write(77,'(1x,26g15.7)') (UU(ITIME,3,NRUN),ITIME=1,IOUT)

      END DO   !END LOOP

	CLOSE(55)
      CLOSE(66)      
      
      write(*,*) 'Level E Program terminated'
       
       STOP
       END
C2345678901234567890123456789012345678901234567890123456789012345678901   
