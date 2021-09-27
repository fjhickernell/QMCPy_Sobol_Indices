	nlayer=2
	tocc=1.0E02
	tmax=1.0E07
	ndc=2
	nel(1)=1
	nel(2)=3
	alamb(1,1)=4.410E-08
	alamb(2,1)=3.240E-07
	alamb(2,2)=4.370E-06
	alamb(2,3)=9.440E-05
C
C	INITIAL INVENTORY OF I-129    MOLS
        C0    ( 1,1)   = 1.000E+02
c	INITIAL INVENTORY OF NP-237   MOLS        
        C0    ( 2,1)   = 1.000E+03
c	INITIAL INVENTORY OF U-233    MOLS       
        C0    ( 2,2)   = 1.000E+02
c	INITIAL INVENTORY OF TH-229   MOLS        
        C0    ( 2,3)   = 1.000E+03
c	DISPERSION LENGTH             METERS         
        DISPC (1)      = 1.000E+01
c	RETENSION COEFFICIENT OF NP        
        RETC  ( 2,1,1) = 1.000E+02
c	RETENSION COEFFICIENT OF U       
        RETC  ( 2,2,1) = 1.000E+01
c	RETENSION COEFFICIENT OF TH        
        RETC  ( 2,3,1) = 1.000E+02
c	DISPERSION LENGTH             METERS       
        DISPC (2)      = 5.000E+00
c	RETENSION COEFFICIENT OF NP      
        RETC  ( 2,1,2) = 1.000E+02
c	RETENSION COEFFICIENT OF U        
        RETC  ( 2,2,2) = 1.000E+01
c       RETENSION COEFFICIENT OF TH
        RETC  ( 2,3,2) = 1.000E+02
c	MOLECULAR DIFFUSIVITY  M**2/A     
        DIFFM          = 0.000E+00
c	INGESTION DOSE FACTOR  I-129  SV/MOL        
        DSF   ( 1,1)   = 5.600E+01
c	INGESTION DOSE FACTOR  NP-237 SV/MOL       
        DSF   ( 2,1)   = 6.800E+03
c	INGESTION DOSE FACTOR  U-233  SV/MOL       
        DSF   ( 2,2)   = 5.900E+03
c	INGESTION DOSE FACTOR  TH-229 SV/MOL       
        DSF   ( 2,3)   = 1.800E+06
c	WATER INGESTION RATE          M**3/A       
	RMW            = 7.300E-01
