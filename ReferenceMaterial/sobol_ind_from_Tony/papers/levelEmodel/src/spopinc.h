                                                                                
C>>>>>>>>>>>>>>INCLUDED COMMONS>>>>>>>>>>>>>                                    
                                                                                
      PARAMETER (NRUNSM=16384,NVARSM=50,NTPM=500,                                
     &          NTIMEM=44,NESTM=14,NVMOSM=10,                                   
     &          NCMBMX=100)                                                     
                                                                                
                                                                                
      COMMON /CONST/                                                            
     &            EPSI,ALPHA,BETA,GAMMA,YTRUN,IOPT                              
                                                                                
      COMMON /VARIAB/                                                           
     &            TIMES(NTIMEM),PTIME(NTPM),                                    
     &            XINP(NRUNSM,NVARSM+1),YVAR(NTPM),                             
     &            JTIME(NTIMEM)                                                 
                                                                                
      COMMON /HIMARR/                                                           
     &            VHIM(NCMBMX,NTIMEM),VHIMR(NCMBMX,NTIMEM),                     
     &            DHAT(NCMBMX,NTIMEM),DHATR(NCMBMX,NTIMEM),                     
     &            DHATER(NCMBMX,NTIMEM,2),                                        
     &            DHARER(NCMBMX,NTIMEM,2),                                        
     &            SCHIM(NTIMEM,3),SCHIMR(NTIMEM,3),                             
     &            NV(NCMBMX,NVARSM),KEL(NCMBMX),                                
     &            HIMCMB,NCMB,HIMS                                              
                                                                                
      COMMON /RUNCON/                                                           
     &            ISELVA(NVARSM),VALNOR,                              
     &            NRUNS,NVARS,NVARP1,NVMP1,NTP,NTIMES,NEST,NVMOST,              
     &            NOUTVA,NVAOUT,ISEL,                                           
     &            SELKEY(NVARSM),                                               
     &            SWPEA,SWSMI,SWTTS,SWSHI,SWHIM,SWPOIN,YESNOT(NESTM),           
     &            SWINTP,WORD,SWTIME,SWNORM,SWREGR                                  
                                                                                
      COMMON /ANAMES/                                                           
     &            VARNAM(NVARSM)                                                
     &           ,ESTNAM(NESTM)                                                 
     &           ,ACTNAM(NESTM)                                                 
     &           ,SCHEME(2)                                                     
                                                                                
      COMMON /VARSTA/                                                           
     &            VAS   (NVARSM,NTIMEM),                                        
     &            VASTAR(NVARSM,NTIMEM),                                        
     &            VMB(NRUNSM,NTIMEM),                                           
     &            VOB(NRUNSM,NTIMEM),                                           
     &            SAVG(NTIMEM),STDV(NTIMEM),                                    
     &            YVAMAX(NRUNSM,NTIMEM),YVAMEA(NTPM),YVASTD(NTPM),              
     &            YVAOUT(NRUNSM,NTIMEM),                                        
     &            NTSCA(NTIMEM),NVASC(NVARSM),NTISCA,NVASCA                     
                                                                                
      COMMON /WORK/                                                             
     &            WORK12(NESTM,NESTM)                                           
     &           ,WORK1(NRUNSM),WORK2(NRUNSM),IWORK1(NRUNSM)                    
     &           ,WORK4(NVARSM+1),WORK0(NRUNSM),WORK3(NRUNSM)                   
     &           ,WORK5(NRUNSM),WORK6(NRUNSM),WORK7(NRUNSM)                     
     &           ,WORK10(NVARSM),WORK11(NVARSM),WORK8(NTIMEM)                   
                                                                                
      COMMON /INTVAR/                                                           
     &            REGREV(8,NTIMEM,NVARSM)                                       
     &           ,EST(NESTM,NVARSM),YMCD(2)                                     
     &           ,ACTEST(NESTM,NVARSM)                                          
     &           ,RANKS(NRUNSM,NVARSM+1)                                        
     &           ,CORR(NVARSM+1,NVARSM+1)                                       
     &           ,NREGRU(2,NTIMEM,NVARSM),NREGVA(2,NTIMEM)                      
     &           ,XPCT,QTT,MPLOT                                                
                                                                                
      CHARACTER*20 SCHEME                                                       
      CHARACTER*10 VARNAM                                                       
      CHARACTER*6  ESTNAM,ACTNAM                                                
      CHARACTER*4  WORD                                                         
      CHARACTER*3  SELKEY                                                       
      CHARACTER*1  YESNOT,SWPOIN,SWPEA,SWSMI,SWTTS,SWSHI,SWHIM                  
     &,SWINTP,SWTIME,SWNORM,SWREGR                                         
                                                                                
C>>>>>>>>>>>>>>>INCLUDE SECTION>>>>>>>>>>>>>                                    
                                                                                
