C
      PARAMETER ( NPTH   =     1 )
      COMMON /COUTP /
     &  DOSRT(NPTH,MXNCHN,MXNELM,MXNTPT)      ,
     &  ADOS(MXNELM,MXNTPT)   , ADOS1(MXNELM,MXNTPT,MXNLYS)   ,
     &  ISAVE(MXNCHN)         , KSAVE(MXNCHN)
