* Written by R;
*  write.foreign(thedata, "c:/Users/TA00/My Documents/thedata.txt",  ;

DATA  rdata ;
INFILE  "c:/Users/TA00/My Documents/thedata.txt" 
     DSD 
     LRECL= 10 ;
INPUT
 x_s
 y_s
;
LABEL  x_s = "x.s" ;
LABEL  y_s = "y.s" ;
RUN;
