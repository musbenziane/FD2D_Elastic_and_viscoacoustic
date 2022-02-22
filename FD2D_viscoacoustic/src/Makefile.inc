
EXT=

#Intel compiler
CC =  mpiicc
FC =  mpiifort
FL =  mpiifort
FF =  mpiifort


#PREPOCESSING COMMANDES :
OPTPRE =   
OPTF = -g -O3 -assume byterecl -qopenmp  $(OPTPRE)
OPTC =  -g -O3
OPTL =  -g -O3 -assume byterecl -qopenmp  $(OPTPRE)


OPTFF = ${OPTF}

# variables pour les librairies

AR= ar
ARFUNCT= cruv
AREXTRA= xo
RANLIB= ranlib


LIBPAR = $(LMKL)  
INCPAR = 



