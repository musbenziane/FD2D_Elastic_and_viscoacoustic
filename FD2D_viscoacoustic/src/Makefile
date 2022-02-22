#
include Makefile.inc
#
#
#
BIN = ../bin

# Toolboxs: Optimization, SU, Interpolation 
#LIB_TOOLS_BOX = -L $(TOOLS_BOX_DIR)/OPTIMIZATION_PAR/lib/ -lSEISCOPE_OPTIM_PAR \
                -L $(TOOLS_BOX_DIR)/SU_IO/lib/ -lSU_IO \
                -L$(TOOLS_BOX_DIR)/INTERPOLATION/lib -lSEISCOPE_INTERP \
                -L$(TOOLS_BOX_DIR)/CARFS/lib -lSEISCOPE_CARFS -lm 

#INC_TOOLS_BOX = -I $(TOOLS_BOX_DIR)/OPTIMIZATION_PAR/COMMON/include \
                -I $(TOOLS_BOX_DIR)/SU_IO/include \
                -I$(TOOLS_BOX_DIR)/CARFS/include


INC= -I../include $(INC_TOOLS_BOX) $(INCPAR)
LIB= $(LIB_TOOLS_BOX) $(LIBPAR)

all  : mod

%.o : %.mod


.SUFFIXES: .f90 .o
.f90.o:
	$(FF) $(OPTFF) $(INC) -c $*.f90

.SUFFIXES: .f .o
.f.o:
	$(FC) $(OPTF) $(INC) -c $*.f
.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) -c $*.c


clean:
	$(RM) *.o  *~ *.mod  \#* *_genmod.f90 *.optrpt


SUBINV = ricker.o readmodelfiles.o read_acqui.o 


#DEPENDANCIES WITH MAKEFILE AND .h include files
$(SUBINV2) $(SUBINV) FD2DVAC.o:  Makefile

mod:$(SUBINV2) $(SUBINV)  FD2DVAC.o
	$(FL) -o $(BIN)/FD2DVAC$(EXT) $(OPTL) $(SUBINV) FD2DVAC.o $(LIB) 

