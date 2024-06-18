ARCH           = LINUXAMD64
ARCHDIR        = lib_$(ARCH)/molfile
FC             = ifort
#FC             = g77 #-fno-second-underscore
#CC             = gcc
CC             = icc
LD             = $(FC)
#OPT            = -O2
CPPFLAGS       = -Iplugins/include -I$(ARCHDIR) -D_F77_FD2_UNDERSCORE
CFLAGS         = -c $(CPPFLAGS) $(OPT)
#FFLAGS         = -c $(OPT) -check all -traceback -openmp -fp-stack-check
#NG     FFLAGS         = -c $(OPT) -check all -traceback -openmp-stubs -fp-stack-check
FFLAGS         = -c $(OPT) -check all -traceback -fp-stack-check
#NETCDFLIB      = -L/usr/local/amber10-intel/src/netcdf/lib
#NETCDFLDFLAGS  = -lnetcdf
#TCLLIB         = -L/usr/lib64/
#TCLLDFLAGS     = -ltcl8.5
LDFLAGS        = -L$(ARCHDIR) -liomp5
LDLIBS         = -lmolfile_plugin -lstdc++
#CODEOPTFLAGS= -O2
LARGEARRAYFLAGS= -mcmodel=large -i-dynamic


OBJECTS = decl_solv.o water_solv.o inp_out_solv.o main_solv.o f77_molfile.o 

default: Solvation

Solvation: $(OBJECTS)
	$(FC) -o $@ $(LDFLAGS) $^ $(LDLIBS)
	
clean:
	@rm -rf *.o 	


# pattern rules
.SUFFIXES:
.SUFFIXES: .c .f90 .f .o

.c.o:
	$(CC) $(CFLAGS) $< -o $@

.f90.o:
	$(FC) $(FFLAGS) $< -o $@

.f.o:
	$(FC) $(FFLAGS) $< -o $@



