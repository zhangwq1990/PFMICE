
# fortran compiler
#FC = mpiifort
FC = mpif90
#FC   = ftn

#FFLAGS := -r8 -fpconstant -O3 -132 -cpp #-ipo
FFLAGS := -O3 -ffixed-line-length-none -mcmodel=large -fdefault-real-8 -cpp -Wall -fcheck=all -fbounds-check
#FFLAGS := -O3 -ffixed-line-length-none -mcmodel=large -fdefault-real-8 -cpp -fbounds-check
## FFLAGS := -r8 -O3 -mcmodel=medium #-Mlarge_arrays -mp #-fopenmp 

BIG :=  -mcmodel=medium
#DBG := -g -traceback
PROF :=# -pg
OMP := #-openmp
# 2decomp&fft
include libs/2decomp_fft/src/Makefile.inc
INCLUDE = -I libs/2decomp_fft/include
#
LIB = -L libs/2decomp_fft/lib -l2decomp_fft -L libs/fft/lib -lfft

TARGET = PFMICE

SRC = constants.f90  \
      vars.f90 \
      init.f90   \
      logo.f90   \
      field.f90  \
      parallel.f90 \
      chem.f90 \
      prefft.f90 \
      bou.f90   \
      restart.f90 \
      loadd.f90  \
      conserve.f90 \
      thomas.f90 \
      phi.f90 \
      phi_im.f90 \
      mapc.f90  \
      phi_RHS.f90 \
      revise.f90 \
      triperiodic.f90 \
      c.f90   \
      cim.f90   \
      c2.f90   \
      checknan.f90 \
      uvwstar.f90 \
      p_RHS.f90 \
      zredistribute.f90 \
      fastsolver.f90 \
      correc.f90 \
      temp.f90 \
      energy.f90 \
      vtk_write.f90 \
      debug_write.f90 \
      main.f90

OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -o $@ $(OBJ) $(LIB)

.PHONY: clean  
clean:
	rm -f *.o *.mod $(TARGET) .depend

veryclean: 
	rm -f *.o *.mod $(TARGET) .depend libs/fft/fft.o /libs/fft/lib/libfft.a; cd libs/2decomp_fft/src/; make clean; cd ../../../; rm ./data/*

libraries:
	cd libs/fft; $(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -c fft.f; ar qc libfft.a fft.o; mv libfft.a lib; cd ../../;
	cd libs/2decomp_fft/src/; make; cd ../../../

%.o: %.f90
	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) $<

.depend dep:
	./.makedepo $(SRC) > .depend

#include .depend
