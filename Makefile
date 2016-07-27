MF= Makefile
FC = ifort 
FFLAGS = -O3 -axSSE4.1 -openmp -ipo -free -I$(MKLROOT)/include/ia32  -I$(MKLROOT)/include #-check all -fp-stack-check -traceback -gen-interfaces -warn interfaces  -fstack-protector -debug-parameters all -debug inline-debug-info -fp-model precise -ftrapuv -g #-fpe0 -vec-report -par-report 
LIBFLAGS = -O3 -axSSE4.1 -openmp -ipo -free -L$(MKLROOT)/lib/ia32 $(MKLROOT)/lib/ia32/libmkl_blas95.a $(MKLROOT)/lib/ia32/libmkl_lapack95.a -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread # -lm #not required?
VPATH = ./src/

EXE = vesemd

SRC = \
        constants.f90 \
        shared_data.f90 \
	functions_module.f90 \
	IO_module.f90 \
	geometry_module.f90 \
	SEM_module.f90 \
	waters_solution_module.f90 \
	boundary_module.f90 \
	initial_conditions_module.f90 \
	OIFS_module.f90 \
	movingSphere_module.f90 \
	pardiso_solver.f90 \
	devss_module.f90 \
	result_analysis.f90 \
	viscoelastic_module.f90 \
	main.f90 \


#
# No need to edit below here
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ= $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $< 

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) -o $@ $(OBJ) ${LIBFLAGS}

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f *.o *.bin *.mod core $(EXE)



