#
#   TRON main directory
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS   = $(OPTFLAGS) 
FFLAGS   = $(OPTFLAGS)

# Libraries.

TRON   = src/tron/$(LIB_NAME)
ICF    = src/icf/$(LIB_NAME)
COLOR  = src/coloring/$(LIB_NAME)
BLAS   = src/blas/$(LIB_NAME)
TPROBS = src/tprobs/$(LIB_NAME)
UTILS   = src/utils/$(LIB_NAME)

LIBS = $(TRON) $(ICF) $(COLOR) $(BLAS) $(TPROBS) $(UTILS) 

install: libs exec

libs: tron_lib icf_lib coloring_lib blas_lib tprobs_lib utils_lib

tron_lib: 
	cd src/tron; make

icf_lib: 
	cd src/icf; make

coloring_lib:
	cd src/coloring; make 

blas_lib: 
	cd src/blas; make

tprobs_lib:
	cd src/tprobs; make

utils_lib:
	cd src/utils; make

# MINPACK-2 Newton's method for large bound-constrained optimization

exec : driver.o $(LIBS)
	$(FC) $(FFLAGS) -o tron driver.o $(LIBS)

clean:
	cd src/tron;     make clean
	cd src/coloring; make clean
	cd src/icf;      make clean
	cd src/blas;     make clean
	cd src/tprobs;   make clean
	cd src/utils;    make clean
	- rm -f *.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f
