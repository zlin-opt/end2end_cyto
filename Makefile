export PETSC_DIR=${HOME}/petsc-3.6.4
export PETSC_ARCH=arch-mumps-opt
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CLEANFILES = *.o *.so *.a

CC=mpicc

CFLAGS   += -O3 -Wall -march=native -fcx-limited-range -fno-exceptions -fPIC
INCFLAGS = -I. ${PETSC_CC_INCLUDES} 
LIBS= $(PETSC_LIB)

all: lib2d.so

maxwell.o: maxwell.c
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

solver.o: solver.c
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

array2vec.o: array2vec.c
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

interface.o: interface.c
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

lib2d.so: maxwell.o solver.o array2vec.o interface.o
	$(CC) -shared -o $@ *.o $(LIBS)


