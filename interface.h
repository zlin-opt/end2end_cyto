#ifndef GUARD_obj_h
#define GUARD_obj_h

#include "petsc.h"
#include "maxwell.h"
#include "solver.h"
#include "array2vec.h"

typedef struct{

  MPI_Comm comm;

  PetscInt Nx;
  PetscInt Nz;
  PetscInt Npmlz0;
  PetscInt Npmlz1;
  PetscReal dz;
  
  Mat DDe;

  PetscReal omega;
  
  KSP ksp;
  PetscInt *its;
  PetscInt maxit;

  Vec Ey;
  
} ops_;

ops_ *initialize(MPI_Comm comm,
		 PetscInt Nx, PetscInt Nz,
		 PetscInt Npmlx0, PetscInt Npmlx1, PetscInt Npmlz0, PetscInt Npmlz1,
		 PetscReal dx, PetscReal dz,
		 PetscReal omega,
		 PetscInt maxit);

void forward(ops_ *ops,
	     double *eps,
	     double *Eyreal, double *Eyimag);

void adjoint(ops_ *ops,
	     double *greal, double *gimag,
	     double *grad);

void petscfinalize();

#endif
