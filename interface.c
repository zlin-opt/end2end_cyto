#include "interface.h"

#undef __FUNCT__
#define __FUNCT__ "initialize"
ops_ *initialize(MPI_Comm comm,
		 PetscInt Nx, PetscInt Nz,
		 PetscInt Npmlx0, PetscInt Npmlx1, PetscInt Npmlz0, PetscInt Npmlz1,
		 PetscReal dx, PetscReal dz,
		 PetscReal omega,
		 PetscInt maxit)
{

  PetscInitialize(NULL,NULL,NULL,NULL);
  
  ops_ *ops = (ops_ *)malloc(sizeof(ops_));

  ops->comm=comm;

  ops->Nx=Nx;
  ops->Nz=Nz;
  ops->Npmlz0=Npmlz0;
  ops->Npmlz1=Npmlz1;
  ops->dz=dz;
  
  create_DDe(comm,&(ops->DDe),
	     Nx,Nz,
	     Npmlx0,Npmlx1,Npmlz0,Npmlz1,
	     dx,dz,
	     omega);
  ops->omega=omega;
  
  ops->maxit=maxit;
  setupKSPDirect(comm, &(ops->ksp), maxit);
  ops->its = (PetscInt *) malloc(sizeof(PetscInt));
  *(ops->its)=100;

  VecCreateMPI(comm,PETSC_DECIDE,Nx*Nz,&(ops->Ey));
  
  return ops;
  
}

#undef __FUNCT__
#define __FUNCT__ "forward"
void forward(ops_ *ops,
	     double *eps,
	     double *Eyreal, double *Eyimag)
{

  PetscInt Nx=ops->Nx, Nz=ops->Nz, Npmlz0=ops->Npmlz0, Npmlz1=ops->Npmlz1;
  PetscReal dz=ops->dz, omega=ops->omega;
  PetscInt Nxz=Nx*Nz;
  PetscInt nrows,ncols;
  MatGetSize(ops->DDe, &nrows,&ncols);
  PetscPrintf(ops->comm, "Nx: %d, Nz: %d, Nxz: %d \n",Nx,Nz,Nxz);
  PetscPrintf(ops->comm, "The dimensions of DDe is %d x %d. \n",nrows,ncols);
  
  Mat M;
  Vec w2eps,b;
  MatDuplicate(ops->DDe,MAT_COPY_VALUES,&M);
  VecDuplicate(ops->Ey,&w2eps);
  VecDuplicate(ops->Ey,&b);

  array2mpi(eps, REAL, w2eps);
  VecScale(w2eps,-omega*omega);
  MatDiagonalSet(M,w2eps,ADD_VALUES);

  PetscScalar *Jy = (PetscScalar *)malloc(Nxz*sizeof(PetscScalar));
  for(PetscInt iz=0;iz<Nz;iz++){
    for(PetscInt ix=0;ix<Nx;ix++){

      PetscInt i = ix + Nx * iz;
      if(iz==Npmlz0+1)
	Jy[i]=1.0/dz;
      else
	Jy[i]=0.0;
      
    }
  }
  array2mpi(Jy, SCAL, b);
  VecScale(b,PETSC_i*omega);

  PetscScalar *Ey = (PetscScalar *)malloc(Nxz*sizeof(PetscScalar));
  SolveMatrixDirect(ops->comm, ops->ksp, M, b, ops->Ey, ops->its, ops->maxit);
  mpi2array(ops->Ey, Ey, SCAL, Nxz);
  if(Eyreal && Eyimag){
    for(PetscInt ix=0;ix<Nx;ix++){
      PetscInt iz = Nz - Npmlz1 - 1;
      PetscInt i = ix + Nx*iz;
      Eyreal[ix]=creal(Ey[i]), Eyimag[ix]=cimag(Ey[i]);
    }      
  }

  MatDestroy(&M);
  VecDestroy(&b);
  VecDestroy(&w2eps);
  free(Jy);
  free(Ey);

  PetscPrintf(ops->comm," Forward problem, (DDe - w^2 eps) Ey = i w Jy, has been solved.\n"); 

}


#undef __FUNCT__
#define __FUNCT__ "adjoint"
void adjoint(ops_ *ops,
	     double *greal, double *gimag,
	     double *grad)
{

  PetscInt Nx=ops->Nx, Nz=ops->Nz, Npmlz1=ops->Npmlz1;
  PetscReal omega=ops->omega;
  PetscInt Nxz=Nx*Nz;
  
  Vec b,u;
  VecDuplicate(ops->Ey,&b);
  VecDuplicate(ops->Ey,&u);
  
  PetscScalar *Jy = (PetscScalar *)malloc(Nxz*sizeof(PetscScalar));
  for(PetscInt iz=0;iz<Nz;iz++){
    for(PetscInt ix=0;ix<Nx;ix++){

      PetscInt i = ix + Nx * iz;
      if(iz==Nz-Npmlz1-1)
	Jy[i]=greal[ix] - PETSC_i * gimag[ix];
      else
	Jy[i]=0.0;
      
    }
  }
  array2mpi(Jy, SCAL, b);

  KSPSolveTranspose(ops->ksp,b,u);
  VecPointwiseMult(u,u,ops->Ey);
  VecScale(u,omega*omega);
  PetscScalar *gradcx = (PetscScalar *)malloc(Nxz*sizeof(PetscScalar));
  mpi2array(u, gradcx, SCAL, Nxz);
  
  if(grad){
    for(PetscInt iz=0;iz<Nz;iz++)
      for(PetscInt ix=0;ix<Nx;ix++)
	grad[ix+Nx*iz] = creal(gradcx[ix+Nx*iz]);
  }


  VecDestroy(&b);
  VecDestroy(&u);
  free(Jy);
  free(gradcx);

  PetscPrintf(ops->comm," Adjoint problem, (DDe - w^2 eps) u = gr - i * gi, has been solved.\n"); 

}


void petscfinalize()
{
    PetscFinalize();
}

