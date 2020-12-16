#ENV["JULIA_MPI_BINARY"]="system";
#using Pkg;
#Pkg.build("MPI", verbose=true);

import MPI;
MPI.Init();
comm=MPI.COMM_WORLD;

Nx=250;
Nz=250;
Npmlx0=50;
Npmlx1=50;
Npmlz0=50;
Npmlz1=50;
dx=0.02;
dz=0.02;
omega=2*pi;
maxit=50;
Nxz=Nx*Nz;

abstract type ops_ end;
ptrops=ccall((:initialize,"/u/zinlin/cytometer/lib/lib2d"), Ptr{ops_},
							    (MPI.MPI_Comm, Int64,Int64,Int64,Int64,Int64,Int64, Float64,Float64,Float64, Int64),
							    comm, Nx,Nz,Npmlx0,Npmlx1,Npmlz0,Npmlz1, dx,dz, omega,maxit);

eps=ones(Float64,Nxz);
Eyreal=zeros(Float64,Nx);
Eyimag=zeros(Float64,Nx);
grad=zeros(Float64,Nxz);

ip=Int((Nx+Nxz)/2);
dp=0.01;
nn=100;
eps[ip]=1.0;

for i=1:nn
    eps[ip]+=dp;
    ccall((:forward,"/u/zinlin/cytometer/lib/lib2d"),Cvoid, (Ptr{ops_}, Ptr{Float64},Ptr{Float64},Ptr{Float64}), ptrops,eps,Eyreal,Eyimag)
    obj=sum(Eyreal.^2 + Eyimag.^2)
    greal=2.0*Eyreal
    gimag=2.0*Eyimag
    ccall((:adjoint,"/u/zinlin/cytometer/lib/lib2d"),Cvoid, (Ptr{ops_}, Ptr{Float64},Ptr{Float64},Ptr{Float64}), ptrops,greal,gimag,grad)
    if MPI.Comm_rank(comm)==0
       println("obj: $(eps[ip]) $obj $(grad[ip])")
    end
end

ccall((:petscfinalize,"/u/zinlin/cytometer/lib/lib2d"),Cvoid, ())





