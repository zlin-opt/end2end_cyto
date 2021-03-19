#ENV["JULIA_MPI_BINARY"]="system";
#using Pkg;
#Pkg.build("MPI", verbose=true);

import MPI;
MPI.Init();
comm=MPI.COMM_WORLD;

Nx=250;
Npmlx0=50;
Npmlx1=50;
Npmlz0=50;
Npmlz1=50;
gap=100;
Mx=Nx-Npmlx0-Npmlx1; #design width
Mz=50;  #design thickness
Nz=Npmlz0+gap+Mz+gap+Npmlz1
dx=0.02;
dz=0.02;
omega=2*pi;
maxit=50;
Nxz=Nx*Nz;

abstract type ops_ end;
ptrops=ccall((:initialize,"/u/zinlin/cytometer/lib/lib2d"), Ptr{ops_},
							    (MPI.MPI_Comm, Int64,Int64,Int64,Int64,Int64,Int64, Float64,Float64,Float64, Int64),
							    comm, Nx,Nz,Npmlx0,Npmlx1,Npmlz0,Npmlz1, dx,dz, omega,maxit);

# to make beads, you can modify the epsbkg
epsbkg=ones(Float64,Nz,Nx);
epsdiff=3.0
dof=zeros(Float64,Mz,Mx)

Eyreal=zeros(Float64,Nx);
Eyimag=zeros(Float64,Nx);
grad=zeros(Float64,Nxz);

ipx=Int(Mx/2)
ipz=Int(Mz/2)
dp=0.01;
nn=100;

for i=1:nn

    dof[ipx,ipy] += dp
    eps = deepcopy(epsbkg)
    eps[Npmlz0+gap+1:Npmlz0+gap+1+Mz,Npmlx0+1:Npmlx0+1+Mx] += dof * epsdiff
    eps = vec(eps) #I assumed last-index-fastest flattening; change the index order if that's not true to match with the C code
    
    ccall((:forward,"/u/zinlin/cytometer/lib/lib2d"),Cvoid, (Ptr{ops_}, Ptr{Float64},Ptr{Float64},Ptr{Float64}), ptrops,eps,Eyreal,Eyimag)
    obj=sum(Eyreal.^2 + Eyimag.^2)
    greal=2.0*Eyreal
    gimag=2.0*Eyimag
    ccall((:adjoint,"/u/zinlin/cytometer/lib/lib2d"),Cvoid, (Ptr{ops_}, Ptr{Float64},Ptr{Float64},Ptr{Float64}), ptrops,greal,gimag,egrad)

    grad = epsdiff * reshape(egrad, (Nz,Nx))[Npmlz0+gap+1:Npmlz0+gap+1+Mz,Npmlx0+1:Npmlx0+1+Mx] #Again I assumed last-index-fastest reshaping; change the index order otherwise
    if MPI.Comm_rank(comm)==0
       println("obj: $(eps[ip]) $obj $(grad[ipx,ipy])")
    end
end

ccall((:petscfinalize,"/u/zinlin/cytometer/lib/lib2d"),Cvoid, ())





