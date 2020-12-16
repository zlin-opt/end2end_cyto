import numpy as np
filename="adj_vs_cendiff.dat"
x=np.loadtxt(filename)
y=np.gradient(x[:,1],x[:,0])
z=np.transpose([x[:,0],x[:,2],y])
np.savetxt(filename,z)
