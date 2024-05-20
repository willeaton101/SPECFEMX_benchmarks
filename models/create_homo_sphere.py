import numpy as np

m = np.loadtxt("prem28.model", skiprows=3)

# Homogenous values:
vals  = [3000, 5500, 3200] # rho, vp, vs

for i in range(3):
    m[:,i+1] = vals[i]


np.savetxt(fname="homosphere.model", X=m, fmt="% 10f")
