import matplotlib.pyplot as plt
import numpy as np

# Load MINEOS
nm = np.loadtxt("poly_model_mineos.txt")
# Load YSPEC
ys = np.loadtxt("prem28_att_long.model", skiprows=3)

fig, ax = plt.subplots(1)


ynm =   6371000 - nm[:,1]

# Density
ax.plot(nm[:,2], ynm)

ax.invert_yaxis()

ax.axhline(220*1000)
ax.axhline(410*1000)
ax.axhline(660*1000)


# Vp
ax.plot(nm[:,3], ynm)

# Vs
ax.plot(nm[:,5], ynm)
"""
# Q_mu
ax[3].plot(nm[:,1], nm[:,8])
ax[3].plot(ys[:,0], ys[:,4])

# Q_kappa
ax[4].plot(nm[:,1], nm[:,9], '-x')
ax[4].plot(ys[:,0], ys[:,5], '-o')"""

plt.show()