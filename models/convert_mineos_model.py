# Use this to delete the extra lines given out by MINEOS in the non-binary files of the model so that NMSYN can read it
# as its model file in the inputs
import numpy as np

m = np.loadtxt("poly_model_mineos.txt")
m = np.delete(m, [7,8,9], 1)

np.savetxt(fname="poly_model.txt", X=m, fmt="%d %10.1f %10.2f %10.2f %10.2f %10.2f %10.2f")
