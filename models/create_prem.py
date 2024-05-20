import numpy as np
'''
# Load data
model = np.loadtxt(fname="./prem28.model", skiprows=3)

att = np.array([[0,      1221.5, 84.6, 1327.7],
                [1221.5, 3480.0, 0.0,  57823],
                [3480.0, 5701.0, 312,  57823],
                [5701.0, 6151.0, 143,  57823],
                [6151.0, 6291.0, 80,   57823],
                [6291.0, 6346.6, 600,  57823],
                [6346.6, 6371.0, 600,  57823]])

for row in range(len(att[:,0])):
    truth = np.array([model[:,0] >= att[row, 0]*1000]) & np.array([model[:,0] <= att[row, 1]*1000])
    model[:, -2][truth[0,:]] = att[row,2]
    model[:, -1][truth[0,:]] = att[row,3]

print("REMEMBER TO MANUALLY CHECK THE BOUNDARIES")

np.savetxt(fname="prem28_att_long.model", X=model, fmt="% 10f")'''

model = np.loadtxt(fname="./prem28_poly.model")

out = np.zeros((len(model[:,0]), 6))

out[:,0] = model[:,1]
out[:,1] = model[:,2]
out[:,2] = model[:,3]
out[:,3] = model[:,5]
out[:,4] = model[:,8]
out[:,5] = model[:,9]

np.savetxt(fname="prem28_poly_YSPEC.model", X=out, fmt="% 10f")
