# Quick script that replaces the time values in a file with a reglarly spaced array
import matplotlib.pyplot as plt
import numpy as np

dt = 0.2394

fname = "dec28.yspec.X1.DN"
path  = "/Users/eaton/Documents/Princeton/SPECFEMX_work/benchmark/HM_test_data/test1/dec20/yspec/new_prem/dec28/split"

data = np.loadtxt(f"{path}/{fname}")
time = data[:,0]
trace = data[:,1]

t = np.arange(len(time))*dt
out = np.transpose(np.array([t, trace]))


fig, ax = plt.subplots()

ax.plot(time, trace, 'k-')
ax.plot(t, trace, 'r-')

plt.show()
