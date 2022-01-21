# Script to linearly interpolate for regular sampling
import numpy as np
import os

# User defined params:
sample_dt = 0.2394

path = "/Users/eaton/Documents/Princeton/SPECFEMX_work/benchmark/HM_test_data/test1/dec20/yspec/new_prem/dec28/split/"

for fname in os.listdir(path):
#for fname in ["conv_dec28.yspec.X171.DE"]:


    yspec = np.loadtxt(f"{path}/{fname}")
    time = yspec[:,0]
    f    = yspec[:,1]


    t_start = time[0]
    t_end = time[-1]

    # Create output time array:
    t = np.arange(t_start, t_end, sample_dt)

    val = []

    for i in range(len(t)):

        lower = np.array(np.where(time <= t[i]))
        lower = lower[0, -1]
        upper = np.array(np.where(time >= t[i]))
        upper = upper[0,0]

        if lower != upper:
            val.append( (f[upper] - f[lower])/(time[upper] - time[lower])*(t[i] - time[lower]) + f[lower])
            #print("interp")
        else:
            assert(f[lower] == f[upper])
            val.append(f[lower])
            #print("equal")


    output = np.transpose(np.array([t, val]))

    if os.path.exists(path + 'resampled/') == False:
        os.mkdir(path + 'resampled/')

    np.savetxt(fname=f"{path}/resampled/{fname}", X=output)
    print(f"Completed {fname}")