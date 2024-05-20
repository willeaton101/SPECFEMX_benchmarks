# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Convolution with a gaussian STF using SPECFEM3D_GLOBE parameters (see manual, pg. 27):
#   https://geodynamics.org/cig/software/specfem3d_globe/specfem3d_globe-manual.pdf
# _______________________________________________________________________________________________________________________
import numpy as np
import matplotlib.pyplot as plt

def gauss_STF_convolve(time, data, half_duration, alpha=1.628, amp=1):
    time_mult = 2.5

    dt        = time[1] - time[0]                                        # Timestep (assumes a regularly spaced dt)

    stf_t     = np.arange(-time_mult*half_duration, time_mult*half_duration+dt, dt)  # Create time array for STF

    fact      = amp* alpha/((np.pi**0.5)*half_duration)                     # Gaussian pre-factor
    stf       = fact*np.exp(-(alpha*stf_t/half_duration)**2)             # Gaussian STF

    conv      = np.convolve(data, stf, mode="full")                      # Convolve signal

    time_conv = np.arange(len(conv))*dt - time_mult*half_duration    + time[0]    # Generate corresponding convolution time array

    output    = np.transpose(np.array([time_conv, conv]))                # Collate data into 2D array
    output    = output[int(np.floor(len(stf_t)/2)):, :]                  # Slice lower end of array
    output    = output[:-len(stf_t), :]                                  # Slice upper end of array

    return output



if __name__== "__main__":
    # Example:
    n = 9000
    x = np.linspace(35, 4000, n)
    y = np.zeros(n)

    y[4532] = 1
    y[532]  = 1
    y[8593] = 1
    y[3424] = -1


    fig, ax = plt.subplots()
    ax.plot(x,y)


    out = gauss_STF_convolve(time=x, data=y, half_duration=5, alpha=1.628, amp=1)

    ax.plot(out[:,0],out[:,1])
    plt.show()