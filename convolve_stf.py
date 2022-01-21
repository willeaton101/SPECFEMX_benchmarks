# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Convolution with a gaussian STF using SPECFEM3D_GLOBE parameters (see manual, pg. 27):
#   https://geodynamics.org/cig/software/specfem3d_globe/specfem3d_globe-manual.pdf
# _______________________________________________________________________________________________________________________
import numpy as np

def gauss_STF_convolve(time, data, half_duration, alpha=1.628):
    dt        = time[1] - time[0]                                        # Timestep (assumes a regularly spaced dt)

    stf_t     = np.arange(-1.5*half_duration, 1.5*half_duration+dt, dt)  # Create time array for STF

    fact      = alpha/(((np.pi)**0.5)*half_duration)                     # Gaussian pre-factor
    stf       = fact*np.exp(-(alpha*stf_t/half_duration)**2)             # Gaussian STF

    conv      = np.convolve(data, stf, mode="full")                      # Convolve signal

    time_conv = np.arange(len(conv))*dt - 1.5*half_duration              # Generate corresponding convolution time array

    output    = np.transpose(np.array([time_conv, conv]))                # Collate data into 2D array
    output    = output[int(np.floor(len(stf_t)/2)):, :]                  # Slice lower end of array
    output    = output[:-len(stf_t), :]                                  # Slice upper end of array

    return output



