# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Rewritten script to convolve NMSYN or SPECFEMX data
# _______________________________________________________________________________________________________________________
import os
import numpy as np
from convolve_stf import gauss_STF_convolve

master_dir = "./PEGS/yspec/anelastic_iso2c/acc/"
raw_dir = master_dir + "/raw/"



# Identifier in file name that indicates this is a file we want to convolve - e.g for synthetics may want to use 'X'
id = 'Z'       #IU
half_duration = 70 # use 60 for benchmarks tests // 70 for Tohoku

conv_dir = master_dir + f"/conv/"

if os.path.exists(conv_dir) == False:
    os.mkdir(conv_dir)

for file in os.listdir(raw_dir):
    if id in file:
        input = np.loadtxt(raw_dir + file)

        conv = gauss_STF_convolve(input[:,0], input[:,1], half_duration=half_duration, alpha=1.628)

        out_fname =  conv_dir+"conv_"+file
        np.savetxt(fname=out_fname, X=conv)
        print(f"Completed {out_fname}")