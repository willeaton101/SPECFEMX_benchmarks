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

master_dir = "./self_gravitation/nmsyn/"
raw_dir = master_dir + "/raw/"
conv_dir = master_dir + "/conv/"

if os.path.exists(conv_dir) == False:
    os.mkdir(conv_dir)

for file in os.listdir(raw_dir):
    if "Y5" in file:
        input = np.loadtxt(raw_dir + file)

        conv = gauss_STF_convolve(input[:,0], input[:,1], half_duration=60, alpha=1.628)

        out_fname =  conv_dir+"conv_"+file
        np.savetxt(fname=out_fname, X=conv)
        print(f"Completed {out_fname}")