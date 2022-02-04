# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    1st Feb 2022
# Notes:
#   Function to process Axisem data including renaming the files into a SPECFEMX-like style and then convolving
# ______________________________________________________________________________________________________________________import os
import numpy as np
from convolve_stf import gauss_STF_convolve
import os
import numpy as np

def process_raw_yspec(channels, master_dir, raw_folder="/raw/", new_folder="/conv/", measure="D"):

    # Create directory strings
    raw_dir = master_dir + raw_folder
    conv_dir = master_dir + new_folder

    # Now we need to load and convolve the data. Note here that I use a convolution script based based off the script
    # convolve_source_timefunction.f90 within SPECFEMGLOBE3D. However there is some scaling error between the two. This
    # is currently not fixed because I normalise all my traces anyway but needs solving at some point.

    # Create directory if required:
    if os.path.exists(conv_dir)==False:
        os.mkdir(conv_dir)

    for file in os.listdir(raw_dir):

        if "X" in file:
            # Load data
            data = np.loadtxt(raw_dir+file)
            time = np.loadtxt(raw_dir+"/data_time.ascii")

            # Loop through each channel
            for i in range(len(channels)):
                trace = data[:, i]                                                          # Isolate channel data

                convolved = gauss_STF_convolve(time=time, data=trace, half_duration=60)     # Convolve data

                out_fname = f"{conv_dir}/conv_{file}.{measure}{channels[i]}"                # Output file name
                np.savetxt(fname=out_fname, X=convolved, delimiter="         ")             # Save file

                print(f"Saved: {out_fname}")


if __name__ == "__main__":
    process_raw_yspec(channels="TPZ",
                  master_dir="./test2/no_grav/axisem/",
                  raw_folder="/raw/",
                  new_folder="/conv/",
                  measure="D")