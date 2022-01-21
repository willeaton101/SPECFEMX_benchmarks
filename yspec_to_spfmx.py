# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Function to process YSPEC data including renaming the files into a SPECFEMX-like style and then convolving
# ______________________________________________________________________________________________________________________import os
import numpy as np
from convolve_stf import gauss_STF_convolve
import os
import numpy as np

def process_raw_yspec(channels, master_dir, rename_raw=True, raw_folder="/raw/", new_folder="/conv/", measure="D"):

    # Catch issue if channel string in wrong order
    if (channels[:3] != "ZNE") and (channels[:3] != "ZTP"):
        raise ValueError("Channels must be ZNE(G) or ZTP(G)" )

    # Create directory strings
    raw_dir = master_dir + raw_folder
    conv_dir = master_dir + new_folder

    if rename_raw==True:
        # Rename all of the raw files - this is sufficiently efficient because no data is loaded so I keep it in its own
        # loop for clarity
        for file in os.listdir(raw_dir):
            if "yspec" in file:

                result = file.find('.', -4)                 # Find index of . in yspec file name
                stn_number = int(file[result+1:])           # Get station number (to right of dot)
                network = f"X{(stn_number-1)*10 + 1}"       # Convert station number to X format

                new_name = file[:result+1] + network        # Replace stn format with new X format
                os.rename(raw_dir+file, raw_dir+new_name)   # Rename raw file


    # Now we need to load and convolve the data. Note here that I use a convolution script based based off the script
    # convolve_source_timefunction.f90 within SPECFEMGLOBE3D. However there is some scaling error between the two. This
    # is currently not fixed because I normalise all my traces anyway but needs solving at some point.

    # Create directory if required:
    if os.path.exists(conv_dir)==False:
        os.mkdir(conv_dir)

    for file in os.listdir(raw_dir):

        if "yspec" in file:
            # Load data
            data = np.loadtxt(raw_dir+file)
            time = data[:,0]

            # Loop through each channel
            for i in range(len(channels)):
                trace = data[:, i+1]                                                        # Isolate channel data

                convolved = gauss_STF_convolve(time=time, data=trace, half_duration=60)     # Convolve data

                out_fname = f"{conv_dir}/conv_{file}.{measure}{channels[i]}"                # Output file name
                np.savetxt(fname=out_fname, X=convolved, delimiter="         ")             # Save file

                print(f"Saved: {out_fname}")


if __name__ == "__main__":
    process_raw_yspec(channels="ZTPG",
                  master_dir="/Users/eaton/Documents/Princeton/SPECFEMX_work/benchmark/HM_test_data/test1/equator_test/yspec/full_modes/",
                  rename_raw=True,
                  raw_folder="/raw/",
                  new_folder="/conv/",
                  measure="D")