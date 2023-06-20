# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    4th Feb 2022
# Notes:
#   HNG SPFMX data comes in the form RTZ but everything else I produce uses TPZ
# ______________________________________________________________________________________________________________________


import os
from copy import copy

dir = "./Tohoku_2011_benchmark/specfemx/raw/"

id = "sem"

for file in os.listdir(dir):
    if id in file:
        # need to convert all of the trans --> phi before then doing radial --> theta
        if "MXT" in file:
            new_file = copy(file)
            new_file = new_file.replace("MXT", "MXP")
            os.rename(f"{dir}/{file}", f"{dir}/{new_file}")
            print(f"{file} --> {new_file}")

# Now do seperate loop for second radial --> theta to avoid conflics:
for file in os.listdir(dir):
    if id in file:
        # need to convert all of the trans --> phi before then doing radial --> theta
        if "MXR" in file:
            new_file = copy(file)
            new_file = new_file.replace("MXR", "MXT")
            os.rename(f"{dir}/{file}", f"{dir}/{new_file}")
            print(f"{file} --> {new_file}")