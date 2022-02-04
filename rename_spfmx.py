# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    4th Feb 2022
# Notes:
#   Hom Naths SPFMX data comes in the form RTZ but everything else I produce uses TPZ
# ______________________________________________________________________________________________________________________


import os
from copy import copy

dir = "./test1/cowling/specfemx/raw/"



for file in os.listdir(dir):
    if "Y5" in file:
        # need to convert all of the trans --> phi before then doing radial --> theta
        if "T" in file:
            new_file = copy(file)
            new_file = new_file.replace("T", "P")
            os.rename(f"{dir}/{file}", f"{dir}/{new_file}")
            print(f"{file} --> {new_file}")

# Now do seperate loop for second radial --> theta to avoid conflics:

for file in os.listdir(dir):
    if "Y5" in file:
        # need to convert all of the trans --> phi before then doing radial --> theta
        if "R" in file:
            new_file = copy(file)
            new_file = new_file.replace("R", "T")
            os.rename(f"{dir}/{file}", f"{dir}/{new_file}")
            print(f"{file} --> {new_file}")