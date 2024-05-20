# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    4th Feb 2022
# Notes:
#   HNG SPFMX data comes in the form RTZ but everything else I produce uses TPZ
# ______________________________________________________________________________________________________________________


import os
from copy import copy

dir = "./correction_terms/specfemx70/FD/raw/OUTPUT_FILES/"

def rename(old, new, id):
    for file in os.listdir(dir):
        if id in file:
            # need to convert all of the trans --> phi before then doing radial --> theta
            if f"MX{old}" in file:
                new_file = copy(file)
                new_file = new_file.replace(f"MX{old}", f"MX{new}")
                os.rename(f"{dir}/{file}", f"{dir}/{new_file}")
                print(f"{file} --> {new_file}")


ID = "sem"
OLD = 'T'
NEW = 'P'
rename(OLD, NEW, ID)

OLD = 'R' #'R'
NEW = 'T' #'T'
rename(OLD, NEW, ID)



old = 'T'
new = 'P'

# Now do seperate loop for second radial --> theta to avoid conflics:
for file in os.listdir(dir):
    if id in file:
        # need to convert all of the trans --> phi before then doing radial --> theta
        if "MXR" in file:
            new_file = copy(file)
            new_file = new_file.replace("MXR", "MXT")
            os.rename(f"{dir}/{file}", f"{dir}/{new_file}")
            print(f"{file} --> {new_file}")