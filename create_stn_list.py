# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    4th Feb 2022
# Notes:
#   Function that creates station list - also part of Metadata class but here for external use
# ______________________________________________________________________________________________________________________
def gen_stn_list(no_stns):
    stn_list = []
    for j in range(no_stns):
        stn_list.append(f"X{str(10 * (j + 1) + 1)}")
    return stn_list