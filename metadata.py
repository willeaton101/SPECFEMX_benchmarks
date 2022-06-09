# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Class to hold metadata for plotting/processing
# ______________________________________________________________________________________________________________________
import numpy as np

class PlotMetadata():

    def __init__(self, chls, out_chls, rotate, geoco, fmin, fmax, ev_coords,colours, tmax, tmin=0, no_stns=17,
                 t_offset_spfmx=0, t_offset_nmsyn=0, t_offset_axisem=0, t_offset_yspec=0, gravity_subtract=False):
        # Simulation types:
        self.types = ["axisem","yspec","nmsyn", "spfmx"]

        # Gravity info:
        self.gravity_subtract = gravity_subtract

        # Bandwidth info
        self.fmin = fmin
        self.fmax = fmax
        self.fmin_mHz = fmin*1000
        self.fmax_mHz = fmax*1000
        self.period_max = 1/fmin
        self.period_min = 1/fmax

        # Event info
        self.src_lat = ev_coords[0]
        self.src_lon = ev_coords[1]
        self.offset  = None
        self.geoco = geoco

        # Time info:
        self.tmin = tmin
        self.tmax = tmax
        self.time_offset = {"spfmx": t_offset_spfmx,
                            "nmsyn": t_offset_nmsyn,
                            "yspec": t_offset_yspec,
                            "axisem": t_offset_axisem}

        # Plots
        self.colours = colours
        self.plot_channels = out_chls

        # Stations:
        self.no_stns  = no_stns
        self.stn_list = self.gen_stn_list()

        # Channel data:
        self.channels = chls
        self.no_chnls = len(chls)
        self.rotate_channels = rotate
        self.chl_options = {"N":  "North (N)",
                            "E":  "East (E)",
                            "Z": r"Vertical (Z, $\mathbf{\hat{r}}$)",
                            "T": r"Theta (T, $\mathbf{\hat{\Theta}}$)",
                            "P": r"Phi (P, $\mathbf{\hat{\Phi}}$)",
                            "G":  "Gravity (G)"
                            }


    def gen_stn_list(self):
        stn_list = []
        for j in range(self.no_stns):
            stn_list.append(f"X{str(10 *(j + 1) + 1)}")
        return stn_list


    def calc_offset(self, stn_lat, stn_lon):
        # Using formula:
        # \Delta = acos(sin(src_lon)*sin(stn_lon) + cos(src_lon)*cos(stn_lon)*con(abs(src_lat - stn_lat)))
        pi = np.pi
        theta_src = (pi / 2) - np.arctan(self.geoco * np.tan(self.src_lat * pi / 180))  # Source colatitude
        theta_stn = (pi / 2) - np.arctan(self.geoco * np.tan(stn_lat * pi / 180))  # Station colatitude

        phi_src = self.src_lon * pi / 180  # Source longitude
        phi_stn = stn_lon * pi / 180  # Station longitude

        # Calculate epicentral distance $\Theta$ (scalar):
        a = np.cos(theta_stn) * np.cos(theta_src)
        b = np.sin(theta_stn) * np.sin(theta_src) * np.cos(phi_stn - phi_src)
        offset = np.arccos(a + b)*180/pi # In degrees

        offset = stn_lat
        print("WARNING USING STATION LAT AS OFFSET")

        return offset

    def set_channels(self, channels):
        self.channels = channels

    def set_gravsubtract(self, gs):
        self.gravity_subtract = gs