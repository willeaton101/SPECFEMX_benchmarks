# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Class to hold metadata for plotting/processing
# ______________________________________________________________________________________________________________________
import numpy as np

class PlotMetadata():

    def __init__(self, chls, out_chls, fmin, fmax, ev_coords,tmax, colours={}, tmin=0, no_stns=17, rotate=False, geoco=1.0,
                 t_offset_spfmx=0, t_offset_globe=0, t_offset_nmsyn=0, t_offset_axisem=0, t_offset_yspec=0,t_offset_qssp=0, gravity_subtract=False,
                 stn_file='Xstations', stn_list=[], network='Y5', attach_coords=True, measure='U', t_offset_real=0,
                 types=["axisem","yspec","nmsyn", "spfmx"], picks_file=None, corrections=None):
        # Simulation types:
        self.types = types

        # Data type: defaults to displacement
        self.measure = measure

        self.corrections = corrections

        if picks_file!=None:
            import pickle
            with open(picks_file, 'rb') as fp:
                self.picks = pickle.load(fp)


        # Gravity info:
        self.gravity_subtract = gravity_subtract

        # Bandwidth info
        self.set_freqlims(fmin, fmax)

        # Event info
        self.src_lat = ev_coords[0]
        self.src_lon = ev_coords[1]
        self.offset  = None
        self.geoco = geoco

        # Time info:
        self.tmin = tmin
        self.tmax = tmax
        self.time_offset = {"spfmx": t_offset_spfmx,
                            "specfemx": t_offset_spfmx,
                            "3Dglobe": t_offset_globe,
                            "nmsyn": t_offset_nmsyn,
                            "yspec": t_offset_yspec,
                            "qssp": t_offset_qssp,
                            "axisem": t_offset_axisem,
                            "real": t_offset_real}

        # Plots
        self.colours = colours
        self.plot_channels = out_chls

        # Stations:
        self.attach_coords = attach_coords
        self.network = network
        self.no_stns  = no_stns
        self.stn_file     = stn_file

        self.custom_networks = False

        if stn_list == []:
            # IF station list is empty then try and load from file or Xstations
            if self.stn_file == 'Xstations':
                # If using benchmark X stations then use this
                self.stn_list = self.gen_Xstn_list()
                self.Xstations = True
            else :
                self.stn_list = self.load_stn_list()
                self.no_stns  = len(self.stn_list)
                self.custom_networks = True
                self.process_stn_list()
                self.Xstations = False
        else:
            self.stn_list = stn_list

        # Channel data:
        self.channels = chls
        self.no_chnls = len(chls)
        self.rotate_channels = rotate
        self.chl_options = {"N":  "North (N)",
                            "E":  "East (E)",
                            "Z": r"Vertical (Z, $\mathbf{\hat{r}}$)",
                            "T": r"Theta (T, $\mathbf{\hat{\Theta}}$)",
                            "P": r"Phi (P, $\mathbf{\hat{\Phi}}$)",
                            "G":  "Gravity (G)",
                            "ZC.PGRAV":        r"Z correction $\nabla \phi^{E1}$ OUTSIDE Earth ",
                            "ZC.PGRAV_INSIDE": r"Z correction $\nabla \phi^{E1}$ INSIDE Earth  ",
                            "ZC.GRAV": r"Z correction $- s \cdot$ $\nabla g$",
                            "ZC.FA": r"Z Free air correction",
                            "PEGS":' Pegs calculation',
                            'GRAVIMETERZ': 'QSSP gravimeter',
                            "ZC.CORIO": 'Vertical coriolis correction',
                            "NC.CORIO": 'North coriolis correction',
                            "EC.CORIO": 'East coriolis correction',
                            "ZZ.STRAIN": 'Vertical Gravity Strain'
                            }


    def set_freqlims(self, fmin, fmax):
        self.fmin = fmin
        self.fmax = fmax
        self.fmin_mHz = fmin*1000
        self.fmax_mHz = fmax*1000
        self.period_max = 1/fmin
        self.period_min = 1/fmax

    def gen_Xstn_list(self):
        stn_list = []
        for j in range(self.no_stns):
            stn_list.append(f"X{str(10 *(j + 1) + 1)}")
        return stn_list

    def load_stn_list(self):
        # Generates list from file
        stn_list = []
        print(f"Generating stn list from: {self.stn_file}")
        with open(self.stn_file) as f:
            # Clean out \n in strings:
            for z in f.readlines():
                stn_list.append(z[:z.find('\n')])
        return stn_list


    def calc_offset(self, stn_lat, stn_lon, make_even=True):
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

        #if make_even ==True:
            #offset = stn_lat
        #    print("WARNING USING STATION LAT AS OFFSET")

        return offset

    def set_channels(self, channels):
        self.channels = channels
        self.no_chnls = len(self.channels)

    def set_gravsubtract(self, gs):
        self.gravity_subtract = gs


    def process_stn_list(self):

        stn_list_tmp = []
        self.network_list = {}
        self.lat_list = {}
        self.lon_list = {}

        for s in range(self.no_stns):
            outlist = self._sort_stn_string(self.stn_list[s])

            stmp = outlist[0]

            stn_list_tmp.append(stmp)

            self.network_list[stmp]   = outlist[1]
            self.lat_list[stmp]       = outlist[2]
            self.lon_list[stmp]       = outlist[3]

        self.stn_list = stn_list_tmp



    def _sort_stn_string(self, s):

        out_lst = []
        i = 10000
        while i!=-1:
            i = s.find(' ')
            if i==0:
                s=s[1:]
            else:
                out_lst.append(s[:i])
                s=s[i:]

        return out_lst