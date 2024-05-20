# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    25th March 2022
# Notes:
#   Plots an example of the significance of including full gravity
# ______________________________________________________________________________________________________________________
import os
import numpy as np
from get_data import get_data, process_stream
from plot_funcs import plot_record, plot_station, plot_spectra
import obspy
from rotate_stream import rotate_stream_data
from metadata import PlotMetadata
import matplotlib.pyplot as plt

def norm(x):
    return x/np.amax(np.abs(x))

# ______________________________________________________________________________________________________________________
# SPECIFY METADATA AND PARAMETERS FOR PROCESSING/PLOTTING:
# Feel free to edit these - name needs to correspond to the label you have given the data when it loads
legend_colours = {
  "spfmx": "royalblue",
  "nmsyn": "r",
  "yspec": "lightseagreen",
  "LDDRK": "salmon",
  "No gravitation": "purple",
  "high-res": "salmon",
  "pure": "y"  ,
  "conv": "r",
  "spfmx_old": "purple",
  "spfmx_new": "lightgreen",
  "Cowling approximation": "orange",
  "Self-gravitating" : "royalblue"
}
meta = PlotMetadata(chls= "ZTP",                      # Channel names to load
                    out_chls = ["Z"],       # Channels for plotting
                    rotate = False,                   # I dont think this actually does anything!
                    gravity_subtract=False,
                    geoco = 1.0,                        # Ellipticity
                    fmin = 0.002,
                    fmax = 0.02,
                    ev_coords = (-13.82, -67.25),     # Geographic source coordinates
                    tmin=0,                           # For slicing time array
                    tmax = 18000,                     # For slicing time array
                    colours = legend_colours,         # As above
                    t_offset_spfmx = 8.999999999e-02,               # Time offset for SPFMX data
                    t_offset_nmsyn = 1.439962400e-02,   # Time offset for NMSYN data
                    t_offset_yspec = 1.439962400e-02,   # Time offset for YSPEC data
                    t_offset_axisem=9.6300000000e-02 ,  # Time offset for AXISEM data
                    no_stns=1)                         # You will need to edit gen_stn_list() in plot_meta_data_class.py
                                                        # Depending on what the stations are called in data file names


# ______________________________________________________________________________________________________________________
# LOAD DATA
# BEST CURRENT RESULTS:
#st_nmsyn_sg  = get_data(type="nmsyn", filepath="./SPECFEMX_benchmarks/test1/self_gravitation/nmsyn/conv", meta=meta,  label="Self-gravitating")
st_nmsyn_c  = get_data(type="nmsyn", filepath="./SPECFEMX_benchmarks/test1/cowling/nmsyn/conv", meta=meta,  label="Cowling Approximation")
st_yspec_ng  = get_data(type="yspec", filepath="./SPECFEMX_benchmarks/test1/no_gravitation/yspec/conv", meta=meta, yspec="conv_nograv_noatt.yspec", label="No Gravity")



# Process the data:
data_list = [st_yspec_ng, st_nmsyn_c]               # List of streams you want to put into some master, processed stream
pdf_name = "test1.pdf"


st = obspy.Stream()                            # This will be the master trace with all data
for stream in data_list:
    st +=  process_stream(stream, meta=meta)   # Check this as it might need editing for your simulation params.


print()


st = st.select(channel="T")

for i in st:
    if i.stats.label == "No Gravity":
        ng = i
    elif i.stats.label == "Cowling Approximation":
        c = i
    elif i.stats.label == "Self-gravitating":
        g = i

fig, ax = plt.subplots(2, sharex=True, sharey=True)

fig.set_tight_layout(True)


legend = []
for i in st:
    print(i.stats.delta)
    if i.stats.label[:2]=="No":
        colour = 'k'
    else:
        colour = 'red'
    ax[0].plot(np.arange(0, i.stats.delta*i.stats.npts, i.stats.delta),  norm(i.data), colour, linewidth=2)
    legend.append(i.stats.label)

ax[0].legend(legend)
ax[1].plot(np.arange(0, st[1].stats.delta*st[1].stats.npts, st[1].stats.delta), norm(st[1].data)-norm(st[0].data[:-2]), 'k-')
ax[1].plot(np.arange(0, st[1].stats.delta*st[1].stats.npts, st[1].stats.delta), norm(st[1].data)-norm(st[0].data[:-2]), ':', color="red")

ax[0].set_ylabel("Normalised displacement", fontsize=16)
ax[1].set_ylabel("Difference", fontsize=16)
ax[1].set_xlabel("Time [s]", fontsize=16)

ax[0].set_xlim([0, 17500])
ax[0].set_ylim([-1.05, 1.05])
plt.show()